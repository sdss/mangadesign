#!/usr/bin/env python
# encoding: utf-8
"""
PlateMags.py

Created by José Sánchez-Gallego on 6 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    9 Feb 2014 J. Sánchez-Gallego
      Initial version
    26 Feb 2014 J. Sánchez-Gallego
      SDSS support: if no NSA images are found, it tries
      to download SDSS images.
    3 March 2014: J. Sánchez-Gallego
      Added documentation.
    25 August 2014: J. Sánchez-Gallego
      General update. Now the input is a plateInput file or a design ID.

"""

import numpy as np
import tempfile
import urllib2
import urllib
import os
import cStringIO
import bz2
import warnings
import shutil as sh
from urlparse import urlparse

from .exceptions import GohanWarning, GohanError
from . import log, config, readPath
from .utils import pywcsgrid2 as pw2
from .utils import getSDSSRun

from sdss.utilities.yanny import yanny, write_ndarray_to_yanny
from astropy import wcs
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
from scipy.misc import imsave
from scipy import interpolate
from astropy.modeling import models, fitting
from astropy import table
from matplotlib import pyplot as plt
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from matplotlib.patches import Ellipse
from matplotlib.collections import EllipseCollection
import gzip
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
matplotlib.rc('font', size=10)
matplotlib.rc('text', usetex=True)
plt.ioff()


MANGA_SAMPLE = table.Table.read(
    os.path.expandvars(config['files']['sciSample']), format='fits')

MANGACORE = readPath(config['mangacore'])

FILTERS = config['plateMags']['filters']
NFILTERS = len(FILTERS)
NSAFILTERPOS = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}

TARGET_SEEING = config['plateMags']['targetSeeing']

SIMBMAP = table.Table(yanny(readPath(config['plateMags']['simbmap']),
                            np=True)['SIMBMAP'])

BOSS_SN_DATA = table.Table.read(
    readPath(config['plateMags']['BOSS_SN']),
    format='ascii.commented_header')

BOSS_SN = interpolate.interp1d(BOSS_SN_DATA['fiber2flux'][::-1],
                               BOSS_SN_DATA['sn_hour'][::-1])

NIFUS = np.sum(config['IFUs'].values())

__ALL__ = ['PlateMags', 'PlateMagsIFU']

BASE_URL = os.path.join(config['nsaImaging']['baseURL'], '{0}/')
BASE_URI = '{uri.scheme}://{uri.netloc}/'.format(uri=urlparse(BASE_URL))

password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_mgr.add_password(realm=None, uri=BASE_URI,
                          user=config['nsaImaging']['username'],
                          passwd=config['nsaImaging']['password'])
auth_handler = urllib2.HTTPBasicAuthHandler(password_mgr)
opener = urllib2.build_opener(auth_handler)
urllib2.install_opener(opener)

PLATE_MAGS_TABLE = table.Table(
    None,
    names=['IFUDESIGN', 'MANGAID', 'FNUMDESIGN',
           'RAOFF', 'DECOFF', 'FIBER2MAG'],
    dtype=[int, 'S20', int, 'f8', 'f8', ('f8', NFILTERS)])

TMP_DIR = tempfile.gettempdir()

SDSS_FRAME = 'http://data.sdss3.org/sas/dr10/boss/photoObj/frames/' + \
    '{rerun:d}/{run:d}/{camcol:d}/frame-{filter}-{run:06d}-{camcol:d}-' + \
    '{field:04d}.fits.bz2'
SDSS_PSFIELD = 'http://data.sdss3.org/sas/dr10/boss/photo/redux/' + \
    '{rerun:d}/{run:d}/objcs/{camcol:d}/psField-{run:06d}-{camcol:d}-' + \
    '{field:04d}.fit'
SDSS_PHOTOFIELD = 'http://data.sdss3.org/sas/dr10/boss/photoObj/' + \
    '{rerun:d}/{run:d}/photoField-{run:06d}-{camcol:d}.fits'
SDSS_FRAME_IRG = 'http://data.sdss3.org/sas/dr10/boss/photoObj/frames/' + \
    '{rerun:d}/{run:d}/{camcol:d}/frame-irg-{run:06d}-{camcol:d}-' + \
    '{field:04d}.jpg'

SCALE = 0.396
IFU_MIN_SIZE = 19
NPIX = 50

DPI = 72


class PlateMags(list):
    """Creates a list of PlateMagsIFU.

    This class is initialised using a PlateInput instance and returns
    a list of `PlateMagsIFU` objects, one for each IFU target in the
    PlateInput.

    Parameters
    ----------
    plateInput : str
        The path to the plateInput file to use.
    kwarg : dict
        Keyword arguments to be passed to PlateMagsIFU.

    """

    def __init__(self, plateInput, designid=None, **kwargs):

        self.plateInput = plateInput
        self.struct1 = table.Table(
            yanny(self.plateInput, np=True)[config['mangaInputStructure']])
        self.kwargs = kwargs

        self.plateMags = None

        if designid is not None:
            self.designid = designid
        else:
            self.designid = int(self.plateInput.split('_')[-1][0:4])

        plateMagsIFUList = [PlateMagsIFU(row, **self.kwargs)
                            for row in self.struct1
                            if row['ifudesignsize'] >= IFU_MIN_SIZE]

        list.__init__(self, plateMagsIFUList)

        self._createPlateMags()

    @property
    def preImages(self):
        """Returns a list of preimage filenames for all IFUs."""
        return [ss.preImage for ss in self if ss.preImage is not None]

    def plot(self, filename=None, toRepo=False):
        """Plots the bundle position and flux for each target.

        This method recursively calls `PlateMags.plot()` for each target
        and saves the result in a single, multipage PDF file with name
        ``filename``. If ``filename=None``, the default name is
        ``plateMags-XXXX.pdf``` with ``XXXX`` the design ID.

        Refer to `PlateMags.plot()` for details on the layout of the plot.
        ``copyToRepo/moveToRepo`` will copy/move the resulting plot to the
        corresponding ``$MANGACORE/plateMags`` path.

        """

        if filename is None:
            filename = 'plateMags-{0:04d}.pdf'.format(self.designid)

        plt.cla()
        plt.clf()

        pp = PdfPages(filename)

        for ifu in self:

            if ifu.fluxTable is None:
                continue

            ifu.plot(None, multipage=pp)

        pp.close()

        log.info('plot saved as {0}'.format(filename))

        if toRepo:
            self._toRepo(filename)

    def _createPlateMags(self):
        """Creates a `Table` with the plateMags columns."""

        plateMags = PLATE_MAGS_TABLE.copy()
        missingIFUs = 0

        for ifu in self:
            if ifu.fluxTable is None:
                missingIFUs += 1
                continue
            ifuDesign = ifu.row['ifudesign']
            mangaID = ifu.row['mangaid']
            for fibre in ifu.fluxTable:
                fNumDesign = fibre['fnumdesign']
                raOff = fibre['raoff']
                decOff = fibre['decoff']
                fluxes = fibre['fluxes']
                fluxesNaN = fluxes.copy()
                fluxesNaN[fluxesNaN < 0.] = np.nan
                mags = 22.5 - 2.5 * np.log10(fluxesNaN)

                plateMags.add_row(
                    (ifuDesign, mangaID,
                     fNumDesign, raOff, decOff, mags))

        self.plateMags = plateMags
        self._missingIFUs = missingIFUs
        del plateMags

    def write(self, filename=None, toRepo=False):
        """Writes a plateMags file to disk.

        Parameters
        ----------
        filename : str or None, optional
            The filename of the plateMags file. If None, the default is
            ``plateMags-XXXX.par``` with ``XXXX`` the design ID.
        copyToRepo/moveToRepo : bool, optional
            If True, copies/moves the plateMags file to the appropiate
            directory in $MANGACORE/plateMags.

        """

        if len(self.plateMags) == 0:
            raise ValueError('plateMags has not values.')

        if filename is None:
            plateMagsFilename = 'plateMags-{0:04d}.par'.format(
                self.designid)
        else:
            plateMagsFilename = filename

        if os.path.exists(plateMagsFilename):
            os.remove(plateMagsFilename)

        write_ndarray_to_yanny(plateMagsFilename, self.plateMags._data,
                               structname='PLATEMAGS')

        if self._missingIFUs > 0:
            warnings.warn('no data for {0} IFUs. '.format(
                self._missingIFUs) + 'The remaining data has been saved.',
                GohanWarning)

        log.info('plateMags saved as {0}'.format(plateMagsFilename))

        if toRepo:
            self._toRepo(plateMagsFilename)

        return plateMagsFilename

    def _toRepo(self, filename):
        """Copies/moves a file to $MANGACORE/plateMags/ """

        plateMagsPath = os.path.join(MANGACORE, 'platedesign/platemags/')
        dDir = 'D00{0:s}XX/'.format('{0:d}'.format(self.designid)[0:2])
        plateMagsPath += dDir

        if not os.path.exists(plateMagsPath):
            os.makedirs(plateMagsPath)

        sh.copy(filename, os.path.join(plateMagsPath, filename))

        log.info('{0} copied to repo.'.format(filename))


class PlateMagsIFU(object):
    """A class that computes the flux for each IFU on a target.

    This class downloads imaging for a target from NSA or SDSS,
    creates preimaging files, bins the data to match a single fibre
    spatial resolution and calculates the flux for each fibre in an
    IFU bundle.

    The class uses the information from the corresponding row in
    plateInput to determine the target and from where the data must be
    downloaded. If IAUNAME, SUBDIR and PID are present and the imaging
    can be found in the NSA atlas, those images are preferred. Otherwise,
    the RA and DEC coordinates of the target are used to identify the
    SDSS field and that imaging is used.

    Preimaging FITS files are created by trimming and concatenating the
    imaging data. The format of those multiextension FITS files is

    #0: empty extension
    #1: imaging for the first filter
    #2: inverse variance for the first filter
    #3: PSF for the imaging data of the first filter
    #4: imaging for the second filter
    ...

    The filters used are defined in `Gohan.defaults.FILTERS` and must
    always be a string containing a substring of 'ugriz'.

    The data is then rebinned to match the spatial resolution of one
    fibre (defined in `TARGET_SEEING`) and fluxes for each filter are
    measured at the position of each fibre. Note that in any case a 127-fibre
    bundle is assumed, regardless of the actual IFUDESGNSIZE.

    The class is not intended to be used directly but to be called
    from `PlateMags`.

    Parameters
    ----------
    struct1Row : `astropy.table.Row`
        A row of a `Table` containing the plateInput information for
        a single target.
    keepImages : bool, optional
        If False (the default), intermediate step imaging is deleted
        when not needed anymore.
    createPreimages : bool, optional
        If True (the default), MaNGA preimaging files are saved to disk.
        These preimages are not deleted even if ``keepImages=True``.
    overwrite : bool, optional
        If True, imaging data is redownloaded and rewritten even if imaging
        with matching filenames is found. Otherwise, the already existing data
        is read.

    """

    def __init__(self, struct1Row, keepImages=False,
                 overwrite=False, createPreimages=True, **kwargs):

        self.row = struct1Row
        self.RA = self.row['ra']
        self.Dec = self.row['dec']
        self.MaNGAID = self.row['mangaid'].strip()

        self.createPreimages = createPreimages
        self.overwrite = overwrite

        if keepImages or self.createPreimages:
            self._tmpDir = '.'
        else:
            self._tmpDir = TMP_DIR

        # Initialises some variables
        self.fluxTable = None
        self.dataOrigin = None
        self._imageName = None
        self.preImage = self._tmpDir + \
            '/preimage-{0}.fits.gz'.format(self.MaNGAID)
        self.preimageIRG = self._tmpDir + \
            '/preimage-{0}_irg.jpg'.format(self.MaNGAID)

        log.debug('MANGAID: {0}'.format(self.MaNGAID))

        preimageData = self._getPreimageData()

        if preimageData is not None:

            self.binnedData = self.binData(preimageData)
            preimageData.close()  # Not needed anymore

            self._calculateFluxes(self.binnedData)

        else:

            log.important(
                'no available data for {0}'.format(self.MaNGAID.strip()))

        if not keepImages:
            self._deleteImages()

    def _deleteImages(self):
        """Deletes unnecessary files."""

        self.removeFile(self._imageName)
        self._imageName = None
        if not self.createPreimages:
            self.removeFile(self.preImage)
            self.removeFile(self.preimageIRG)
            self.preImage = None
            self.preimageIRG = None

    def removeFile(self, file):
        if file is not None and os.path.exists(file):
            os.remove(file)

    def _getPreimageData(self):
        """Creates the preimage file and downloads IRG imaging."""

        preimageData = None

        if os.path.exists(self.preImage) and not self.overwrite:

            preimageData = fits.open(self.preImage)
            log.debug('using previously saved preimage data.')
            if 'DSOURCE' in preimageData[1].header:
                self.dataOrigin = preimageData[1].header['DSOURCE']
            return preimageData

        else:

            nsaSuccess = self._getFromNSA()

            if nsaSuccess is False:
                sdssSuccess = self._getFromSDSS()

                if sdssSuccess is False:
                    return None

        data = fits.open(self._imageName)
        preimageData = self._cropHDUs(data)

        if self.createPreimages:

            tmpNonGzip = tempfile.NamedTemporaryFile(
                dir=self._tmpDir, suffix='.fits', delete=False)
            preimageData.writeto(tmpNonGzip)
            preimageData.close()

            preimageUnit = gzip.open(self.preImage, 'wb')
            tmpData = open(tmpNonGzip.name, 'rb')
            preimageUnit.writelines(tmpData)
            tmpData.close()
            preimageUnit.close()

            os.remove(tmpNonGzip.name)

        # It has to be r because the astrometry of the IRG image
        # in SDSS matches the one in the r frame. In NSA all frames
        # are aligned
        rFilt = FILTERS.index('r')
        rIdx = 3*rFilt+1
        self._getIRGImage(data[rIdx].header)

        data.close()

        return preimageData

    def _getNSAParams(self, data):
        """Determines if enough data is present to download NSA imaging."""

        cols = [col for col in data.colnames]
        neededCols = ['iauname', 'subdir', 'pid']

        returnValues = []
        for neededCol in neededCols:

            if neededCol not in cols:
                return None

            if data[neededCol] in [-999, '-999', '-999.', 'NULL']:
                return None
            returnValues.append(data[neededCol])

        return returnValues

    def _getFromNSA(self):
        """Main routine to download NSA data for the target."""

        log.debug('... Trying to get NSA data.')

        nsaParams = self._getNSAParams(self.row)

        if nsaParams is not None:

            self.IAUName = nsaParams[0]
            self.subDir = nsaParams[1]
            self.pID = nsaParams[2]

        else:

            log.debug('{0} not found in NSA.'.format(self.MaNGAID.strip()))
            return False

        self._imageName = self._tmpDir + '/{0}_NSA.fits'.format(self.IAUName)

        if not os.path.exists(self._imageName) or self.overwrite:

            nsaImageSuccess = self.getNSAImage(
                self.IAUName, self.subDir, self.pID)

            if nsaImageSuccess is False:
                log.debug('{0} has no NSA data.'.format(self.MaNGAID))
                return False

        else:

            log.debug('{0} already exists.'.format(
                os.path.basename(self._imageName)))

        self.dataOrigin = 'NSA'

        return True

    def getNSAImage(self, IAUName, subDir, pID):
        """Downloads NSA imaging and creates the parent FITS file."""

        baseURL = BASE_URL.format(subDir)

        parentURL = baseURL + \
            '/atlases/{0}/{1}-parent-{0}.fits.gz'.format(pID, IAUName)
        varURL = baseURL + \
            '/atlases/{0}/{1}-ivar-{0}.fits.gz'.format(pID, IAUName)
        bpsfURL = [baseURL + '/{0}-{1}-bpsf.fits.gz'.format(IAUName, band)
                   for band in FILTERS]

        hdus = []
        for url in [parentURL] + [varURL] + bpsfURL:
            hdus.append(self._getNSAHDU(url))
            if None in hdus:
                return False

        fullHDU = self.concatenateNSAHDUs(hdus)

        fullHDU.writeto(self._imageName, clobber=True)
        fullHDU.close()

        for ii in hdus:
            ii.close()
        del hdus

        return True

    def _getNSAHDU(self, url):
        """Gets a NSA file."""

        ff = tempfile.NamedTemporaryFile(
            dir=self._tmpDir, suffix='.fits.gz', delete=False)

        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError:
            log.debug('URL not found.')
            self.removeFile(ff.name)
            return None

        ff.write(response.read())
        ff.close()

        hduTmp = fits.open(ff.name)

        self.removeFile(ff.name)

        return hduTmp

    @staticmethod
    def concatenateNSAHDUs(hdus):
        """Concatenates the various NSA data in a single HDUList."""

        # Concatenates the images, vars and PSFs in a single FITS.
        fullHDU = fits.HDUList([fits.PrimaryHDU()])
        for ii, filt in enumerate(FILTERS):
            nn = NSAFILTERPOS[filt]
            fullHDU.append(hdus[0][nn].copy())    # image
            fullHDU.append(hdus[1][nn].copy())    # variance
            fullHDU.append(hdus[2+ii][0].copy())  # psf

        for ii in range(NFILTERS):
            fullHDU[ii*3+1].header['EXTNAME'] = '{0} img'.format(FILTERS[ii])
            fullHDU[ii*3+1].header['FLUXUNIT'] = 'nanomaggies'
            fullHDU[ii*3+2].header['EXTNAME'] = \
                '{0} ivar'.format(FILTERS[ii])
            fullHDU[ii*3+2].header['FLUXUNIT'] = 'nanomaggies'
            fullHDU[ii*3+3].header['EXTNAME'] = \
                '{0} psf'.format(FILTERS[ii])

        fullHDU.update_extend()

        return fullHDU

    def _cropHDUs(self, hdu):
        """Trims and HDU to the size of an IFU."""

        log.debug('... Creating preimage.')

        newHDU = self.copyHDUList(hdu)

        ra = self.row['ra']
        dec = self.row['dec']

        for ii in range(NFILTERS):

            imgIdx = ii*3 + 1
            varIdx = ii*3 + 2
            psfIdx = ii*3 + 3

            ww = wcs.WCS(newHDU[imgIdx].header)
            xx, yy = ww.wcs_world2pix(ra, dec, 0)

            xmin = int(xx) - NPIX
            xmax = int(xx) + NPIX
            ymin = int(yy) - NPIX
            ymax = int(yy) + NPIX

            for nn in [imgIdx, varIdx, psfIdx]:
                newHDU[nn].header.set(
                    'DSOURCE', self.dataOrigin, 'The source of the data')

            for nn in [imgIdx, varIdx]:
                newHDU[nn].data = np.array(
                    newHDU[nn].data[ymin:ymax, xmin:xmax])
                newHDU[nn].header['CRPIX1'] -= xmin
                newHDU[nn].header['CRPIX2'] -= ymin

        return newHDU

    def _getFromSDSS(self):
        """Determines the SDSS field for the target and gets the data."""

        log.debug('... Trying to get SDSS data.')

        self.SDSSField = getSDSSRun(self.RA, self.Dec)
        run, rerun, camcol, field = self.SDSSField

        if None in [run, rerun, camcol, field]:
            log.debug('galaxy not in SDSS footprint.')
            return False

        self._imageName = self._tmpDir + '/Field_{0}_{1}_{2}_SDSS.fits'.format(
            run, camcol, field)

        if os.path.exists(self._imageName) and not self.overwrite:
            log.debug('{0} already exists.'.format(
                os.path.basename(self._imageName)))

        else:

            try:
                self._getSDSSImage(run, rerun, camcol, field)
            except:
                raise GohanError('images could not be downloaded.')

        self.dataOrigin = 'SDSS'

        return True

    def _getSDSSImage(self, run, rerun, camcol, field):
        """Download the SDSS images.

        SDSS frames are downloaded for each filter, as well as psField files
        used for PSF reconstruction. The final image is a multiextension FITS
        file with the same format as the preimage.

        """

        fieldDic = {'run': run, 'rerun': rerun,
                    'camcol': camcol, 'field': field}

        frameFilenames = []
        for filt in FILTERS:
            fieldDic.update({'filter': filt})
            frameFilename = self._getTmpFilename()
            urllib.urlretrieve(SDSS_FRAME.format(**fieldDic), frameFilename)
            frameFilenames.append(frameFilename)

        psFieldFilename = self._getTmpFilename(suffix='.fit')
        urllib.urlretrieve(SDSS_PSFIELD.format(**fieldDic), psFieldFilename)

        photoFieldFilename = self._getTmpFilename(suffix='.fits')
        urllib.urlretrieve(
            SDSS_PHOTOFIELD.format(**fieldDic), photoFieldFilename)

        image = fits.HDUList([fits.PrimaryHDU()])
        for frameFN in frameFilenames:

            frame, iVar = self._readSDSSFrame(frameFN, photoFieldFilename)

            filt = frame.header['FILTER'].lower()

            ww = wcs.WCS(frame.header)
            jj, ii = ww.wcs_world2pix(self.RA, self.Dec, 0)
            psf = self._reconstructPSF(psFieldFilename, filt, ii, jj)

            image.append(frame)
            image.append(iVar)
            image.append(psf)

        image.writeto(self._imageName, clobber=True)

        # Removes temporary files.
        for img in frameFilenames + [photoFieldFilename, psFieldFilename]:
            self.removeFile(img)

        image.close()

        return True

    def _getTmpFilename(self, suffix='.fits.bz2'):
        """Creates a temporary file and returns its filename."""

        tmpFile = tempfile.NamedTemporaryFile(
            dir=self._tmpDir, suffix='.fits.bz2', delete=True)
        filename = tmpFile.name
        tmpFile.close()
        return filename

    def _readSDSSFrame(self, filename, photoFieldFilename):
        """Reads the photoField file and return imaging and inverse variance.

        This method opens the photoField frame and returns two HDUs, one
        with the imaging and the other with inverse variance data.

        """

        hdu = fits.open(bz2.BZ2File(filename))
        frame = fits.ImageHDU(data=hdu[0].data, header=hdu[0].header)
        calib = np.tile(hdu[1].data, (frame.data.shape[0], 1))
        sky = self.getSky(hdu[2].data)
        hdu.close()

        filt = frame.header['FILTER'].lower()
        filterIdx = 'ugriz'.index(filt)

        frame.header['EXTNAME'] = '{0} img'.format(filt)
        frame.header['FLUXUNIT'] = 'nanomaggies'

        photoField = fits.open(photoFieldFilename)
        gain = photoField[1].data['gain'][0][filterIdx]
        darkVariance = photoField[1].data['dark_variance'][0][filterIdx]
        photoField.close()

        dn = frame.data / calib + sky
        # nElec = dn * gain
        dnErr = np.sqrt(dn / gain + darkVariance)
        imgErr = dnErr * calib

        iVar = (1. / imgErr)**2
        iVarHDU = fits.ImageHDU(data=iVar, header=frame.header)
        iVarHDU.header['EXTNAME'] = '{0} ivar'.format(filt)
        iVarHDU.header['FLUXUNIT'] = 'nanomaggies'

        return frame, iVarHDU

    @staticmethod
    def getSky(tt):
        """Returns an sky image from input SDSS sky data."""

        xInterp = tt['XINTERP'][0]
        yInterp = tt['YINTERP'][0]

        zz = tt['ALLSKY'][0]
        xx = np.arange(0, zz.shape[0])
        yy = np.arange(0, zz.shape[1])
        ff = interpolate.interp2d(xx, yy, zz, kind='cubic')

        return ff(xInterp, yInterp)

    @staticmethod
    def _reconstructPSF(psFieldFilename, filt, row, col):
        """Reconstruct the PSF image at a certain position.

        Parameters
        ----------
        psFieldFilename : str
            The path of the psField file.
        filt : str
            The filter for which the PSF is being reconstructed.
            Must be one of: ``'u'``, ``'g'``, ``'r'``, ``'i'`` or ``'z'``
        row, col : int
            The row and column of the original image for which the PSF
            is being reconstructed.

        """

        filtIdx = 'ugriz'.index(filt) + 1
        psField = fits.open(psFieldFilename)
        pStruct = psField[filtIdx].data

        nrow_b = pStruct['nrow_b'][0]
        ncol_b = pStruct['ncol_b'][0]

        rnrow = pStruct['rnrow'][0]
        rncol = pStruct['rncol'][0]

        nb = nrow_b * ncol_b
        coeffs = np.zeros(nb.size, float)
        ecoeff = np.zeros(3, float)
        cmat = pStruct['c']

        rcs = 0.001
        for ii in range(0, nb.size):
            coeffs[ii] = \
                (row * rcs)**(ii % nrow_b) * (col * rcs)**(ii / nrow_b)

        for jj in range(0, 3):
            for ii in range(0, nb.size):
                ecoeff[jj] = ecoeff[jj] + \
                    cmat[ii / nrow_b, ii % nrow_b, jj] * coeffs[ii]

        psf = pStruct['rrows'][0] * ecoeff[0] + \
            pStruct['rrows'][1] * ecoeff[1] + \
            pStruct['rrows'][2] * ecoeff[2]

        psf = np.reshape(psf, (rnrow, rncol))
        # psf = psf[10:40, 10:40]  # Trim non-zero regions.

        psfHDU = fits.ImageHDU(data=psf)
        psfHDU.header['EXTNAME'] = '{0} psf'.format(filt)

        return psfHDU

    def binData(self, data):
        """Bins the data to TARGET_SEEING using a Gaussian kernel."""

        log.debug('... Rebinning images.')

        rebin = self.copyHDUList(data)

        for ii in range(NFILTERS):

            img = data[3*ii+1]
            psf = data[3*ii+3]

            # ww = wcs.WCS(img.header)
            # scale = np.abs(ww.wcs.cd[0, 0]) * 3600
            scale = SCALE
            seeingPix = self.getSeeing(psf)
            seeingArcSec = scale * seeingPix

            if seeingArcSec < TARGET_SEEING:

                log.debug(
                    '...... Band ' +
                    '{0} from seeing {1:.3f} arcsec to {2} arcsec.'.format(
                        FILTERS[ii], seeingArcSec, TARGET_SEEING))

                targetSigma = TARGET_SEEING / (2*np.sqrt(2*np.log(2))) / scale
                sourceSigma = seeingPix / (2*np.sqrt(2*np.log(2)))

                rebinnedData = self.rebinImage(img.data, sourceSigma,
                                               targetSigma)

            else:

                log.debug('...... Band {0}: '.format(FILTERS[ii]) +
                          'target seeing is lower than the image seeing.')

                rebinnedData = img.data

            rebin[3*ii+1].data = rebinnedData

        return rebin

    @staticmethod
    def copyHDUList(hdu):
        """Returns a new HDUList that is a copy of the input one."""
        newHDUs = fits.HDUList([])
        for hh in hdu:
            newHDUs.append(hh.copy())
        newHDUs.update_extend()
        return newHDUs

    @staticmethod
    def getSeeing(psf):
        """Calculates the seeing of an image from its PSF."""
        # Returns seeing in pixels

        shape = psf.data.shape

        xx, yy = np.mgrid[:shape[0], :shape[1]]

        xMean = shape[0] / 2.
        yMean = shape[1] / 2.

        gInit = models.Gaussian2D(amplitude=1., x_mean=xMean, y_mean=yMean,
                                  x_stddev=1., y_stddev=1.)
        f2 = fitting.LevMarLSQFitter()
        gg = f2(gInit, xx, yy, psf.data)

        fwhm = np.sqrt(gg.x_stddev.value * gg.y_stddev.value) * 2 * \
            np.sqrt(2 * np.log(2))

        return fwhm

    @staticmethod
    def rebinImage(data, sigma1, sigma2):

        kernel = np.sqrt(sigma2**2 - sigma1**2)
        rebinnedData = gaussian_filter(data, kernel)
        return rebinnedData

    def _calculateFluxes(self, data):
        """Calculates the flux for each fibre and updates self.fluxTable."""

        log.debug('... Calculating fluxes.')

        factor = 10  # Factor to resample the image to improve photometry.

        fluxTable = table.Table(
            None, names=['fnumdesign', 'raoff', 'decoff',
                         'x', 'y', 'fluxes'],
            dtype=[int, float, float,
                   ('float', len(FILTERS)), ('float', len(FILTERS)),
                   ('float', len(FILTERS))])

        fibreWorld = self._getFibreCoordinates()

        # Pixel scale of the resampled image.
        scaleZoomed = SCALE / factor

        # Radius of one fibre in pixels in the resampled image.
        radius = TARGET_SEEING / 2. / scaleZoomed

        fluxArray = np.zeros((len(fibreWorld), len(FILTERS)), float)
        fibrePixXArray = np.zeros((len(fibreWorld), len(FILTERS)), float)
        fibrePixYArray = np.zeros((len(fibreWorld), len(FILTERS)), float)

        for jj in range(NFILTERS):

            img = data[jj*3+1]
            ww = wcs.WCS(img.header)
            fibrePix = ww.wcs_world2pix(fibreWorld, 0)

            zoomed = self.zoomImage(img, factor=factor)
            wwZoomed = wcs.WCS(zoomed.header)
            fibrePixZoomed = wwZoomed.wcs_world2pix(fibreWorld, 0)

            fluxData = self.phot(zoomed.data, fibrePixZoomed, radius)

            for ii, ff in enumerate(fluxData):
                fluxArray[ii, jj] = ff
                fibrePixXArray[ii, jj] = fibrePix[ii, 0]
                fibrePixYArray[ii, jj] = fibrePix[ii, 1]

            del zoomed

        for ii, row in enumerate(fibrePix):

            flux = fluxArray[ii]
            fibrePixX = fibrePixXArray[ii]
            fibrePixY = fibrePixYArray[ii]

            fluxTable.add_row(
                (SIMBMAP['fnumdesign'][ii],
                 SIMBMAP['raoff'][ii], SIMBMAP['decoff'][ii],
                 fibrePixX, fibrePixY, flux))

        self.fluxTable = fluxTable

        del img
        return

    @staticmethod
    def zoomImage(img, factor=10):
        """Resamples an HDU by a certain factor."""

        def zoom(aa, factor):
            return np.repeat(np.repeat(aa, factor, axis=0), factor, axis=1)

        zoomed = img.copy()
        zoomed.data = zoom(zoomed.data, factor)
        zoomed.data /= (factor * factor)
        zoomed.header['CD1_1'] /= factor
        zoomed.header['CD1_2'] /= factor
        zoomed.header['CD2_1'] /= factor
        zoomed.header['CD2_2'] /= factor
        zoomed.header['CRPIX1'] *= factor
        zoomed.header['CRPIX2'] *= factor

        return zoomed

    @staticmethod
    def phot(data, centres, aperture):
        """Performs aperture photometry on the target."""

        rMax = int(aperture) + 1

        ii = np.arange(-rMax, rMax+1)
        jj = np.arange(-rMax, rMax+1)

        coords = np.dstack(np.meshgrid(ii, jj)).reshape(-1, 2)
        distances = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
        validCoords = coords[distances < aperture]

        fluxes = []
        for xx, yy in np.round(centres):
            vv = validCoords + np.array([yy, xx], np.int)
            fluxes.append(data[(vv[:, 0], vv[:, 1])].sum())

        return fluxes

    def _getFibreCoordinates(self):
        """Returns the RA,DEC coordinates for each fibre."""

        raCen = self.row['ra']
        decCen = self.row['dec']

        fibreWorld = np.zeros((len(SIMBMAP), 2), float)

        raOffDeg = SIMBMAP['raoff'] / 3600. / np.cos(decCen * np.pi / 180.)
        decOffDeg = SIMBMAP['decoff'] / 3600.

        fibreWorld[:, 0] = raCen + raOffDeg
        fibreWorld[:, 1] = decCen + decOffDeg

        return fibreWorld

    def plot(self, filename, multipage=False):
        """Plots the target, bundle and measured flux.

        This method creates a plot that includes images of the target
        (in composite colour and for each of the filters used), the position
        of the fibres in the bundle and estimations of the measured flux and
        cumulated S/N.

        The outputted plot is a grid with three columns and as many rows as
        filters have been defined, plus an extra, top row for the IRG colour
        composite image. The left panel in each row displays the target, while
        the middle panel shows the position of each fibre on top of the target.
        Fibres plotted with a white line are fibres includes in an IFU of the
        corresponding IFUDESIGNSIZE for the target. The remaining fibres up to
        a 127-fibre bundles are plotted in red. If the information is
        available, black concentric circles indicate the location of 1.5Reff
        and 2.5Reff. The right panel displays, in greyscale, the measured flux
        for each fibre. If the r-band image is included, an estimation of the
        cumulated S/N in one hour is also plotted. The right panel for the IRG
        image is different and contains text displaying the name of the target
        and, if available, the redshift, stellar mass, Reff and whether the
        target belong to the primary or secondary MaNGA samples.

        Parameters
        ----------
        filename : str
            The filename for the plot.
        multipage : False or `PdfPages` object, optional
            If False, the plot is saved to the path defined by ``filename``.
            If it's an instance of `PdfPages`, the plot is added to the
            multipage PDF file. The latter case is mostly intended to be used
            by `PlateMags.plot()` when producing multipage PDF plots with one
            target per page.

        """

        plt.cla()
        plt.clf()

        hh = self.binnedData[1].header
        ww = wcs.WCS(hh)

        # Backwards compatibility for the sake of pywcsgrid2
        ww.wcs_pix2sky = ww.wcs_pix2world
        ww.wcs_sky2pix = ww.wcs_world2pix

        scale = SCALE

        fig = plt.figure(figsize=(8.5, 11), tight_layout=True)
        grid_helper = pw2.GridHelper(wcs=ww)
        grid_helper.locator_params(nbins=3)

        grid = axes_grid.ImageGrid(
            fig, 111,
            nrows_ncols=(NFILTERS+1, 3),
            ngrids=None,
            direction='row',
            axes_pad=0.04, add_all=True,
            share_all=True, aspect=True,
            label_mode='L', cbar_mode=None,
            axes_class=(pw2.Axes, dict(grid_helper=grid_helper)))

        irgImg = plt.imread(self.preimageIRG, format='jpg')

        nRows = NFILTERS + 1
        nCols = 3

        nn = 0
        for ii in range(nRows):

            imgData = self.binnedData[3*(ii-1)+1].data

            for jj in range(nCols):

                filt = FILTERS[ii-1]

                if ii == 0:

                    if jj == 0:
                        grid[nn].imshow(irgImg, origin='lower')

                    elif jj == 1:
                        grid[nn].imshow(irgImg, origin='lower')
                        filterIdx = FILTERS.index('r')
                        self._plotBundle(grid[nn], scale,
                                         filterIdx=filterIdx)

                    else:
                        self._addCaptions(grid[nn])
                        grid[nn].axis[
                            'left', 'right',
                            'top', 'bottom'].set_visible(False)

                elif ii > 0:

                    if jj != 2:
                        grid[nn].imshow(imgData, origin='lower')

                        if jj == 0:
                            grid[nn].add_compass(loc=1)
                        elif jj == 1:
                            filterIdx = FILTERS.index(filt)
                            self._plotBundle(grid[nn], scale,
                                             filterIdx=filterIdx)

                    else:
                        grid[nn].set_axis_bgcolor('k')
                        self._plotFlux(grid[nn], scale, ii-1)
                        grid[nn].text(1.1, 0.5, r'${0}$'.format(filt),
                                      horizontalalignment='center',
                                      verticalalignment='center',
                                      transform=grid[nn].transAxes,
                                      fontsize=14)

                nn += 1

        if not multipage:
            plt.savefig(filename, format='pdf', dpi=DPI)
        else:
            multipage.savefig(dpi=DPI)

        del imgData
        del grid
        plt.close('all')

    def _plotBundle(self, ax, scale, filterIdx=0):
        """Plots the location of each fibre in the bundle on the target."""

        re = self._getRe()  # Re in arcsec

        if re is not None:

            eff15Pix = re / scale * 1.5
            eff25Pix = re / scale * 2.5

            centralFibre = self.fluxTable[
                self.fluxTable['fnumdesign'] == 1]['x', 'y'][0]
            for eff in [eff15Pix, eff25Pix]:
                ell = Ellipse(
                    (centralFibre['x'][filterIdx],
                     centralFibre['y'][filterIdx]),
                    eff * 2, eff * 2,
                    facecolor='None', edgecolor='k',
                    linewidth=0.5, linestyle='solid', zorder=10)
                ax.add_artist(ell)

        ifuSize = self.row['ifudesignsize']

        width = TARGET_SEEING / scale
        coords = self.fluxTable['x', 'y', 'fnumdesign']

        xy = [(coords[ii]['x'][filterIdx], coords[ii]['y'][filterIdx])
              for ii in range(len(coords))]
        xy = np.array(xy)
        colours = ['w' if fnum <= ifuSize else 'r'
                   for fnum in coords['fnumdesign']]

        cc = EllipseCollection(
            width, width, 0.0, units='xy', transOffset=ax.transData,
            offsets=xy, linewidth=0.7, facecolor='None', edgecolor=colours,
            zorder=20)
        ax.add_collection(cc)

    def _getIRGImage(self, refHeader):
        """Downloads and saves to disk the IRG colour image of the target."""

        if self.dataOrigin == 'NSA':
            baseURL = BASE_URL.format(self.subDir, self.IAUName)
            irgURL = baseURL + \
                '/atlases/{0}/{1}-parent-{0}-irg.jpg'.format(
                    self.pID, self.IAUName)
        elif self.dataOrigin == 'SDSS':
            run, rerun, camcol, field = self.SDSSField
            fieldDic = {'run': run, 'rerun': rerun,
                        'camcol': camcol, 'field': field}
            irgURL = SDSS_FRAME_IRG.format(**fieldDic)

        ww = wcs.WCS(refHeader)
        xx, yy = ww.wcs_world2pix(self.RA, self.Dec, 0)

        ff = urllib2.urlopen(irgURL)

        data = plt.imread(
            cStringIO.StringIO(ff.read()),
            format='jpg')[::-1, :, :]
        data = data[yy-NPIX:yy+NPIX, xx-NPIX:xx+NPIX]

        imsave(self.preimageIRG, data)

    def _plotFlux(self, ax, scale, band):
        """Displays the flux on each fibre using a greyscale."""

        width = TARGET_SEEING / scale
        bandName = FILTERS[band]

        fluxes = self.fluxTable['fluxes'][:, band]
        fluxes[fluxes < 0.] = 0.0
        fluxesNorm = fluxes / np.max(fluxes)
        fluxesStr = ['{0:.3f}'.format(flux) for flux in fluxesNorm]

        coords = self.fluxTable['x', 'y']
        xy = [(coords[ii]['x'][band], coords[ii]['y'][band])
              for ii in range(len(coords))]
        xy = np.array(xy)

        cc = EllipseCollection(
            width, width, 0.0, units='xy', transOffset=ax.transData,
            offsets=xy, linewidth=0.7, edgecolor=fluxesStr,
            facecolor=fluxesStr)
        ax.add_collection(cc)

        if bandName.lower() == 'r':
            for ii in range(len(coords)):
                xx = coords['x'][ii][band]
                yy = coords['y'][ii][band]
                sn = float(BOSS_SN(fluxes[ii]))
                ax.text(xx, yy, '{0:.1f}'.format(sn),
                        color='r', fontsize=3,
                        horizontalalignment='center',
                        verticalalignment='center')

        return

    def _getRe(self):
        """Calculates the effective radius for the target, if available."""

        # Returns the effective radius of the galaxy, if exists.

        if config['reffField'].lower() in self.row.colnames:
            re = self.row[config['reffField'].lower()]
        elif len(MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] ==
                 self.MaNGAID]) > 0:
            re = MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] ==
                              self.MaNGAID][0][config['reffField']]
        else:
            re = None

        return re

    def _addCaptions(self, ax):
        """Adds information about the target to the plot."""

        ax.text(
            0.5, 0.9, '{0}'.format(self.MaNGAID),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes, fontsize=17)

        ax.text(
            0.5, 0.035, r'Scaled fiber2flux \& S/N per hour',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes, fontsize=6)

        if 'stellar_mass' in self.row.colnames:
            sMass = self.row['stellar_mass']
        elif len(MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] == self.MaNGAID]) > 0:
            sMass = MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] ==
                                 self.MaNGAID][0]['STELLAR_MASS']
        else:
            sMass = None

        if 'z' in self.row.colnames:
            zz = self.row['z']
        elif 'redshift' in self.row.colnames:
            zz = self.row['redshift']
        elif len(MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] == self.MaNGAID]) > 0:
            zz = MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] == self.MaNGAID][0]['Z']
        else:
            zz = None

        targetSubsample = None
        # if 'MANGA_TARGET1' in self.row.colnames and \
        #         'MANGA_TARGET2' in self.row.colnames:
        #     if self.row['MANGA_TARGET1'] >= 1:
        #         targetSubsample = 'Primary'
        #     elif self.row['MANGA_TARGET2'] >= 1:
        #         targetSubsample = 'Secondary'

        re = self._getRe()

        text = ''
        if sMass is not None:
            text += r'$\log\,M_{{\rm *}}={0:.3f}$'.format(sMass) + '\n'
        if zz is not None:
            text += r'$z={0:.3f}$'.format(zz) + '\n'
        if re is not None:
            text += r'$R_{{\rm e}}={0:.2f}\,{{\rm arcsec}}$'.format(re) + \
                '\n'
        if targetSubsample is not None:
            text += '{0}'.format(targetSubsample)

        text = text.strip('\n')
        if text != '':
            ax.text(
                0.5, 0.5, text,
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
