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
from urlparse import urlparse

from Gohan.exceptions import GohanWarning, GohanError
from Gohan import log, config, readPath
from Gohan.utils import pywcsgrid2 as pw2
from Gohan.utils import getSDSSRun
from Gohan.utils.yanny import yanny, write_ndarray_to_yanny

from astropy import wcs
from astropy.io import fits
from astropy import table
from astropy.modeling import models, fitting

from scipy.ndimage.filters import gaussian_filter
from scipy.misc import imsave
from scipy import interpolate

import matplotlib
from matplotlib import pyplot as plt
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from matplotlib.patches import Ellipse
from matplotlib.collections import EllipseCollection
from matplotlib.backends.backend_pdf import PdfPages


__ALL__ = ['PlateMags']

matplotlib.rc('font', size=10)
matplotlib.rc('text', usetex=True)
plt.ioff()

FILTERS = config['plateMags']['filters']
NFILTERS = len(FILTERS)
NSAFILTERPOS = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}

try:
    MANGA_SAMPLE = table.Table.read(readPath(config['catalogues']['science']))
except:
    MANGA_SAMPLE = None
    warnings.warn('no MaNGA sample file loaded', GohanWarning)

SIMBMAP = table.Table(yanny(readPath(config['plateMags']['simbmap']),
                            np=True)['SIMBMAP'])

BOSS_SN_DATA = table.Table.read(
    readPath(config['plateMags']['BOSS_SN']),
    format='ascii.commented_header')

BOSS_SN = interpolate.interp1d(BOSS_SN_DATA['fiber2flux'][::-1],
                               BOSS_SN_DATA['sn_hour'][::-1])

NIFUS = np.sum(config['IFUs'].values())

NSA_BASE_URL = os.path.join(config['nsaImaging']['baseURL'], '{0}/')
NSA_BASE_URI = '{uri.scheme}://{uri.netloc}/'.format(
    uri=urlparse(NSA_BASE_URL))

SDSS_BASE_URL = config['sdssImaging']['baseURL']
SDSS_BASE_URI = '{uri.scheme}://{uri.netloc}/'.format(
    uri=urlparse(SDSS_BASE_URL))

# Creates authentication openers.
password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()

for imagingMode in ['nsaImaging', 'sdssImaging']:
    if config[imagingMode]['password'] != 'None':
        uri = NSA_BASE_URI if imagingMode == 'nsaImaging' else SDSS_BASE_URI
        password_mgr.add_password(realm=None, uri=uri,
                                  user=config[imagingMode]['username'],
                                  passwd=config[imagingMode]['password'])

auth_handler = urllib2.HTTPBasicAuthHandler(password_mgr)
opener = urllib2.build_opener(auth_handler)
urllib2.install_opener(opener)

SDSS_FRAME = os.path.join(
    SDSS_BASE_URL,
    'photoObj/frames/{rerun:d}/{run:d}/{camcol:d}/' +
    'frame-{filter}-{run:06d}-{camcol:d}-{field:04d}.fits.bz2')
SDSS_PSFIELD = os.path.join(
    SDSS_BASE_URL,
    'photo/redux/{rerun:d}/{run:d}/objcs/{camcol:d}/psField-{run:06d}-' +
    '{camcol:d}-{field:04d}.fit')
SDSS_PHOTOFIELD = os.path.join(
    SDSS_BASE_URL,
    'photoObj/{rerun:d}/{run:d}/photoField-{run:06d}-{camcol:d}.fits')
SDSS_FRAME_IRG = os.path.join(
    SDSS_BASE_URL,
    'photoObj/frames/{rerun:d}/{run:d}/{camcol:d}/frame-irg-{run:06d}-' +
    '{camcol:d}-{field:04d}.jpg')

IFU_MIN_SIZE = np.min(config['IFUs'].keys())
REFF = config['plateMags']['reffField'].lower()


def shortenPath(path):
    return path[:15] + '...' + path[-45:] if len(path) > 70 else path


class PlateMags(list):
    """Creates a list of PlateMagsIFU.

    This class is initialised using a PlateInput instance and returns
    a list of `PlateMagsIFU` objects, one for each IFU target in the
    PlateInput.

    Parameters
    ----------
    plateInputFile : str
        The path to the plateInput file to use.
    designid : int, optional
        The designid of the design. If None, the value is taken from the
        plateInput `designid` pair.
    kwarg : dict
        Keyword arguments to be passed to PlateMagsIFU.

    """

    def __init__(self, plateInputFile, designid=None, **kwargs):

        assert os.path.exists(plateInputFile)

        self.plateInputFile = plateInputFile

        yn = yanny(self.plateInputFile, np=True)
        self.struct1 = table.Table(
            yn[config['plateInputs']['mangaInputStructure']])
        self.designid = int(yn['designid']) if designid is None else designid

        log.info('plateInput: {0}'.format(
            os.path.basename(self.plateInputFile)))

        plateMagsIFUList = [PlateMagsIFU(row, designid=self.designid,
                                         **kwargs)
                            for row in self.struct1
                            if row['ifudesignsize'] >= IFU_MIN_SIZE]

        list.__init__(self, plateMagsIFUList)

        self._createPlateMags()

    @property
    def preImages(self):
        """Returns a list of preimage filenames for all IFUs."""
        return [ss.preImagePath for ss in self if ss.preImagePath is not None]

    def plot(self, overwrite=False, **kwargs):
        """Plots the bundle position and flux for each target.

        This method recursively calls `PlateMags.plot()` for each target
        and saves the result in a single, multipage PDF file with name
        ``filename``. If ``filename=None``, the default name is
        ``plateMags-XXXX.pdf``` with ``XXXX`` the design ID.

        Refer to `PlateMags.plot()` for details on the layout of the plot.

        Parameters
        ----------
        filename : str or None, optional
            The filename of the plot file. If None, the default is
            ``plateMags-XXXX.pdf``` with ``XXXX`` the design ID.
        useRepo : bool, optional
            If True, the plateMags plot file is create in the appropiate
            directory in $MANGACORE/plateMags. filename, if defined,
            is ignored in this case.
        repoPath : str, optional
            The path of the repository. If not defined, the mangacore value
            in configuration file is used.

        """

        plateMagsPlotPath = self._getFilePath(extension='.pdf', **kwargs)

        if os.path.exists(plateMagsPlotPath):
            if overwrite:
                os.remove(plateMagsPlotPath)
            else:
                log.info('plot {0} aready exists'.format(
                         shortenPath(plateMagsPlotPath)))
                return

        plt.cla()
        plt.clf()

        pp = PdfPages(plateMagsPlotPath)

        for ifu in self:

            if ifu.fluxTable is None:
                continue

            ifu.plot(None, multipage=pp)

        pp.close()

        log.info('plot saved as {0}'.format(shortenPath(plateMagsPlotPath)))

    def _createPlateMags(self):
        """Creates a `Table` with the plateMags columns."""

        plateMags = table.Table(
            None,
            names=['IFUDESIGN', 'MANGAID', 'FNUMDESIGN',
                   'RAOFF', 'DECOFF', 'FIBER2MAG'],
            dtype=[int, 'S20', int, 'f8', 'f8', ('f8', NFILTERS)])

        missingIFUs = 0

        for ifu in self:
            if ifu.fluxTable is None:
                missingIFUs += 1
                continue
            ifuDesign = ifu._plateInputRow['ifudesign']
            mangaID = ifu._plateInputRow['mangaid']
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

    def write(self, **kwargs):
        """Writes a plateMags file to disk.

        Parameters
        ----------
        filename : str or None, optional
            The filename of the plateMags file. If None, the default is
            ``plateMags-XXXX.par``` with ``XXXX`` the design ID.
        useRepo : bool, optional
            If True (the default), the plateMags file is create in the
            appropiate directory in $MANGACORE/plateMags. filename, if defined,
            is ignored in this case.
        repoPath : str, optional
            The path of the repository. If not defined, the mangacore value
            in configuration file is used.

        """

        if len(self.plateMags) == 0:
            raise ValueError('plateMags has not values.')

        plateMagsPath = self._getFilePath(extension='.par', **kwargs)
        if os.path.exists(plateMagsPath):
            os.remove(plateMagsPath)

        write_ndarray_to_yanny(plateMagsPath, self.plateMags.as_array(),
                               structname='PLATEMAGS')

        if self._missingIFUs > 0:
            warnings.warn('no data for {0} IFUs. '.format(
                self._missingIFUs) + 'The remaining data has been saved.',
                GohanWarning)

        log.info('plateMags saved as {0}'.format(shortenPath(plateMagsPath)))

        return plateMagsPath

    def _getFilePath(self, filename=None, useRepo=True,
                     repoPath=config['mangacore'], extension='.par', **kwargs):
        """Returns the path to use by PlateMags.write and PlateMags.plot."""

        if useRepo:
            plateMagsPath = os.path.join(
                readPath(repoPath),
                'platedesign/platemags/'
                'D00{0:s}XX/plateMags-{1:04d}{2}').format(
                    str(self.designid)[0:2], self.designid, extension)
        else:
            if filename is not None:
                plateMagsPath = filename
            else:
                plateMagsPath = 'plateMags-{0:04d}{1}'.format(self.designid,
                                                              extension)
        plateMagsPath = os.path.realpath(plateMagsPath)

        if not os.path.exists(os.path.dirname(plateMagsPath)):
            os.makedirs(os.path.dirname(plateMagsPath))

        return plateMagsPath


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

    The filters used are defined in `Gohan.defaults['plateMags']['filters']`
    and must always be a list containing a subset of [u, g, r, i, z].

    The data is then rebinned to match the spatial resolution of one
    fibre (defined in  `Gohan.defaults['plateMags']['targetSeeing']`)
    and fluxes for each filter are measured at the position of each fibre.
    Note that in any case a 127-fibre bundle is assumed,
    regardless of the actual IFUDESGNSIZE.

    The class is not intended to be used directly but to be called
    from `PlateMags`.

    Parameters
    ----------
    struct1Row : `astropy.table.Row`
        A row of a `Table` containing the plateInput information for
        a single target.
    useRepo : bool, optional
        If True, the default, preimaging is created at the appropiate location
        in the repositiory defined by `repoPath`.
    repoPath : bool, optional
        The path of the repository. Deafaults to the value of
        `Gohan.default['preimaging']`.
    designid : int, optional
        The designid of the design this IFU belongs to. If useRepo is True,
        this field is mandatory.
    keepTemporary : bool, optional
        If False (the default), intermediate step imaging is deleted
        when not needed anymore.
    overwrite : bool, optional
        If True, imaging data is redownloaded and rewritten even if imaging
        with matching filenames is found. Otherwise, the already existing data
        is read.

    """

    def __init__(self, struct1Row, useRepo=True, designid=None,
                 repoPath=config['preimaging'], **kwargs):

        self._plateInputRow = struct1Row
        colnames = [cc.lower() for cc in self._plateInputRow.colnames]

        self.RA = self._plateInputRow[colnames.index('ra')]
        self.Dec = self._plateInputRow[colnames.index('dec')]
        self.mangaid = self._plateInputRow[colnames.index('mangaid')].strip()

        log.debug('MANGAID: {0}'.format(self.mangaid))

        if useRepo:
            if designid is None:
                raise GohanError('useRepo=True but designid=None.')
            preImageDir = os.path.join(
                readPath(repoPath),
                'D00{0:s}XX/{1:04d}').format(str(designid)[0:2], designid)
        else:
            preImageDir = './'

        preimageData = self._getPreimageData(preImageDir=preImageDir,
                                             **kwargs)

        if preimageData is not None:
            self.binnedData = self.binData(preimageData)
            preimageData.close()  # Not needed anymore
            self.fluxTable = self._calculateFluxes(self.binnedData)
        else:
            self.fluxTable = None
            log.important(
                'no available data for {0}'.format(self.mangaid.strip()))

    def removeFile(self, file):
        if file is not None and os.path.exists(file):
            os.remove(file)

    def _getPreimageData(self, **kwargs):
        """Creates the preimage file and downloads IRG imaging."""

        self.preImagePath = None
        self.preImageIRG = None
        self.source = None

        preImageDir = kwargs.get('preImageDir', './')
        overwrite = kwargs.get('overwrite', False)
        keepTemporary = kwargs.get('keepTemporary', False)

        _preImagePath = os.path.join(
            preImageDir, 'preimage-{0}.fits.gz'.format(self.mangaid))
        _preImageIRG = os.path.join(
            preImageDir, 'preimage-{0}_irg.jpg'.format(self.mangaid))

        if os.path.exists(_preImagePath) and not overwrite:

            self.preImagePath = _preImagePath
            preimageData = fits.open(self.preImagePath)
            log.debug('using previously saved preimage data.')

            if 'DSOURCE' in preimageData[1].header:
                self.source = preimageData[1].header['DSOURCE']
            if os.path.exists(_preImageIRG):
                self.preImageIRG = _preImageIRG
            else:
                raise GohanError('IRG image not downloaded. This is probably '
                                 'due to half-downloaded data. Remove the '
                                 'preimaging directory for this design and '
                                 'try again.')

            return preimageData

        else:

            output = self._getFromNSA(**kwargs)
            source = 'NSA'
            if output is False:
                output = self._getFromSDSS(**kwargs)
                source = 'SDSS'
                if output is False:
                    return None

        data = fits.open(output)
        preimageData = self._cropHDUs(data, source=source)

        if not keepTemporary:
            self.removeFile(output)

        if not os.path.exists(preImageDir):
            os.makedirs(preImageDir)

        preimageData.writeto(_preImagePath)

        self.preImagePath = _preImagePath
        self.preImageIRG = _preImageIRG
        self.source = source

        # It has to be r because the astrometry of the IRG image
        # in SDSS matches the one in the r frame. In NSA all frames
        # are aligned
        rFilt = FILTERS.index('r')
        rIdx = 3*rFilt+1
        self._getIRGImage(_preImageIRG, data[rIdx].header)

        data.close()
        return preimageData

    def _getNSAParams(self):
        """Determines if enough data is present to download NSA imaging."""

        data = self._plateInputRow
        cols = [col.lower() for col in data.colnames]
        neededCols = ['iauname', 'subdir', 'pid']

        returnValues = []
        for neededCol in neededCols:

            if neededCol not in cols:
                return None

            value = data[cols.index(neededCol)]
            if value in [-999, '-999', '-999.', 'NULL']:
                return None

            if neededCol == 'pid':
                value = int(value)

            returnValues.append(value)

        return returnValues

    def _getFromNSA(self, **kwargs):
        """Main routine to download NSA data for the target."""

        log.debug('... Trying to get NSA data.')

        preImageDir = kwargs.get('preImageDir', './')
        overwrite = kwargs.get('overwrite', False)

        nsaParams = self._getNSAParams()

        if nsaParams is not None:

            self.IAUName = nsaParams[0]
            self.subDir = nsaParams[1]
            self.pID = nsaParams[2]

        else:

            log.debug('...... {0} not found in NSA.'.format(
                      self.mangaid.strip()))
            return False

        imageName = os.path.join(preImageDir,
                                 '{0}_NSA.fits'.format(self.IAUName))

        if not os.path.exists(imageName) or overwrite:

            output = self.getNSAImage(
                imageName, self.IAUName, self.subDir, self.pID)

            if output is False:
                warnings.warn(
                    '{0} has no NSA data.'.format(self.mangaid),
                    GohanWarning)
                return False

        else:

            log.debug('...... {0} already exists.'.format(
                      shortenPath(imageName)))

        return imageName

    def getNSAImage(self, imageName, IAUName, subDir, pID, **kwargs):
        """Downloads NSA imaging and creates the parent FITS file."""

        baseURL = NSA_BASE_URL.format(subDir)

        parentURL = baseURL + \
            '/atlases/{0}/{1}-parent-{0}.fits.gz'.format(pID, IAUName)
        varURL = baseURL + \
            '/atlases/{0}/{1}-ivar-{0}.fits.gz'.format(pID, IAUName)
        bpsfURL = [baseURL + '/{0}-{1}-bpsf.fits.gz'.format(IAUName, band)
                   for band in FILTERS]

        hdus = []
        for url in [parentURL] + [varURL] + bpsfURL:
            hdus.append(self._getNSAHDU(url, **kwargs))
            if None in hdus:
                return False

        if not os.path.exists(os.path.dirname(imageName)):
            os.makedirs(os.path.dirname(imageName))

        fullHDU = self.concatenateNSAHDUs(hdus)

        fullHDU.writeto(imageName, clobber=True)
        fullHDU.close()

        for ii in hdus:
            ii.close()
        del hdus

        return True

    def _getNSAHDU(self, url, **kwargs):
        """Gets a NSA file."""

        preImageDir = kwargs.get('preImageDir', './')
        keepTemporary = kwargs.get('keepTemporary', False)

        ff = tempfile.NamedTemporaryFile(
            dir=preImageDir, suffix='.fits.gz', delete=False)

        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError:
            log.debug('...... URL not found.')
            self.removeFile(ff.name)
            return None

        ff.write(response.read())
        ff.close()

        hduTmp = fits.open(ff.name)

        if not keepTemporary:
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

    def _cropHDUs(self, hdu, source='UNKNWN'):
        """Trims and HDU to the size of an IFU."""

        log.debug('... Creating preimage.')

        newHDU = self.copyHDUList(hdu)

        ra = self.RA
        dec = self.Dec

        for ii in range(NFILTERS):

            imgIdx = ii*3 + 1
            varIdx = ii*3 + 2
            psfIdx = ii*3 + 3

            ww = wcs.WCS(newHDU[imgIdx].header)
            xx, yy = ww.wcs_world2pix(ra, dec, 0)

            nPix = config['plateMags']['nPix']
            xmin = int(xx) - nPix
            # xmax = int(xx) + nPix
            ymin = int(yy) - nPix
            # ymax = int(yy) + nPix

            for nn in [imgIdx, varIdx, psfIdx]:
                newHDU[nn].header.set(
                    'DSOURCE', source, 'The source of the data')

            for nn in [imgIdx, varIdx]:
                # Creates a frame for the image for the case its shape < 2*nPix
                data = newHDU[nn].data.copy()
                tmpData = np.zeros((data.shape[0]+nPix*2,
                                    data.shape[1]+nPix*2), np.float)
                tmpData[nPix:tmpData.shape[0]-nPix,
                        nPix:tmpData.shape[1]-nPix] = data

                tmpData = tmpData[yy:yy+2*nPix, xx:xx+2*nPix]

                newHDU[nn].data = np.array(tmpData)

                newHDU[nn].header['CRPIX1'] -= xmin
                newHDU[nn].header['CRPIX2'] -= ymin

        return newHDU

    def _getFromSDSS(self, **kwargs):
        """Determines the SDSS field for the target and gets the data."""

        log.debug('... Trying to get SDSS data.')

        preImageDir = kwargs.get('preImageDir', './')
        overwrite = kwargs.get('overwrite', False)

        self.SDSSField = getSDSSRun(self.RA, self.Dec)

        if self.SDSSField is not None:
            run, rerun, camcol, field = self.SDSSField
        else:
            warnings.warn('galaxy not in SDSS footprint.', GohanWarning)
            return False

        imageName = os.path.join(preImageDir,
                                 'Field_{0}_{1}_{2}_SDSS.fits'.format(
                                     run, camcol, field))

        if os.path.exists(imageName) and not overwrite:
            log.debug('...... {0} already exists.'.format(
                      shortenPath(imageName)))

        else:

            try:
                self._getSDSSImage(imageName, run, rerun, camcol, field,
                                   **kwargs)
            except:
                raise GohanError('...... images could not be downloaded.')

        return imageName

    def _getSDSSImage(self, imageName, run, rerun, camcol, field, **kwargs):
        """Download the SDSS images.

        SDSS frames are downloaded for each filter, as well as psField files
        used for PSF reconstruction. The final image is a multiextension FITS
        file with the same format as the preimage.

        """

        preImageDir = kwargs.get('preImageDir', './')
        keepTemporary = kwargs.get('keepTemporary', False)

        fieldDic = {'run': run, 'rerun': rerun,
                    'camcol': camcol, 'field': field}

        frameFilenames = []
        for filt in FILTERS:
            fieldDic.update({'filter': filt})
            frameFilename = self._getTmpFilename(dir=preImageDir)
            urllib.urlretrieve(SDSS_FRAME.format(**fieldDic), frameFilename)
            frameFilenames.append(frameFilename)

        psFieldFilename = self._getTmpFilename(dir=preImageDir, suffix='.fit')
        urllib.urlretrieve(SDSS_PSFIELD.format(**fieldDic), psFieldFilename)

        photoFieldFilename = self._getTmpFilename(dir=preImageDir,
                                                  suffix='.fits')
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

        if not os.path.exists(os.path.dirname(imageName)):
            os.makedirs(os.path.dirname(imageName))

        image.writeto(imageName, clobber=True)

        # Removes temporary files.
        if not keepTemporary:
            for img in frameFilenames + [photoFieldFilename, psFieldFilename]:
                self.removeFile(img)

        image.close()

        return True

    def _getTmpFilename(self, dir='./', suffix='.fits.bz2'):
        """Creates a temporary file and returns its filename."""

        tmpFile = tempfile.NamedTemporaryFile(
            dir=dir, suffix='.fits.bz2', delete=True)
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
        """Bins the data to targetSeeing using a Gaussian kernel."""

        log.debug('... Rebinning images.')

        rebin = self.copyHDUList(data)

        for ii in range(NFILTERS):

            img = data[3*ii+1]
            psf = data[3*ii+3]

            if self.source == 'NSA':
                scale = config['nsaImaging']['scale']
            elif self.source == 'SDSS':
                scale = config['sdssImaging']['scale']

            seeingPix = self.getSeeing(psf)
            seeingArcSec = scale * seeingPix
            targetSeeing = config['plateMags']['targetSeeing']

            if seeingArcSec < targetSeeing:

                log.debug(
                    '...... Band ' +
                    '{0} from seeing {1:.3f} arcsec to {2} arcsec.'.format(
                        FILTERS[ii], seeingArcSec, targetSeeing))

                targetSigma = targetSeeing / (2*np.sqrt(2*np.log(2))) / scale
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
        if self.source == 'NSA':
            scale = config['nsaImaging']['scale']
        elif self.source == 'SDSS':
            scale = config['sdssImaging']['scale']
        scaleZoomed = scale / factor

        # Radius of one fibre in pixels in the resampled image.
        radius = config['plateMags']['targetSeeing'] / 2. / scaleZoomed

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

        del img
        return fluxTable

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

        raCen = self._plateInputRow['ra']
        decCen = self._plateInputRow['dec']

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

        if self.source == 'NSA':
            scale = config['nsaImaging']['scale']
        elif self.source == 'SDSS':
            scale = config['sdssImaging']['scale']

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

        irgImg = plt.imread(self.preImageIRG, format='jpg')

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
            plt.savefig(filename, format='pdf', dpi=config['plateMags']['DPI'])
        else:
            multipage.savefig(dpi=config['plateMags']['DPI'])

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

        ifuSize = self._plateInputRow['ifudesignsize']

        width = config['plateMags']['targetSeeing'] / scale
        coords = self.fluxTable['x', 'y', 'fnumdesign']

        xy = [(coords[ii]['x'][filterIdx], coords[ii]['y'][filterIdx])
              for ii in range(len(coords))]
        xy = np.array(xy)
        colours = ['w' if fnum <= ifuSize else 'r'
                   for fnum in coords['fnumdesign']]

        cc = EllipseCollection(
            width, width, 0.0, units='x', offsets=xy,
            transOffset=ax.transData,
            linewidth=0.7, facecolor='None', edgecolor=colours,
            zorder=20)
        ax.add_collection(cc)

    def _getIRGImage(self, imageName, refHeader):
        """Downloads and saves to disk the IRG colour image of the target."""

        if self.source == 'NSA':
            baseURL = NSA_BASE_URL.format(self.subDir, self.IAUName)
            irgURL = baseURL + \
                '/atlases/{0}/{1}-parent-{0}-irg.jpg'.format(
                    self.pID, self.IAUName)
        elif self.source == 'SDSS':
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

        nPix = config['plateMags']['nPix']

        # Creates a frame for the image for the case its shape < 2*nPix
        tmpData = np.zeros((data.shape[0]+nPix*2,
                            data.shape[1]+nPix*2, 3), np.float)

        tmpData[nPix:tmpData.shape[0]-nPix,
                nPix:tmpData.shape[1]-nPix, :] = data

        tmpData = tmpData[yy:yy+2*nPix, xx:xx+2*nPix, :]

        imsave(imageName, tmpData)

    def _plotFlux(self, ax, scale, band):
        """Displays the flux on each fibre using a greyscale."""

        width = config['plateMags']['targetSeeing'] / scale
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
            width, width, 0.0, units='x', transOffset=ax.transData,
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

        if REFF in self._plateInputRow.colnames:
            return self._plateInputRow[REFF]

        if MANGA_SAMPLE is None:
            return None

        sampleRow = MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] == self.mangaid]

        if len(sampleRow) == 0:
            return None
        if REFF not in sampleRow.colnames:
            return None

        return sampleRow[REFF]

    def _addCaptions(self, ax):
        """Adds information about the target to the plot."""

        if MANGA_SAMPLE is None:
            sampleRow = None
        else:
            sampleRow = MANGA_SAMPLE[MANGA_SAMPLE['MANGAID'] == self.mangaid]
            if len(sampleRow) == 0:
                sampleRow = None

        ax.text(
            0.5, 0.9, '{0}'.format(self.mangaid),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes, fontsize=17)

        ax.text(
            0.5, 0.035, r'Scaled fiber2flux \& S/N per hour',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes, fontsize=6)

        if 'stellar_mass' in self._plateInputRow.colnames:
            sMass = self._plateInputRow['stellar_mass']
        elif sampleRow is not None and 'STELLAR_MASS' in sampleRow.colnames:
            sMass = sampleRow['STELLAR_MASS']
        else:
            sMass = None

        if 'z' in self._plateInputRow.colnames:
            zz = self._plateInputRow['z']
        elif 'redshift' in self._plateInputRow.colnames:
            zz = self._plateInputRow['redshift']
        elif sampleRow is not None and 'Z' in sampleRow.colnames:
            zz = sampleRow['Z']
        else:
            zz = None

        re = self._getRe()

        text = ''
        if sMass is not None:
            text += r'$\log\,M_{{\rm *}}={0:.3f}$'.format(sMass) + '\n'
        if zz is not None:
            text += r'$z={0:.3f}$'.format(zz) + '\n'
        if re is not None:
            text += r'$R_{{\rm e}}={0:.2f}\,{{\rm arcsec}}$'.format(re) + \
                '\n'

        text += r'ifudesign={0:d}'.format(self._plateInputRow['ifudesign']) + \
            '\n'

        text = text.strip('\n')
        if text != '':
            ax.text(
                0.5, 0.5, text,
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
