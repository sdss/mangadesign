#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Oct 11, 2017
# @Filename: ic342.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import pathlib
import re
import warnings

import astropy.io.fits
import astropy.table
import astropy.visualization
import astropy.wcs

import photutils

import click
import numpy as np
import PIL

from Gohan import log, config
from Gohan import utils
from Gohan.utils import yanny

import Gohan.utils.image


class PS1Imaging(object):
    """Creates PS1 preimaging.

    Parameters:
        coords (tuple):
            A tuple with the *(RA, Dec)* coordinates of the centre of the
            preimage.
        ps1_data (str or `~pathlib.Path`):
            The path to a directory containing PS1 imaging. The directory
            must contain, for each skycell, the stack, weight stack, mask, and
            psf images.
        width (int):
            The width of the preimage, in pixels. PS1 pixels are 0.25 arcsec.
        verbose (bool):
            If *True*, outputs to the log.

    """

    bands = 'griz'

    def __init__(self, coords, ps1_data, width=150, verbose=False):

        self.ps1_path = pathlib.Path(ps1_data)
        assert self.ps1_path.exists()

        self.coords = coords
        self.verbose = verbose

        self._width = width

        self.preimage = self.create_preimage()

        # Creates the IRG image. Uses Q=100, stretch=0.5 to highlight galaxy features
        self.irg = PIL.Image.fromarray(
            astropy.visualization.make_lupton_rgb(self.preimage['i img'].data,
                                                  self.preimage['r img'].data,
                                                  self.preimage['g img'].data,
                                                  Q=10, stretch=0.5))

    def create_preimage(self):
        """Creates a |HDUList| with the preimage information."""

        ra, dec = self.coords

        preimaging_footprint = self.get_footprint(ra, dec, self._width)
        stack_image_base = self.get_stack_image(preimaging_footprint)

        if self.verbose:
            log.debug('using stack image {}.'.format(stack_image_base))

        preimage = astropy.io.fits.HDUList([astropy.io.fits.PrimaryHDU()])

        for band in self.bands:

            stack_image = self.ps1_path / (stack_image_base + f'.{band}.unconv.fits')
            mask_image = self.ps1_path / (stack_image_base + f'.{band}.unconv.mask.fits')
            wt_image = self.ps1_path / (stack_image_base + f'.{band}.unconv.wt.fits')
            # psf_file = self.ps1_path / (stack_image_base + f'.{band}.target.psf')

            assert stack_image.exists() and wt_image.exists()  # and psf_file.exists()

            stack_cutout = self.get_cutout(stack_image, (ra, dec), self._width)
            mask_cutout = self.get_cutout(mask_image, (ra, dec), self._width)
            wt_cutout = self.get_cutout(wt_image, (ra, dec), self._width)

            stack_preimage = self.hdu_to_nanomaggies(stack_cutout)
            stack_preimage.header['FILTNAME'] = band
            stack_preimage.header['FILTER'] = self.bands.index(band) + 1
            stack_preimage.header['EXTNAME'] = f'{band} img'
            stack_preimage.header['DSOURCE'] = 'PS1'

            preimage.append(stack_preimage)

            ivar_preimage = self.hdu_to_nanomaggies(wt_cutout)

            # Converts variance to ivar
            ivar_preimage.data = 1. / ivar_preimage.data

            # Makes ivar zero for masked values
            for ii in [1, 2, 3, 4, 5, 6, 8, 13]:
                ivar_preimage.data[mask_cutout.data == 1 << ii] = 0.

            ivar_preimage.header['FILTNAME'] = band
            ivar_preimage.header['FILTER'] = self.bands.index(band) + 1
            ivar_preimage.header['EXTNAME'] = f'{band} ivar'
            ivar_preimage.header['DSOURCE'] = 'PS1'

            preimage.append(ivar_preimage)

            # We assume the PSF is 1 arcsec.
            gaussian_2D_kernel = Gohan.utils.image.gaussian_kernel_from_fwhm(1., pixel_scale=0.25)

            psf_preimage = astropy.io.fits.ImageHDU(data=gaussian_2D_kernel.array)
            psf_preimage.header['EXTNAME'] = f'{band} psf'
            psf_preimage.header['DSOURCE'] = 'PS1'
            preimage.append(psf_preimage)

        return preimage

    def get_footprint(self, ra, dec, width, pixel_scale=0.25):
        """Returns the preimaging footprint around some coordinates."""

        width_arcsec = width * pixel_scale / 3600.
        width_ra = width_arcsec / np.cos(np.deg2rad(dec))

        return np.array([[ra - width_ra / 2., dec - width_arcsec / 2.],
                         [ra + width_ra / 2., dec - width_arcsec / 2.],
                         [ra + width_ra / 2., dec + width_arcsec / 2.],
                         [ra - width_ra / 2., dec + width_arcsec / 2.]])

    def get_stack_image(self, preimaging_footprint):
        """Checks PS1 imaging to find an image that contains the footprint."""

        stacks = self.ps1_path.glob('rings.v3.skycell.*.stk.g.unconv.fits')

        for stack in stacks:

            with warnings.catch_warnings():

                warnings.simplefilter('ignore', astropy.wcs.FITSFixedWarning)

                stack_hdu = astropy.io.fits.open(stack)
                stack_wcs = astropy.wcs.WCS(stack_hdu[1].header)

                # Pixel positions in the stack image of the preimaging footprint
                stack_image_pix = stack_wcs.wcs_world2pix(preimaging_footprint, 0)

                if (np.any(stack_image_pix[:, 0] < 0) or
                        np.any(stack_image_pix[:, 0] > stack_hdu[1].data.shape[1]) or
                        np.any(stack_image_pix[:, 1] < 0) or
                        np.any(stack_image_pix[:, 1] > stack_hdu[1].data.shape[0])):
                    continue

                stack_prefix = re.match('.+(rings.v3.skycell.[0-9]+.[0-9]+.stk).*',
                                        str(stack)).groups(0)[0]

                return stack_prefix

        ra, dec = self.coords

        raise ValueError('failed to find stack image for coords ({:.4f}, {:.4f}).'.format(ra, dec))

    def get_cutout(self, path, coords, width):
        """Returns a cropped HDU."""

        hdulist = astropy.io.fits.open(path)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', astropy.wcs.FITSFixedWarning)
            wcs = astropy.wcs.WCS(hdulist[1].header)

        center_pix = wcs.wcs_world2pix([coords], 0)[0]

        xmin = int(center_pix[0] - width / 2.)
        ymin = int(center_pix[1] - width / 2.)

        new_hdu, new_wcs = Gohan.utils.image.crop_hdu(hdulist[1], xmin, xmin + width,
                                                      ymin, ymin + width)

        assert new_wcs is not None

        new_hdu = Gohan.utils.image.replace_wcs(new_hdu, new_wcs)

        # PS1 images use compressed HDUs. Let's use normal ones.
        image_hdu = astropy.io.fits.ImageHDU(data=new_hdu.data, header=new_hdu.header)

        return image_hdu

    def hdu_to_nanomaggies(self, hdu):
        """Converts the data in the hdu to nanomaggies."""

        aa = 2.5 / np.log(10)
        xx = hdu.data / aa
        flux = hdu.header['BOFFSET'] + hdu.header['BSOFTEN'] * 2 * np.sinh(xx)

        # Converts to nanomaggies
        flux /= (10 * hdu.header['EXPTIME'])

        for key in ['BOFFSET', 'BSOFTEN', 'BSCALE', 'BZERO', 'BLANK']:
            if key in hdu.header:
                del hdu.header[key]

        hdu.data = flux
        hdu.header['FLUXUNIT'] = 'nanomaggies'

        return hdu

    def write(self, path, overwrite=True):
        """Write the preimage |HDUList| to disk."""

        path = pathlib.Path(path)

        if path.exists():
            path.unlink()

        self.preimage.writeto(str(path))

    def get_platemags(self, target_fwhm=2.5):
        """Returns an `~astropy.table.Table` of plateMags for this target.

        Rebins the data to a ``target_fwhm`` by measuring the PSF for each
        band and convolving the image with a Gaussian kernel. If
        ``target_fwhm=False``, not rebinning is applied.

        Note that the plateMags returned does not contain the ``MANGAID`` or
        ``IFUDESIGN`` columns. Those must be added manually.

        """

        platemags = astropy.table.Table(
            None, names=['FNUMDESIGN', 'RAOFF', 'DECOFF', 'FIBER2MAG'],
            dtype=[int, 'f8', 'f8', ('f8', len(self.bands))])

        raoff, decoff = utils.get_fibre_offsets().T
        flux = np.zeros((len(raoff), len(self.bands)), dtype=np.float)

        for band in self.bands:

            wcs = astropy.wcs.WCS(self.preimage[f'{band} img'].header)
            data = self.preimage[f'{band} img'].data

            # Zeros out pixels with zero ivar
            data[self.preimage[f'{band} ivar'].data == 0] = 0.

            fibre_coords = utils.get_fibre_offsets(self.coords)
            fibre_xy = wcs.wcs_world2pix(fibre_coords, 0)

            psf = self.preimage[f'{band} psf'].data

            pixel_scale = np.abs(wcs.wcs.cdelt[0]) * 3600

            if target_fwhm:

                gauss = Gohan.utils.image.fit_gaussian(psf)
                sigma = np.sqrt(gauss.x_stddev * gauss.y_stddev)
                fwhm = Gohan.utils.image.sigma_to_fwhm(sigma)

                if self.verbose:
                    log.debug(f'{band}-band: measured FWHM '
                              '{:.2f} arcsec.'.format(fwhm * pixel_scale))

                if fwhm < target_fwhm:

                    target_sigma = Gohan.utils.image.fwhm_to_sigma(target_fwhm) / pixel_scale
                    kernel_sigma = np.sqrt(target_sigma**2 - sigma**2)
                    data = Gohan.utils.image.gaussian_filter(kernel_sigma, data)

            radius = config['plateMags']['targetSeeing'] / 2. / pixel_scale

            circular_apertures = photutils.CircularAperture(fibre_xy, r=radius)
            phot_table = photutils.aperture_photometry(data, circular_apertures)

            ap_flux = phot_table['aperture_sum']

            # Sets <= zero values to nan. Uses a mask because there may already be nans.
            ap_flux_mask = ~np.isnan(ap_flux)
            ap_flux_mask[ap_flux_mask] &= ap_flux[ap_flux_mask] <= 0
            ap_flux[ap_flux_mask] = np.nan

            flux[:, self.bands.index(band)] = ap_flux * 1e-9  # In maggies

        for ii in range(len(raoff)):
            psfmag = -2.5 * np.log10(flux[ii])
            platemags.add_row((ii + 1, np.round(raoff[ii], 4), np.round(decoff[ii], 4), psfmag))

        return platemags


def get_plate_coordinates(plateid):
    """Returns a dictionary {mangaid: (ra, dec, ifudesign)}."""

    plate_holes = utils.getPlateHolesSortedData(plateid)[0]

    targets = plate_holes[(plate_holes['targettype'] == 'science') &
                          (plate_holes['mangaid'] != '')]

    assert len(targets) == 17, 'invalid number of targets'

    return_dict = {}
    for target in targets:
        return_dict[target['mangaid']] = (target['target_ra'], target['target_dec'],
                                          target['ifudesign'])

    return return_dict


@click.group()
def ic342(verbose=False):
    """Handles IC342 plates."""

    pass


@ic342.command()
@click.argument('platerun', type=click.STRING)
@click.argument('ps1_data', type=click.Path(exists=True))
@click.option('-m', '--platemags', is_flag=True, help='Generates plateMags.')
@click.pass_obj
def preimaging(obj, platerun, ps1_data, platemags=True):
    """Generates preimaging for IC342.

    A PLATERUN and the PS1_DATA path to the PS1 imaging directory must
    be provided as arguments.

    """

    verbose = obj['verbose']

    plates = utils.getFromPlatePlans(platerun, column='plateid')

    assert len(plates) > 0, 'no plates found for platerun {!r}'.format(platerun)

    for plate in plates:

        designid = utils.getDesignID(plate)

        if verbose:
            log.important('creating PS1 preimaging for plate {}.'.format(plate))

        targets = get_plate_coordinates(plate)

        platemags_list = []

        preimdir = pathlib.Path(f'./{designid:d}')
        preimdir.mkdir(exist_ok=True)

        for mangaid in targets:

            if verbose:
                log.info('running mangaid={!r}.'.format(mangaid))

            coords = targets[mangaid][0:2]
            ifudesign = targets[mangaid][2]

            ps1_preimage = PS1Imaging(coords, ps1_data, verbose=verbose)
            ps1_preimage_path = preimdir / f'preimage-{mangaid}.fits.gz'
            ps1_preimage.write(str(ps1_preimage_path))

            if verbose:
                log.info(f'saved preimage as {ps1_preimage_path!r}')

            # Saves the IRG image
            ps1_preimage.irg.save(ps1_preimage_path.with_suffix('._irg.jpg'))

            platemags = ps1_preimage.get_platemags()
            platemags.add_columns(
                [astropy.table.Column([ifudesign] * len(platemags), name='IFUDESIGN'),
                 astropy.table.Column([mangaid] * len(platemags), name='MANGAID')], [0, 0])

            platemags_list.append(platemags)

        platemags = astropy.table.vstack(platemags_list)

        platemags_path = pathlib.Path(f'plateMags-{designid}.par')
        if platemags_path.exists():
            platemags_path.unlink()

        yanny.write_ndarray_to_yanny(str(platemags_path), platemags.as_array(),
                                     structname='PLATEMAGS',
                                     hdr={'preimaging_version ': 'v2'})

        if verbose:
            log.info(f'saved plateMags as {platemags_path!s}')
