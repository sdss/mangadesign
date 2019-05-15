#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2017-10-012
# @Filename: general.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
#
# @Last modified by: José Sánchez-Gallego (gallegoj@uw.edu)
# @Last modified time: 2018-08-17 14:38:20


from __future__ import absolute_import, division, print_function

import itertools
import re
import warnings

import astropy.convolution
import astropy.modeling
import astropy.modeling.models
import astropy.wcs
import numpy


__all__ = ['crop_hdu', 'replace_wcs', 'fwhm_to_sigma', 'sigma_to_fwhm',
           'gaussian_kernel_from_fwhm', 'gaussian_filter', 'fit_gaussian',
           'CCD', 'SyntheticImage']


def crop_hdu(hdu, xmin, xmax, ymin, ymax, return_wcs=True, ignore_warnings=True):
    """Crops a HDU.

    Trims the input HDU and, optionally, calculates the WCS of the resulting
    image. It returns a new `~astropy.io.fits.ImageHDU` object with the cropped
    image. If ``return_wcs=True`` (the default), also returns the WCS
    definition for the cropped image.

    Parameters:
        hdu (~astropy.io.fits.ImageHDU):
            The `~astropy.io.fits.ImageHDU` that will be cropped.
        xmin,xmax,ymin,ymax (int):
            The section of ``hdu`` to crop.
        return_wcs (bool):
            If *True*, and the input HDU contains WCS information, will return
            a `~astropy.wcs.WCS` object with the WCS definition for the cropped
            image. Returns `None` if the input HDU does not contain WCS
            information.
        ignore_warnings (bool):
            If *True*, warnings raised during the creation of the WCS object
            will be silenced.

    """

    new_hdu = hdu.copy()
    hdu_shape = new_hdu.data.shape

    assert xmin > 0 and ymin > 0 and xmax < hdu_shape[1] and ymax < hdu_shape[0], \
        'invalid crop region.'

    data = new_hdu.data.copy()
    data = data[ymin:ymax, xmin:xmax]

    new_hdu.data = numpy.array(data)

    if return_wcs is False:
        return new_hdu

    with warnings.catch_warnings():

        if ignore_warnings:
            warnings.simplefilter('ignore', astropy.wcs.FITSFixedWarning)

        # Checks whether this image has WCS information
        wcs_list = astropy.wcs.find_all_wcs(new_hdu.header)
        if len(wcs_list) == 0:
            return new_hdu, None
        else:
            new_wcs = wcs_list[0].deepcopy()

    new_wcs.wcs.crpix[0] -= xmin
    new_wcs.wcs.crpix[1] -= ymin

    return new_hdu, new_wcs


def replace_wcs(hdu, wcs):
    """Replaces WCS information in the header.

    Removes the current WCS information in the input
    `~astropy.io.fits.ImageHDU` and replace it with new one.

    Parameters
    ----------
    hdu : `~astropy.io.fits.ImageHDU`:
        The `~astropy.io.fits.ImageHDU` whose WCS definition we will
        replace.
    wcs : `~astropy.wcs.WCS`:
        A `~astropy.wcs.WCS` object containing the new WCS definition.

    """

    # Checks for old WCS keys in the form PC001002
    pc_old_pattern = re.compile('PC0*[0-9]{1}0*[0-9]{1}')
    header_keys = hdu.header.keys()
    pc_old_in_header = filter(pc_old_pattern.match, header_keys)

    wcs_keys = wcs.to_header().keys()

    for key in itertools.chain(wcs_keys, pc_old_in_header):
        if key in hdu.header:
            del hdu.header[key]

    # Adds the new WCS header to the hdu
    hdu.header.extend(wcs.to_header().cards)

    return hdu


def fwhm_to_sigma(fwhm):
    """Returns the sigma for a FWHM."""

    return fwhm / 2 / numpy.sqrt(2 * numpy.log(2))


def sigma_to_fwhm(sigma):
    """Returns the FWHM for a sigma."""

    return sigma * 2 * numpy.sqrt(2 * numpy.log(2))


def gaussian_kernel_from_fwhm(fwhm, pixel_scale=1, **kwargs):
    """Returns a Gaussian kernel for a FWHM.

    Parameters
    ----------
    fwhm : `float`
        The FWHM (seeing) of the Gaussian kernel, in arcsec.
    pixel_scale : `float`
        The pixels scale, in arcsec.
    kwargs : `dict`
        Other parameters to be passed to
        `~astropy.convolution.Gaussian2DKernel`.

    Returns
    -------
    kernel : `~astropy.convolution.Gaussian2DKernel`
        An astropy `~astropy.convolution.Gaussian2DKernel` kernel for the
        input FHWM.

    """

    stddev = fwhm_to_sigma(fwhm) / pixel_scale

    return astropy.convolution.Gaussian2DKernel(stddev, **kwargs)


def gaussian_filter(stddev, array):
    """Convolves an array with a Gaussian filter."""

    return astropy.convolution.convolve(
        array, astropy.convolution.Gaussian2DKernel(stddev))


def fit_gaussian(array):
    """Fits a 2D gaussian to an array of data."""

    shape = array.shape
    xmean, ymean = numpy.array(shape) / 2.

    xx, yy = numpy.mgrid[:shape[0], :shape[1]]

    g_init = astropy.modeling.models.Gaussian2D(amplitude=1., x_mean=xmean, y_mean=ymean,
                                                x_stddev=1., y_stddev=1.)

    f2 = astropy.modeling.fitting.LevMarLSQFitter()

    gg = f2(g_init, xx, yy, array)

    return gg


class CCD(object):
    """A class representing the parameters that define a CCD chip.

    Parameters
    ----------
    shape : `tuple`
        The shape of the image to generate.
    pixel_size : `float`
        The pixel size, in microns. Assumes the pixel is square.
    read_noise : `float`
        The RMS of the read noise, in electrons.
    gain : `float`
        The gain in electrons per ADU.
    name : `str` or `None`
        A string with the name of the CCD chip (e.g., its model or SN).

    """

    def __init__(self, shape, pixel_size, read_noise=1.0, gain=1.0, name=None):

        self.shape = shape
        self.pixel_size = pixel_size
        self.read_noise = read_noise
        self.gain = gain
        self.name = name


class SyntheticImage(object):
    """Creates and image with Gaussian features, bias, and noise.

    Parameters
    ----------
    ccd : `.CCD`
        A `.CCD` object describing the chip that produces this image.
    bias : `float` or `None`
        The bias level of the image. If `None`, no bias level will be added.
    cosmic_p : float
        The p factor for the binomial distribution used to model cosmic rays.
    exp_time : float
        The exposure time, in seconds. Used as a multiplicative factor for the
        log-normal distribution to estimate the total dark current.
    dark_sigma : float
        The sigma of the log-normal dark current distribution.
    dtype : `numpy.dtype`
        The dtype for the `.signal` and `.noise` arrays.

    Attributes
    ----------
    bias : float
        The bias level.
    noise : `numpy.ndarray`
        The array representing the image noise.
    signal : `numpy.ndarray`
        The array representing the image signal.
    sources : `list`
        A list of `~astropy.modeling.functional_models.Gaussian2D` objects that
        have been added to the image.

    """

    def __init__(self, ccd, bias=400., cosmic_p=0.005, exp_time=1., dark_sigma=5.,
                 dtype=numpy.uint16):

        self.ccd = ccd
        self.dtype = dtype
        self.signal = numpy.zeros(self.ccd.shape[::-1], dtype=self.dtype)
        self.noise = numpy.zeros(self.ccd.shape[::-1], dtype=self.dtype)

        # Add a bias level
        self.bias = self.dtype(0.0)
        if bias is not None:
            self.add_bias_level(self.dtype(bias))

        if self.ccd.read_noise is not None:
            self.noise += self.get_read_noise().astype(self.dtype)

        self.sources = []

    @property
    def image(self):
        """Returns the signal plus its associated noise."""

        return self.signal + self.noise

    @property
    def snr(self):
        """Returns the signal to noise ratio for each element in the image."""

        return self.signal.astype(numpy.float32) / self.noise.astype(numpy.float32)

    def get_read_noise(self):
        """Returns an array of read noise assuming a normal distribution.

        Readout noise is returned as an array of floats and must be converted
        to the type of the image.

        """

        read_noise_adu = self.ccd.read_noise / self.ccd.gain
        return numpy.random.normal(scale=read_noise_adu, size=self.image.shape)

    def add_bias_level(self, bias):
        """Adds a bias level to the image.

        A ``read_noise`` noise is added to the ``bias`` value. The attribute
        `~bias` is updated with the bias level added by this method.

        Parameters
        ----------
        bias : float
            The bias level to add.

        """

        self.bias += bias
        self.signal += bias

    def add_sources(self, model, params, sample_box=100):
        """Adds multiple sources to the image following a certain model.

        Parameters
        ----------
        model : `class` or `str`
            A class with a ``__call__`` method that accept two 2D arrays with
            the x and y positions, and returns a 2D array of the same shape
            with the value of the model evaluated at each point.
            The model must be initialised with a list of arguments that
            determine its parameters. If ``model='gaussian'``, the model for
            a `2D Gaussian <astropy.modeling.functional_models.Gassian2D>`
            will be used.
        params : `list` or `numpy.ndarray`
            Either a list of list, or a 2D array. In the former case, each
            item must be a list with the parameters to instantiate the
            ``model``. In the latter case, each row in the array are the
            parameters for the model. For a 2D Gaussian, the parameters must
            be ``(amplitude, x_mean, y_mean, sigma_x, sigma_y)``.
        sample_box : `int` or `None`
            The length of the box used to sample the source. If `None`,
            defaults to the shape of the CCD chip. If ``sample_box`` is not
            `None`, the ``model`` must have attributes ``x_mean`` and
            ``y_mean`` that will be used to determine the centroid of the
            evaluated function.

        Returns
        -------
        models : list
            A list of the models used to add the sources, each one instantiated
            with the arguments defined in ``params``.

        """

        if model == 'gaussian':
            model = astropy.modeling.functional_models.Gaussian2D

        # We use a small grid to make the computing of the source faster.
        if sample_box is None:
            sample_box_x, sample_box_y = self.ccd.shape
        else:
            sample_box_x = sample_box_y = sample_box

        y, x = numpy.mgrid[0:sample_box_y, 0:sample_box_x]

        if isinstance(params, numpy.ndarray):
            assert params.ndim == 2, 'invalid number of dimensions in params'
            params = params.tolist()

        for nn in range(len(params)):
            model_n = model(*params[nn])
            self.add_source(model_n, x, y)

    def add_source(self, model, x, y):
        """Adds a source to the image.

        Evaluates the model at a series of ``(x, y)`` positions and adds the
        results to `.signal`. Poisson noise is added to `.noise`. If the
        coordinate arrays ``x`` and ``y`` are smaller than the shape of the
        `.CCD` chip, the model is evaluated in a box defined by the
        input coordinate arrays centred around
        ``(model.x_mean, model.y_mean)``.

        Parameters
        ----------
        model : callable
            A function or instance with a ``__call__`` method that evaluates
            the models at the positions given by ``x`` and ``y``.
        x : numpy.ndarray
            A 2D array with the x positions on which the model will be
            evaluated. Usually the output of `numpy.mgrid`.
        y : numpy.ndarray
            As ``x`` for the y axis.


        """

        x_sample = numpy.min([x.shape[1], self.ccd.shape[0]])
        y_sample = numpy.min([x.shape[0], self.ccd.shape[1]])

        if x_sample == self.ccd.shape[0] and y_sample == self.ccd.shape[1]:

            source_data = model(x, y)
            noise = source_data - numpy.random.poisson(source_data)

            self.signal += source_data.astype(self.dtype)
            self.noise += noise.astype(self.dtype)

        else:

            x_mean = model.x_mean \
                if not isinstance(model.x_mean, astropy.modeling.Parameter) else model.x_mean.value
            y_mean = model.y_mean \
                if not isinstance(model.x_mean, astropy.modeling.Parameter) else model.y_mean.value

            x_offset = int(x_mean) - int(x_sample / 2) if x_mean > x_sample / 2 else 0
            y_offset = int(y_mean) - int(y_sample / 2) if y_mean > y_sample / 2 else 0

            model.x_mean = x_mean - x_offset
            model.y_mean = y_mean - y_offset

            source_data = model(x, y)

            noise = source_data - numpy.random.poisson(source_data)

            self.signal[y_offset:y_offset + y_sample,
                        x_offset:x_offset + x_sample] += source_data.astype(self.dtype)
            self.noise[y_offset:y_offset + y_sample,
                       x_offset:x_offset + x_sample] += noise.astype(self.dtype)

            model.x_mean += x_offset
            model.y_mean += y_offset

        self.sources.append(model)
