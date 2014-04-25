#!/usr/bin/env python
# encoding: utf-8
"""
Plot fiber bundles and/or individual fibers on images.
"""
# Created by José Sánchez-Gallego on 18 Apr 2014.
# Licensed under a 3-clause BSD license.

# Revision history:
#     18 Apr 2014 J. Sánchez-Gallego
#       Initial version

from __future__ import division
from __future__ import print_function

import argparse
import os

from matplotlib.collections import EllipseCollection
import matplotlib.cm as cm
from matplotlib import pylab as plt

import numpy as np
from astropy import wcs
from astropy.io import fits

import bundle
import zscale

try:
    import aplpy
    __aplpy__ = True
except:
    __aplpy__ = False

class PlotBundle(bundle.Bundle):
    """Load FITS files and overplot an IFU bundle and/or fibers."""
    def _getImage(self, fitsFile, ext):

        if isinstance(fitsFile, fits.HDUList):
            image = fitsFile[ext]
            fileName = os.path.basename(fitsFile.filename())
        elif isinstance(fitsFile, (fits.PrimaryHDU, fits.ImageHDU)):
            image = fitsFile
            fileName = os.path.basename(fitsFile.fileinfo()['file'].name)
        elif isinstance(fitsFile, basestring):
            image = fits.open(fitsFile)[ext]
            fileName = os.path.basename(fitsFile)
        else:
            raise TypeError('fitsFile is not a fits file or pyfits HDU.')

        return image, fileName

    def plot(self, fitsFile, ext=0, outputFile=None, useAPLpy=False,
             save=True, **kwargs):

        if useAPLpy and not __aplpy__:
            warnings.warn('no APLpy module found. Reverting to matplotlib.')
            useAPLpy = False

        image, fileName = self._getImage(fitsFile, ext)
        if outputFile is None:
            outputFile = fileName + '.pdf'

        wcsImage = wcs.WCS(image.header)
        scale = np.abs(wcsImage.wcs.cd[0, 0]) * 3600
        fibrePix = wcsImage.wcs_world2pix(self.fibres[:, [1, 2]], 0)

        zValues = zscale.zscale(image.data)
        imShowArgs = dict(cmap=cm.Greys_r, vmax=zValues[1], vmin=zValues[0])
        imShowArgs.update(kwargs)

        if not useAPLpy:
            fig, ax = plt.subplots()
            ax.set_xlabel(r'x [pixels]')
            ax.set_ylabel(r'y [pixels]')
            ax.imshow(image.data, origin='lower', **imShowArgs)
            ax.autoscale(False)
            self._plotBundle(ax, fibrePix, scale, **kwargs)
            if save:
                plt.savefig(outputFile, bbox_inches='tight')
            return ax

        else:
            del imShowArgs['cmap']
            aplFig = aplpy.FITSFigure(image)
            aplFig.show_grayscale(**imShowArgs)
            self._plotBundleAPLPy(aplFig, scale, **kwargs)
            if save:
                aplFig.save(outputFile)
            return aplFig

    def _plotBundle(self, ax, fibrePix, scale, **kwargs):
        """Plots the location of each fibre in the bundle on the target."""

        ellipseCollectionArgs = dict(edgecolor='r',
                                     facecolor='None', linewidth=0.7)
        ellipseCollectionArgs.update(kwargs)

        width = 2. / scale

        cc = EllipseCollection(
            width, width, 0.0, units='xy', transOffset=ax.transData,
            offsets=fibrePix, **ellipseCollectionArgs)

        ax.add_collection(cc)

    def _plotBundleAPLPy(self, fig, scale, **kwargs):

        showCirclesArgs = dict(edgecolor='r', facecolor='None', linewidth=0.7)
        showCirclesArgs.update(kwargs)
        radius = 1. / 3600
        fig.show_circles(self.fibres[:, 1], self.fibres[:, 2], radius,
                         **showCirclesArgs)

    def plotHexagon(self, fitsFile, outputFile=None, ext=0,
                    useAPLpy=False, save=True, **kwargs):

        if useAPLpy and not __aplpy__:
            warnings.warn('no APLpy module found. Reverting to matplotlib.')
            useAPLpy = False

        image, fileName = self._getImage(fitsFile, ext)
        if outputFile is None:
            outputFile = fileName + '.pdf'

        wcsImage = wcs.WCS(image.header)
        hexagonPix = wcsImage.wcs_world2pix(self.hexagon, 0)

        # reconnect the last point to the first point.
        hexagon = np.concatenate((self.hexagon, [self.hexagon[0]]), axis=0)
        hexagonPix = np.concatenate((hexagonPix, [hexagonPix[0]]), axis=0)

        zValues = zscale.zscale(image.data)
        imShowArgs = dict(cmap=cm.Greys_r, vmax=zValues[1], vmin=zValues[0])
        imShowArgs.update(kwargs)

        plotArgs = dict(linewidth=0.7, color='r', linestyle='solid')
        plotArgs.update(kwargs)

        if not useAPLpy:
            fig, ax = plt.subplots()
            ax.imshow(image.data, origin='lower', **imShowArgs)
            ax.autoscale(False)
            ax.plot(hexagonPix[:, 0], hexagonPix[:, 1], **plotArgs)
            ax.set_xlim(0.0, ax.get_xlim()[1])
            ax.set_ylim(0.0, ax.get_ylim()[1])
            ax.set_xlabel(r'x [pixels]')
            ax.set_ylabel(r'y [pixels]')
            if save:
                plt.savefig(outputFile, bbox_inches='tight')
            return ax

        else:
            del imShowArgs['cmap']
            aplFig = aplpy.FITSFigure(image)
            aplFig.show_grayscale(**imShowArgs)
            aplFig.show_lines([hexagon.T], **plotArgs)
            if save:
                aplFig.save(outputFile)
            return aplFig


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bundle plotting tool.')
    parser.add_argument('RA', metavar='RA', type=float,
                        help='RA of the centre of the bundle in degrees.')
    parser.add_argument('DEC', metavar='DEC', type=float,
                        help='Dec of the centre of the bundle in degrees.')
    parser.add_argument('--file', '-f', metavar='FILE', type=fits.open,
                        help='the FITS file to plot on.')
    parser.add_argument('--extension', '-e', dest='ext', type=int, default=0,
                        help='the extiension of the FITS file to use')
    parser.add_argument('--simbmap', '-s', dest='simbmap', type=str,
                        default='', help='the simbmap to use. If undefined, ' +
                        'the one in $MANGACORE_DIR will be used.')
    parser.add_argument('--size', '-z', dest='size', type=int,
                        default=127, help='the size of IFU to plot.')
    parser.add_argument('--print', '-p', dest='printFibres',
                        help='prints the RA and Dec of the fibres to stdout.',
                        default=False, action='store_true')
    parser.add_argument('--hexagon', '-x', dest='hexagon',
                        help='plots the hexagonal shape of the bundle.',
                        default=False, action='store_true')
    parser.add_argument('--aplpy', '-a', dest='aplpy',
                        help='uses APLpy if available.',
                        default=False, action='store_true')
    parser.add_argument('--ds9', '-d', action='store_true',
                        help='creates a file with DS9 regions for each fibre.')

    args = parser.parse_args()

    if args.aplpy and not __aplpy__:
        warnings.warn('no APLpy module found. Reverting to matplotlib.')

    bundle = PlotBundle(args.RA, args.DEC, simbmap=args.simbmap, size=args.size)

    if args.file is not None:
        if not args.hexagon:
            bundle.plot(args.file, ext=args.ext, useAPLpy=args.aplpy)
        else:
            bundle.plotHexagon(args.file, ext=args.ext, useAPLpy=args.aplpy)

    if args.ds9:
        bundle.createDS9Regions()

    if args.printFibres:
        if not args.hexagon:
            bundle.printBundle()
        else:
            bundle.printHexagon()
