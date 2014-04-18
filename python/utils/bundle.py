#!/usr/bin/env python
# encoding: utf-8
"""
plotBundle.py

Created by José Sánchez-Gallego on 18 Apr 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
from astropy.io import fits
from astropy import wcs
from astropy import table
from yanny import read_yanny
import numpy as np
import argparse
from matplotlib.collections import EllipseCollection
import matplotlib.cm as cm
from zscale import zscale


from matplotlib import pylab as plt
try:
    import aplpy
    __use_aplpy__ = True
except:
    __use_aplpy__ = False


class Bundle(object):

    def __init__(self, RA, Dec, size=127, simbmap=None, **kwargs):

        self.RA = float(RA)
        self.dec = float(Dec)

        if simbmap is None or simbmap == '':
            self.simbMapFile = os.path.join(
                os.path.expanduser(os.environ['MANGACORE_DIR']),
                'metrology/fiducial/manga_simbmap_127.par')

        if not os.path.exists(self.simbMapFile):
            raise NameError('not simbmap found.')

        try:
            self.simbMap = table.Table(read_yanny(self.simbMapFile)['SIMBMAP'])
        except:
            raise IOError('cannot read simbmap.')

        self.fibres = self.getFibreCoordinates()
        self.fibres = self.fibres[self.fibres[:, 0] <= size]

    def getFibreCoordinates(self):
        """Returns the RA, Dec coordinates for each fibre."""

        fibreWorld = np.zeros((len(self.simbMap), 3), float)

        raOffDeg = self.simbMap['raoff'] / 3600. / np.cos(self.dec *
                                                          np.pi / 180.)
        decOffDeg = self.simbMap['decoff'] / 3600.

        fibreWorld[:, 0] = self.simbMap['fnumdesign']
        fibreWorld[:, 1] = self.RA + raOffDeg
        fibreWorld[:, 2] = self.dec + decOffDeg

        return fibreWorld

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

    def plot(self, fitsFile, ext=0, outputFile=None, **kwargs):

        image, fileName = self._getImage(fitsFile, ext)
        if outputFile is None:
            outputFile = fileName + '.pdf'

        wcsImage = wcs.WCS(image.header)
        scale = np.abs(wcsImage.wcs.cd[0, 0]) * 3600
        fibrePix = wcsImage.wcs_world2pix(self.fibres[:, [1, 2]], 0)

        zValues = zscale(image.data)
        imShowArgs = dict(cmap=cm.Greys_r, vmax=zValues[1], vmin=zValues[0])
        imShowArgs.update(kwargs)

        if not __use_aplpy__:
            fig, ax = plt.subplots()
            ax.set_xlabel(r'x [pixels]')
            ax.set_ylabel(r'y [pixels]')
            ax.imshow(image.data, origin='lower', **imShowArgs)
            self._plotBundle(ax, fibrePix, scale, **kwargs)
            plt.savefig(outputFile)
            return ax

        else:
            del imShowArgs['cmap']
            aplFig = aplpy.FITSFigure(image)
            aplFig.show_grayscale(**imShowArgs)
            self._plotBundleAPLPy(aplFig, scale, **kwargs)
            aplFig.save(outputFile)
            return aplFig

    def _plotBundle(self, ax, fibrePix, scale, **kwargs):
        """Plots the location of each fibre in the bundle on the target."""

        ellipseCollectionArgs = dict(edgecolor='r',
                                     facecolor='None', linewidth=0.7)
        ellipseCollectionArgs.update(kwargs)

        width = 2. / scale
        xy = [(fibrePix[ii][0], fibrePix[ii][1])
              for ii in range(len(fibrePix[:, 0]))]

        cc = EllipseCollection(
            width, width, 0.0, units='xy', transOffset=ax.transData,
            offsets=xy, **ellipseCollectionArgs)

        ax.add_collection(cc)

    def _plotBundleAPLPy(self, fig, scale, **kwargs):

        showCirclesArgs = dict(edgecolor='r', facecolor='None', linewidth=0.7)
        showCirclesArgs.update(kwargs)
        radius = 1. / 3600
        fig.show_circles(self.fibres[:, 1], self.fibres[:, 2], radius,
                         **showCirclesArgs)

    def createDS9Regions(self, outputFile=None):

        if outputFile is None:
            outputFile = 'bundlePositionsDS9.reg'

        radius = 1. / 3600.

        template = """
global color=green font="helvetica 10 normal roman" wcs=wcs
        """
        template.strip().replace('\n', ' ')

        for fibre in self.fibres:
            template += ('\nfk5;circle({0:.10f},{1:.10f},'
                         '{2:.10f}) #text = {{{3}}}').format(
                fibre[1], fibre[2], radius, int(fibre[0]))

        out = open(outputFile, 'w')
        out.write(template)
        out.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Bundle plotting tool.')
    parser.add_argument('file', metavar='FILE', type=fits.open,
                        help='the FITS file to plot on.')
    parser.add_argument('RA', metavar='RA', type=float,
                        help='RA of the centre of the bundle in degrees.')
    parser.add_argument('DEC', metavar='DEC', type=float,
                        help='Dec of the centre of the bundle in degrees.')
    parser.add_argument('--extension', '-e', dest='ext', type=int, default=0,
                        help='the extiension of the FITS file to use')
    parser.add_argument('--simbmap', '-s', dest='simbmap', type=str,
                        default='', help='the simbmap to use. If undefined, ' +
                        'the one in $MANGACORE_DIR will be used.')
    parser.add_argument('--size', '-z', dest='size', type=int,
                        default=127, help='the size of IFU to plot.')

    args = parser.parse_args()

    bundle = Bundle(args.RA, args.DEC, simbmap=args.simbmap, size=args.size)
    bundle.plot(args.file, ext=args.ext)
    bundle.createDS9Regions()
