#!/usr/bin/env python
# encoding: utf-8
"""
Identify the locations and sizes of MaNGA IFU bundles and individual fibers.

Created by José Sánchez-Gallego on 18 Apr 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Apr 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

import os
import numpy as np
import warnings

from astropy import table
from astropy.io import ascii

from yanny import read_yanny

class Bundle(object):
    """The location, size, and shape of a MaNGA IFU bundle."""
    def __init__(self, RA, Dec, size=127, simbmap=None, **kwargs):
        """A bundle of a given size at the RA/Dec coordinates."""

        self.RA = float(RA)
        self.dec = float(Dec)
        self.size = size

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

        self._calculateHexagon()

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

    def _calculateHexagon(self):
        """Calculates the vertices of the bundle hexagon."""

        if self.size not in [19, 37, 61, 91, 127]:
            self.hexagon = None
            return

        simbMapSize = self.simbMap[self.simbMap['fnumdesign'] <= self.size]

        middleRow = simbMapSize[simbMapSize['decoff'] == 0]
        topRow = simbMapSize[simbMapSize['decoff'] ==
                             np.max(simbMapSize['decoff'])]
        bottopRow = simbMapSize[simbMapSize['decoff'] ==
                                np.min(simbMapSize['decoff'])]

        vertice0 = middleRow[middleRow['raoff'] == np.max(middleRow['raoff'])]
        vertice3 = middleRow[middleRow['raoff'] == np.min(middleRow['raoff'])]

        vertice1 = topRow[topRow['raoff'] == np.max(topRow['raoff'])]
        vertice2 = topRow[topRow['raoff'] == np.min(topRow['raoff'])]

        vertice5 = bottopRow[bottopRow['raoff'] == np.max(bottopRow['raoff'])]
        vertice4 = bottopRow[bottopRow['raoff'] == np.min(bottopRow['raoff'])]

        hexagonOff = np.array(
            [[vertice0['raoff'][0], vertice0['decoff'][0]],
             [vertice1['raoff'][0], vertice1['decoff'][0]],
             [vertice2['raoff'][0], vertice2['decoff'][0]],
             [vertice3['raoff'][0], vertice3['decoff'][0]],
             [vertice4['raoff'][0], vertice4['decoff'][0]],
             [vertice5['raoff'][0], vertice5['decoff'][0]]])

        raOffDeg = hexagonOff[:, 0] / 3600. / np.cos(hexagonOff[:, 1] *
                                                     np.pi / 180.)
        decOffDeg = hexagonOff[:, 1] / 3600.

        hexagon = hexagonOff.copy()
        hexagon[:, 0] = self.RA + raOffDeg
        hexagon[:, 1] = self.dec + decOffDeg

        self.hexagon = hexagon

        return hexagon


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

        template += '\n'

        out = open(outputFile, 'w')
        out.write(template)
        out.close()

    def printBundle(self):

        tt = table.Table(self.fibres, names=['fnumdesign', 'RA', 'Dec'],
                         dtype=[int, float, float])
        ascii.write(tt, format='fixed_width_two_line',
                    formats={'RA': '{0:.12f}', 'Dec': '{0:.12f}'})

    def printHexagon(self):

        tt = table.Table(self.hexagon, names=['RA', 'Dec'],
                         dtype=[float, float])
        ascii.write(tt, format='fixed_width_two_line',
                    formats={'RA': '{0:.12f}', 'Dec': '{0:.12f}'})

