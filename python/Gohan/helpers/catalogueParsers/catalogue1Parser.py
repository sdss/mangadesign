#!/usr/bin/env python
# encoding: utf-8
"""
catalogue1Parser.py

Created by José Sánchez-Gallego on 18 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from astropy import table
import fitsio as fits
from Gohan.exceptions import GohanError
import numpy as np
from collections import OrderedDict
from ..runMaNGAPostDesign import log
from .utils import getSampleCatalogue, getPlateHolesSortedPath, \
    getPlateInputData, getPointing
from sdss.utilities.yanny import yanny


def parseCatalogID_1(plateMaNGAID, plateTargetsPath):

    plateTargets = table.Table(yanny(plateTargetsPath, np=True)['PLTTRGT'])
    columns = plateTargets.dtype.names

    log.info('Loading catalogue data ...')
    catalogue = fits.FITS(getSampleCatalogue(1))[1]

    conversions = {
        'redshift': 'z',
        'mstar': 'stellar_mass',
        'petro_th50': 'petroth50',
        'plateid': 'plateId'
    }

    newTargets = table.Table(
        None,
        names=plateTargets.colnames,
        dtype=[plateTargets.dtype[ii]
               for ii in range(len(plateTargets.dtype))])

    for plate in plateMaNGAID:

        log.info('Doing plateid={0} ...'.format(plate))

        mangaids = np.atleast_1d(plateMaNGAID[plate])

        log.debug(
            'Loading plateHolesSorted for plateid={0} ...'.format(plate))
        plateHolesSortedPath = getPlateHolesSortedPath(plate)
        plateHolesSorted = yanny(plateHolesSortedPath, np=True)
        plateHolesSortedStruct = table.Table(
            plateHolesSorted['STRUCT1'],
            names=[nn.lower()
                   for nn in plateHolesSorted['STRUCT1'].dtype.names])

        for mangaid in mangaids:

            plateInputData = getPlateInputData(mangaid, plateHolesSorted)

            plateHolesSortedTarget = plateHolesSortedStruct[
                plateHolesSortedStruct['mangaid'] == mangaid]

            catTarget = catalogue[int(mangaid.split('-')[1])]

            newRow = OrderedDict()

            pointingDict = getPointing(plateInputData['MANGAINPUT'])
            newRow.update(pointingDict)

            for col in columns:

                if col.lower() in conversions:
                    qCol = conversions[col.lower()]
                else:
                    qCol = col

                if col.lower() in newRow:
                    continue

                if qCol == 'mangaid':
                    newRow['mangaid'] = mangaid
                    continue

                if qCol == 'inclination':
                    newRow[col.lower()] = np.round(
                        np.rad2deg(np.arccos(catTarget['BA90'])), 4)
                    continue

                if qCol == 'vdisp':
                    newRow[col.lower()] = -999.
                    continue

                if qCol.lower() in plateHolesSortedTarget.colnames:
                    newRow[col.lower()] = \
                        plateHolesSortedTarget[qCol.lower()][0]
                    continue

                elif qCol in plateHolesSorted.keys():
                    newRow[col.lower()] = plateHolesSorted[qCol]

                elif qCol.lower() in plateInputData['MANGAINPUT'].colnames:
                    newRow[col.lower()] = plateInputData[
                        'MANGAINPUT'][qCol.lower()]
                    continue

                elif qCol.lower() in plateInputData.keys():
                    newRow[col.lower()] = plateInputData[qCol.lower()]
                    continue

                elif qCol.upper() in catTarget.dtype.names:
                    newRow[col.lower()] = catTarget[qCol.upper()]
                    continue

                else:

                    raise GohanError('field {0} not found for mangaid={1}'
                                     .format(col, mangaid))

            newTargets.add_row(newRow)

    return newTargets
