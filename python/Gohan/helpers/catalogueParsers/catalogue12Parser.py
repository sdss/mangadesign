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
import os
import numpy as np
from collections import OrderedDict
from ..runMaNGAPostDesign import log
from .utils import getPlateHolesSortedPath, getPlateInputData, getPointing, \
    getSampleCatalogue
from sdss.utilities.yanny import yanny


def parseCatalogID_12(plateMaNGAID, plateTargetsPath):

    plateTargets = table.Table(yanny(plateTargetsPath, np=True)['PLTTRGT'])
    columns = plateTargets.dtype.names

    log.info('Loading catalogue data ...')
    sample = fits.FITS(
        os.path.join(
            os.environ['MANGASAMPLE_DIR'],
            'MaNGA_targets_nsa_ran_tiles3_adap_1.5_bun_2_4_4_2_5_Remaj_' +
            'primplus_Ng_30000_np_1553_absmag.fits'))[1]

    catalogue = fits.FITS(getSampleCatalogue(1))[1]

    conversions = {
        'redshift': 'z',
        'mstar': 'stellar_mass',
        'petro_th50': 'petroth50',
        'plateid': 'plateId',
        'ifutargetsize': 'ifusizetarget'
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

            sampleTarget = sample[
                sample.where('MANGAID == \'{0}\''.format(mangaid))]
            if len(sampleTarget) == 0:
                raise GohanError('mangaid={0} not found in parent sample'
                                 .format(mangaid))

            IAUName = plateInputData['MANGAINPUT']['iauname']
            catTarget = catalogue[
                catalogue.where('IAUNAME == \'{0}\''.format(IAUName))]

            if len(catTarget) == 0:
                raise GohanError('mangaid={0} not found in catalogue'
                                 .format(mangaid))

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

                elif qCol.upper() in sampleTarget.dtype.names:
                    newRow[col.lower()] = sampleTarget[qCol.upper()]
                    continue

                else:
                    raise GohanError('field {0} not found for mangaid={1}'
                                     .format(col, mangaid))

            newTargets.add_row(newRow)

    return newTargets
