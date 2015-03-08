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
from Gohan.exceptions import GohanError, GohanUserWarning
import numpy as np
from collections import OrderedDict
from Gohan.helpers.runMaNGAPostDesign import log, config, readPath
from Gohan.helpers.utils import getSampleCatalogue, getPlateHolesSortedPath, \
    getPlateInputData, getPointing, getTargetFix
import warnings
from sdss.utilities.yanny import yanny


def parseCatalogID_1(plateMaNGAID, plateTargetsPath):

    plateTargets = table.Table(yanny(plateTargetsPath, np=True)['PLTTRGT'])
    columns = plateTargets.dtype.names

    neverObserve = table.Table.read(
        readPath(config['plateTargets']['neverobserve']),
        format='ascii.no_header', names=['designid'])

    log.info('Loading catalogue data ...')
    catalogue = fits.FITS(getSampleCatalogue(1))[1]

    conversions = {
        'nsa_redshift': 'z',
        'nsa_mstar': 'mass',
        'nsa_petro_th50': 'petroth50',
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

        plateHolesSortedPath = getPlateHolesSortedPath(plate)
        plateHolesSorted = yanny(plateHolesSortedPath, np=True)
        plateHolesSortedStruct = table.Table(
            plateHolesSorted['STRUCT1'],
            names=[nn.lower()
                   for nn in plateHolesSorted['STRUCT1'].dtype.names])

        targetFixPath = getTargetFix(plate)
        targetFix = None if targetFixPath is None else \
            yanny(targetFixPath, np=True)['OPTARFIX']

        designid = int(plateHolesSorted['designid'])

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
                elif 'nsa_' in col.lower():
                    qCol = col.lower()[4:]
                else:
                    qCol = col

                if col.lower() in newRow:
                    continue

                if qCol == 'mangaid':
                    newRow['mangaid'] = mangaid.strip()
                    continue

                if qCol == 'neverobserve':
                    if designid in neverObserve['designid']:
                        nevObs = 1
                    else:
                        nevObs = 0
                    newRow[col.lower()] = nevObs
                    continue

                if 'nsa' in col.lower():
                    if qCol == 'inclination':
                        newRow[col.lower()] = np.round(
                            np.rad2deg(np.arccos(catTarget['BA90'])), 4)
                        continue
                    elif qCol == 'vdisp':
                        newRow[col.lower()] = -999.
                        continue
                    else:
                        newRow[col.lower()] = catTarget[qCol.upper()]
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

                elif qCol.lower() == 'ifudesignwrongsize':
                    # If ifudesignwrongsize cannot be found, we calculate it
                    # ourselves but raise a warning

                    warnings.warn('calculating ifudesignwrongsize from ' +
                                  'ifudesignsize and ifutargetsize.',
                                  GohanUserWarning)

                    ifudesignsize = plateInputData['MANGAINPUT'][
                        'ifudesignsize']
                    ifutargetsize = plateInputData['MANGAINPUT'][
                        'ifutargetsize']

                    newRow[col.lower()] = 0

                    if ifutargetsize > 0:
                        if ifutargetsize > 127:
                            if ifudesignsize < 127:
                                newRow[col.lower()] = 1
                        else:
                            if ifutargetsize > ifudesignsize:
                                newRow[col.lower()] = 1

                    continue

                else:

                    raise GohanError('field {0} not found for mangaid={1}'
                                     .format(col, mangaid))

            mangaid = mangaid.strip()  # Strip it for comparison with targetfix
            if targetFix is not None:
                targetFix_mangaid = targetFix[targetFix['mangaid'] == mangaid]
                if len(targetFix_mangaid) > 0:
                    log.important('Applying target fix to mangaid={0}'
                                  .format(mangaid))
                    for row in targetFix_mangaid:
                        if row['keyword'] in newRow:
                            oldValue = newRow[row['keyword']]
                            newRow[row['keyword']] = row['value']
                            log.debug('mangaid={0}: {1}={2} -> {3}'.format(
                                      mangaid, row['keyword'], oldValue,
                                      newRow[row['keyword']]))

            newTargets.add_row(newRow)

    return newTargets
