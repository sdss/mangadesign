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
from ..runMaNGAPostDesign import log, config, readPath
from ..utils import getPlateHolesSortedPath, getPlateInputData, getPointing, \
    getTargetFix
from sdss.utilities.yanny import yanny


def parseCatalogID_12(plateMaNGAID, plateTargetsPath):

    plateTargets = table.Table(yanny(plateTargetsPath, np=True)['PLTTRGT'])
    columns = plateTargets.dtype.names

    neverObserve = table.Table.read(
        readPath(config['plateTargets']['neverobserve']),
        format='ascii.no_header', names=['designid'])

    catalogue = fits.FITS(
        os.path.join(
            os.environ['MANGASAMPLE_DIR'],
            '12-nsa_v1b_0_0_v2.fits.gz'))[1]

    conversions = {
        'nsa_redshift': 'z',
        'nsa_mstar': 'mass',
        'nsa_petro_th50': 'petroth50',
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

        targetFixPath = getTargetFix(plate)
        targetFix = None if targetFixPath is None else \
            yanny(targetFixPath, np=True)['OPTARFIX']

        designid = int(plateHolesSorted['designid'])

        log.info('Loading catalogue data ...')

        if int(plate) != 7443:
            sample = fits.FITS(
                os.path.join(
                    os.environ['MANGASAMPLE_DIR'],
                    'MaNGA_targets_nsa_ran_tiles3_adap_1.5_' +
                    'bun_2_4_4_2_5_Remaj_primplus_Ng_30000_np_1553' +
                    '_absmag.fits'))[1]
        else:
            sample = fits.FITS(readPath('+etc/targets-7443.fits'))[1]

        for mangaid in mangaids:

            plateInputData = getPlateInputData(mangaid, plateHolesSorted)

            plateHolesSortedTarget = plateHolesSortedStruct[
                plateHolesSortedStruct['mangaid'] == mangaid]

            try:
                sampleTarget = sample[
                    sample.where('MANGAID == \'{0}\''.format(mangaid.strip()))]
            except:
                raise GohanError('mangaid={0} not found in sample'
                                 .format(mangaid))

            IAUName = plateInputData['MANGAINPUT']['iauname']
            IAUName = IAUName if '/' not in IAUName else IAUName.split('/')[-1]
            try:
                catTarget = catalogue[
                    catalogue.where('IAUNAME == \'{0}\''.format(
                        IAUName.strip()))]
            except:
                raise GohanError('mangaid={0} not found in NSA catalogue'
                                 .format(mangaid))

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

                if 'nsa' in col:
                    if qCol == 'vdisp' or qCol == 'zdist':
                        newRow[col.lower()] = -999.
                        continue
                    elif qCol == 'inclination':
                        newRow[col.lower()] = np.round(
                        np.rad2deg(np.arccos(catTarget['BA90'])), 4)
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

                elif qCol.upper() in sampleTarget.dtype.names:
                    newRow[col.lower()] = sampleTarget[qCol.upper()]
                    continue

                elif qCol.lower() in ['ifusizetarget', 'ifudesignwrongsize']:
                    newRow[col.lower()] = -999
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
