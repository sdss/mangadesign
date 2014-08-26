#!/usr/bin/env python
# encoding: utf-8
"""
runMaNGAPostDesign.py

Created by José Sánchez-Gallego on 18 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.utilities.yanny import yanny
from astropy import table
import warnings
import os
from Gohan import log, config, readPath
from Gohan.helpers.catalogueParsers import utils
from Gohan.exceptions import GohanUserWarning
from collections import OrderedDict
from shutil import copy
import numpy as np


def addTargets(plateTargetsPath, newPlateTargetsTable, overwrite=False):

    plateTargets = yanny(plateTargetsPath, np=True)
    plateTargetsTable = table.Table(plateTargets['PLTTRGT'])

    if overwrite is False:
        rowsToDelete = []
        for nn, newTarget in enumerate(newPlateTargetsTable):

            plateid = newTarget['plateid']
            mangaid = newTarget['mangaid']
            designid = newTarget['designid']

            tmpRow = plateTargetsTable[
                (plateTargetsTable['plateid'] == plateid) &
                (plateTargetsTable['designid'] == designid) &
                (plateTargetsTable['mangaid'] == mangaid)]

            if len(tmpRow) > 0:
                warnings.warn('mangaid={0} is already present in plateTargets'
                              'with the same plateid and designid. Skipping'
                              .format(mangaid))
                rowsToDelete.append(nn)

        newPlateTargetsTable.remove_rows(rowsToDelete)

    plateTargets.append({'PLTTRGT': newPlateTargetsTable})
    log.info('Appended {0} targets to {1}'.format(
        len(newPlateTargetsTable), os.path.basename(plateTargetsPath)))


def runMaNGAPostDesign(plateRun, overwrite=False):
    """Accepts a drillrun and performs a series of actions. This function
    is intended to be run after the SDSS-IV plate design process has been
    finished and the $PLATELIST/plates directory has been populated."""

    log.info('Identifying mangaids ... ')

    plates = utils.getPlates(plateRun)
    nBundles = sum(config['IFUs'].values())

    catPlateMangaID = OrderedDict()
    for plate in plates:

        newMaNGAids = utils.getMaNGAIDs(plate)

        if len(newMaNGAids) != nBundles:
            warnings.warn(
                'number of science targets is not {0} for plate {1}'.format(
                    nBundles, plate), GohanUserWarning)

        for mangaid in newMaNGAids:

            catID = int(mangaid.split('-')[0])

            if catID in catPlateMangaID:
                pass
            else:
                catPlateMangaID[catID] = OrderedDict()

            if plate in catPlateMangaID[catID]:
                pass
            else:
                catPlateMangaID[catID][plate] = []

            catPlateMangaID[catID][plate].append(mangaid)

    log.info('Sorting mangaids by catalogid ... ')

    for catID in catPlateMangaID.keys():

        log.info('Doing catalogid={0} now ...'.format(catID))

        parsingFunction = utils.getParsingFunction(catID)

        if parsingFunction is None:
            nSkipped = np.sum([len(mIDs)
                               for mIDs in catPlateMangaID[catID].values()])
            warnings.warn('no parsing function for catID={0}. '
                          'Skipping {1} mangaids'.format(catID, nSkipped),
                          GohanUserWarning)
            continue

        plateTargetsPath = utils.getPlateTargetsPath(catID)
        if plateTargetsPath is None:
            log.important('not plateTargets-{0}.par found. A template '
                          'is necessary. Skipping catID={0}}.'.format(catID))
            continue

        plateTargetsTable = parsingFunction(
            catPlateMangaID[catID], plateTargetsPath)

        addTargets(plateTargetsPath, plateTargetsTable, overwrite=False)

    log.info('Copying plateHolesSorted to mangacore ...')

    for plate in plates:
        plateHolesSortedPath = utils.getPlateHolesSortedPath(plate)
        mangacorePath = os.path.join(
            readPath(config['mangacore']),
            'platedesign/plateholes/{0}XX/'.format(
                '{0:06d}'.format(plate)[0:4]))

        if not os.path.exists(mangacorePath):
            os.makedirs(mangacorePath)

        destinationPath = os.path.join(
            mangacorePath, os.path.basename(plateHolesSortedPath))

        if os.path.exists(destinationPath) and not overwrite:
            warnings.warn('{0} already exists in magacore. Not overwritting.'
                          .format(os.path.basename(plateHolesSortedPath)),
                          GohanUserWarning)
            continue

        if os.path.exists(destinationPath):
            warnings.warn('Overwritting {0} in mangacore.'
                          .format(os.path.basename(plateHolesSortedPath)),
                          GohanUserWarning)

        copy(plateHolesSortedPath, destinationPath)
        log.info('{0} copied to mangacore'.format(
            os.path.basename(plateHolesSortedPath)))

    return


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(
        description='Performs MaNGA post-design actions.')
    parser.add_argument('--force', '-f', action='store_true',
                        help='forces overwrite of existing records.')
    parser.add_argument('plateRun', help='the plate run to process.')

    args = parser.parse_args()

    runMaNGAPostDesign(args.plateRun, overwrite=args.force)
