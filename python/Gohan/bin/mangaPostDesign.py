#!/usr/bin/env python
# encoding: utf-8
"""
mangaPostDesign.py

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
from Gohan.helpers import utils
from Gohan.exceptions import GohanUserWarning, GohanError
from collections import OrderedDict
import shutil
import numpy as np


def addTargets(plateTargetsPath, newPlateTargetsTable, removeFirst=False,
               overwrite=False):

    plateTargets = yanny(plateTargetsPath, np=True)
    plateTargetsTable = table.Table(plateTargets['PLTTRGT'], copy=True)

    if removeFirst:
        plateTargetsTable.remove_row(0)

    header = []
    for row in plateTargets._contents.split('\n'):
        if row.strip() == '':
            header.append('')
        elif row.strip()[0] == '#':
            header.append(row)
        else:
            break

    header = '\n'.join(map(str, header)).strip() + '\n\n'  # Strips empty lines

    nNewLines = 0
    for newTarget in newPlateTargetsTable:

        plateid = newTarget['plateid']
        mangaid = newTarget['mangaid']

        currentRecords = np.where(
            (plateTargetsTable['plateid'] == plateid) &
            (plateTargetsTable['mangaid'] == mangaid))[0]

        if len(currentRecords) == 1:
            if overwrite:
                plateTargetsTable.remove_row(currentRecords[0])
            else:
                log.debug('mangaid={0} in plateid={1} skipped because already'
                          ' exists in {2} and overwrite=False'.format(
                              os.path.basename(mangaid), plateid,
                              plateTargetsPath))
                continue

        elif len(currentRecords) > 1:
            raise GohanError('{0} has more than one row with plateid={1} and '
                             'mangaid={2}'.format(
                                 os.path.basename(plateTargetsPath),
                                 plateid, mangaid))

        plateTargetsTable.add_row(newTarget)
        nNewLines += 1

    os.remove(plateTargetsPath)
    plateTargets['PLTTRGT'] = plateTargetsTable
    plateTargets.write(plateTargetsPath, comments=header)

    # plateTargets.append({'PLTTRGT': newPlateTargetsTable})
    log.info('Appended {0} targets to {1}'.format(
        nNewLines, os.path.basename(plateTargetsPath)))


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

    if len(catPlateMangaID.keys()) == 0:
        log.important('No targets found for that platerun.')
        return

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

        template = False
        plateTargetsPath = utils.getPlateTargetsPath(catID)

        if not os.path.exists(plateTargetsPath):
            log.important('no plateTargets-{0}.par found. Using a template.'
                          .format(catID))

            templatePath = os.path.join(
                readPath(config['plateTargets']['templates']),
                'plateTargets-{0}.template'.format(catID))

            if not os.path.exists(templatePath):
                warnings.warn('template cannot be found. Skipping.')
                continue
            else:
                shutil.copy(templatePath, plateTargetsPath)
                template = True

        log.info('{0} plates found for catalogid={1}'.format(
            len(catPlateMangaID[catID]), catID))

        plateTargetsTable = parsingFunction(
            catPlateMangaID[catID], plateTargetsPath)

        addTargets(plateTargetsPath, plateTargetsTable,
                   removeFirst=template, overwrite=overwrite)

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

        shutil.copy(plateHolesSortedPath, destinationPath)
        log.info('{0} copied to mangacore'.format(
            os.path.basename(plateHolesSortedPath)))

    return


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(
        description='Performs MaNGA post-design actions.')
    parser.add_argument('--overwrite', '-o', action='store_true',
                        help='forces overwrite of existing records.')
    parser.add_argument('plateRuns', type=str, nargs='+',
                        help='the plate run to process.')

    args = parser.parse_args()

    for plateRun in args.plateRuns:
        log.info('Plate run: ' + plateRun)
        runMaNGAPostDesign(plateRun, overwrite=args.overwrite)
