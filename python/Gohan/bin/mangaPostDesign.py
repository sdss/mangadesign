# !/usr/bin/env python
# encoding: utf-8
"""
mangaPostDesign.py

Created by José Sánchez-Gallego on 18 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    18 Jul 2014 J. Sánchez-Gallego
      Initial version
    05 Apr 2015 J. Sánchez-Gallego
      Major redesign.

"""

from __future__ import division
from __future__ import print_function

from Gohan import log, config, readPath
from Gohan import PlateTargets
from Gohan.utils import utils
from Gohan.exceptions import GohanPostDesignError, GohanPostDesignWarning
from Gohan.utils.yanny import yanny

from astropy import table
import numpy as np

from collections import OrderedDict
import warnings
import os
import shutil
import argparse


nBundles = sum(config['IFUs'].values())


# def addTargets(plateTargetsPath, newPlateTargetsTable, removeFirst=False,
#                overwrite=False):

#     plateTargets = yanny(plateTargetsPath, np=True)
#     plateTargetsTable = table.Table(plateTargets['PLTTRGT'], copy=True)

#     if removeFirst:
#         plateTargetsTable.remove_row(0)

#     header = []
#     for row in plateTargets._contents.split('\n'):
#         if row.strip() == '':
#             header.append('')
#         elif row.strip()[0] == '#':
#             header.append(row)
#         else:
#             break

#     header = '\n'.join(map(str, header)).strip() + '\n\n'  # Strips empty lines

#     nNewLines = 0
#     for newTarget in newPlateTargetsTable:

#         plateid = newTarget['plateid']
#         mangaid = newTarget['mangaid']

#         currentRecords = np.where(
#             (plateTargetsTable['plateid'] == plateid) &
#             (plateTargetsTable['mangaid'] == mangaid))[0]

#         if len(currentRecords) == 1:
#             if overwrite:
#                 plateTargetsTable.remove_row(currentRecords[0])
#             else:
#                 log.debug('mangaid={0} in plateid={1} skipped because already'
#                           ' exists in {2} and overwrite=False'.format(
#                               os.path.basename(mangaid), plateid,
#                               plateTargetsPath))
#                 continue

#         elif len(currentRecords) > 1:
#             raise GohanPostDesignError(
#                 '{0} has more than one row with plateid={1} and '
#                 'mangaid={2}'.format(os.path.basename(plateTargetsPath),
#                                      plateid, mangaid))

#         plateTargetsTable.add_row(newTarget)
#         nNewLines += 1

#     os.remove(plateTargetsPath)
#     plateTargets['PLTTRGT'] = plateTargetsTable
#     plateTargets.write(plateTargetsPath, comments=header)

#     # plateTargets.append({'PLTTRGT': newPlateTargetsTable})
#     log.info('Appended {0} targets to {1}'.format(
#         nNewLines, os.path.basename(plateTargetsPath)))


def runMaNGAPostDesign(plateids, overwrite=False):
    """Runs MaNGA post-design procedure.

    This function accepts a list of plateids and updates the necessary
    plateTargets files for all the mangaids in the list of plates.
    It also copies the appropriate plateHolesSorted files to mangacore.

    The routine is intended to be run once the SDSS-IV plate design process
    is finished and all the necessary files have been generated.

    Parameters
    ----------
    plateids : int or list of ints
        A list of plateids that will be processed. It can also be a single
        plateid, as an integer.
    overwrite : bool, optional
        If True, the routine will overwrite the plateTargets lines for the
        targets and replace the plateHolesSorted files.

    Returns
    -------
    result : dictionary of `PlateTargets` instances
        Returns a dictionary with keys the unique catalgids of the targets in
        the plates, and values the `PlateTargets` instances for those
        catalogids.

    """

    plateids = np.atleast_1d(plateids, dtype=int)

    log.info('Identifying mangaids ... ')

    catPlateMangaID = OrderedDict()
    for plate in plates:

        mangaids = utils.getMaNGAIDs(plate)

        if len(mangaids) != nBundles:
            warnings.warn(
                'number of science targets is not {0} for plate {1}'.format(
                    nBundles, plate), GohanPostDesignWarning)

        # Sorts mangaids by catalogid.
        for mangaid in mangaids:

            catID = int(mangaid.split('-')[0])

            if catID not in catPlateMangaID:
                catPlateMangaID[catID] = OrderedDict()

            if plate not in catPlateMangaID[catID]:
                catPlateMangaID[catID][plate] = []

            catPlateMangaID[catID][plate].append(mangaid)

    log.info('Sorting mangaids by catalogid ... ')

    if len(catPlateMangaID.keys()) == 0:
        log.important('No targets found for that platerun.')
        return OrderedDict()

    returnDict = OrderedDict()

    for catID in catPlateMangaID.keys():

        log.info('Doing catalogid={0} ...'.format(catID))

        # parsingFunction = utils.getParsingFunction(catID)

        # if parsingFunction is None:
        #     nSkipped = np.sum([len(mIDs)
        #                        for mIDs in catPlateMangaID[catID].values()])
        #     warnings.warn('no parsing function for catID={0}. '
        #                   'Skipping {1} mangaids'.format(catID, nSkipped),
        #                   GohanUserWarning)
        #     continue

        # template = False
        # plateTargetsPath = utils.getPlateTargetsPath(catID)

        # if not os.path.exists(plateTargetsPath):
        #     log.important('no plateTargets-{0}.par found. Using a template.'
        #                   .format(catID))

        #     templatePath = utils.getPlateTargetsTemplate(catID)

        #     if not os.path.exists(templatePath):
        #         warnings.warn('template cannot be found. Skipping catalogid.')
        #         continue
        #     else:
        #         shutil.copy(templatePath, plateTargetsPath)
        #         template = True

        # log.info('{0} plates found for catalogid={1}'.format(
        #     len(catPlateMangaID[catID]), catID))

        for plate in catPlateMangaID[catID]:
            plateTargetsTable = PlateTargets(
                plate, mangaids=catPlateMangaID[catID][plate], catalogid=catID)
            plateTargets.write()

    log.info('Copying plateHolesSorted to mangacore ...')

    plateHolesDir = os.path.join(readPath(config['mangacore']),
                                 'platedesign/plateholes/')

    if not os.path.exists(plateHolesDir):
        raise GohanPostDesignError('not plateholes dir found in mangacore')

    for plate in plates:
        plateHolesSortedPath = utils.getPlateHolesSortedPath(plate)
        plateHolesPlatePath = os.path.join(
            plateHolesDir, '{0}XX/'.format('{0:06d}'.format(plate)[0:4]))

        if not os.path.exists(plateHolesPlatePath):
            os.makedirs(plateHolesPlatePath)

        destinationPath = os.path.join(
            plateHolesPlatePath, os.path.basename(plateHolesSortedPath))

        if os.path.exists(destinationPath) and not overwrite:
            warnings.warn('{0} already exists in magacore. Not overwritting.'
                          .format(os.path.basename(plateHolesSortedPath)),
                          GohanPostDesignWarning)
            continue

        if os.path.exists(destinationPath):
            warnings.warn('Overwritting {0} in mangacore.'
                          .format(os.path.basename(plateHolesSortedPath)),
                          GohanPostDesignWarning)

        shutil.copy(plateHolesSortedPath, destinationPath)
        log.info('{0} copied to mangacore'.format(
            os.path.basename(plateHolesSortedPath)))

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Performs MaNGA post-design actions.')
    parser.add_argument('--overwrite', '-o', action='store_true',
                        help='forces overwrite of existing records.')
    parser.add_argument('--plateid', '-p', action='store_true',
                        help='if present, the input values will be considered '
                        'plateids instead of plate runs.')
    parser.add_argument('plateRuns', metavar='plateRuns/plateids',
                        type=str, nargs='+', help='the plate run or plateids '
                        'to process.')

    args = parser.parse_args()

    if args.plateid:
        plateids = [int(plateid) for plateid in args.plateRuns]
        runMaNGAPostDesign(plateids, overwrite=args.overwrite)
    else:
        plateRuns = args.plateRuns
        for plateRun in plateRuns:
            log.info('Plate run: ' + plateRun)
            plates = utils.getFromPlatePlans(plateRun, column='plateid')
            runMaNGAPostDesign(plates, overwrite=args.overwrite)
