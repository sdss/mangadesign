#!/usr/bin/env python
# encoding: utf-8
#
# post_design.py
#
# Created by José Sánchez-Gallego on 5 Jul 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import shutil
import os

from collections import OrderedDict

from Gohan import log, config, readPath
from Gohan.exceptions import GohanPostDesignWarning, GohanPostDesignError
from Gohan.PlateTargets import PlateTargets
from Gohan.PlateMags import PlateMags
from Gohan.utils import utils
from Gohan.StarPlateTargets import StarPlateTargets
from Gohan.StandardPlateTargets import StandardPlateTargets

import numpy as np


def post_design_apogee_led(plateids, overwrite=False):
    """Adds targets to starPlateTargets and standardPlateTargets.

    This function accepts a list of plateids and updates starPlateTargets and
    standardPlateTargets for all the mangaids in the list of plates.

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
    result : `starPlateTargets` instance
        Returns the up-to-date starPlateTargets object.

    """

    plateids = np.atleast_1d(plateids)

    log.info('Identifying mangaids ... ')

    starPlateTargets = StarPlateTargets()
    standardPlateTargets = StandardPlateTargets()

    for plateid in plateids:
        log.info('Adding targets for plate_id={0}'.format(plateid))
        starPlateTargets.addTargets(plateid=plateid)
        standardPlateTargets.addTargets(plateid)

    starPlateTargetsPath, nAppended = starPlateTargets.write()
    log.info('{0} saved'.format(os.path.basename(starPlateTargetsPath)))
    log.info('Appended {0} targets to {1}'.format(nAppended, starPlateTargetsPath))

    standardPlateTargetsPath, nAppended = standardPlateTargets.write()
    log.info('{0} saved'.format(os.path.basename(standardPlateTargetsPath)))
    log.info('Appended {0} targets to {1}'.format(nAppended, standardPlateTargetsPath))

    return starPlateTargets, standardPlateTargets


def post_design_manga_led(plateids, overwrite=False, skipPlateHolesSorted=False):
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
    skipPlateHolesSorted : bool, optional
        If True, skips copying the plateHolesSorted files to mangacore.

    Returns
    -------
    result : dictionary of `PlateTargets` instances
        Returns a dictionary with keys the unique catalgids of the targets in
        the plates, and values the `PlateTargets` instances for those
        catalogids.

    """

    plateids = np.atleast_1d(plateids)

    nBundles = sum(config['IFUs'].values())

    log.info('Identifying mangaids ... ')

    # Gets managids for each plate and creates a dictionary in which the
    # key is the catalogid and the values are dictionaries of plate-mangaids
    # values
    catPlateMangaID = OrderedDict()
    for plate in plateids:

        mangaids = utils.getMaNGAIDs(plate)

        if len(mangaids) != nBundles:
            log.warning('number of science targets is not {0} for plate {1}'
                        .format(nBundles, plate), GohanPostDesignWarning)

        # Sorts mangaids by catalogid.
        for mangaid in mangaids:

            catID = int(mangaid.split('-')[0])

            if catID not in catPlateMangaID:
                catPlateMangaID[catID] = OrderedDict()

            if plate not in catPlateMangaID[catID]:
                catPlateMangaID[catID][plate] = []

            catPlateMangaID[catID][plate].append(mangaid)

    if len(catPlateMangaID.keys()) == 0:
        log.important('No targets found for that list of plates.')
        return OrderedDict()

    # Creates the returned dictionary {catID: plateTargets}
    returnDict = OrderedDict()

    for catID in sorted(catPlateMangaID.keys()):

        if catID == 0:
            continue

        addedRows = 0

        log.info('Doing catalogid={0} ...'.format(catID))

        if catID not in returnDict:
            # If returnDict does not contain a PlateTargets instance for
            # catId creates it.
            returnDict[catID] = PlateTargets(catID)

        plateTargets = returnDict[catID]

        # Adds the targets with catalogid=catID for each plate
        for plateid in catPlateMangaID[catID]:
            log.info('Adding targets for plate_id={0}'.format(plateid))

            plateMangaids = [mangaid.strip() for mangaid in catPlateMangaID[catID][plateid]]

            if len(plateMangaids) > 0:
                newRows = plateTargets.addTargets(plateMangaids,
                                                  plateid=plateid,
                                                  overwrite=overwrite)
                addedRows += len(newRows)

        # Logs some information
        if addedRows > 0:
            plateTargetsPath, nAppended = plateTargets.write()
            log.info('{0} saved'.format(os.path.basename(plateTargetsPath)))
            log.info('Appended {0} targets to {1}'.format(nAppended, plateTargetsPath))
        else:
            log.info('no targets added to {0}'.format(
                     os.path.basename(plateTargets.path)))

    if skipPlateHolesSorted:
        log.warning('skipping copying plateHolesSorted files to mangacore',
                      GohanPostDesignWarning)
        return returnDict

    # Copies plateHolesSorted to mangacore
    log.info('Copying plateHolesSorted to mangacore ...')

    plateHolesDir = os.path.join(readPath(config['mangacore']), 'platedesign/plateholes/')

    if not os.path.exists(plateHolesDir):
        raise GohanPostDesignError('not plateholes dir found in mangacore')

    for plate in plateids:
        plateHolesSortedPath = utils.getPlateHolesSortedPath(plate)
        plateHolesPlatePath = os.path.join(plateHolesDir,
                                           '{0}XX/'.format('{0:06d}'.format(plate)[0:4]))

        if not os.path.exists(plateHolesPlatePath):
            os.makedirs(plateHolesPlatePath)

        destinationPath = os.path.join(
            plateHolesPlatePath, os.path.basename(plateHolesSortedPath))

        if os.path.exists(destinationPath) and not overwrite:
            log.warning('{0} already exists in magacore. Not overwritting.'
                          .format(os.path.basename(plateHolesSortedPath)),
                          GohanPostDesignWarning)
            continue

        if os.path.exists(destinationPath):
            log.warning('Overwritting {0} in mangacore.'
                          .format(os.path.basename(plateHolesSortedPath)), GohanPostDesignWarning)

        shutil.copy(plateHolesSortedPath, destinationPath)
        log.info('{0} copied to mangacore'.format(os.path.basename(plateHolesSortedPath)))

    return returnDict


def create_plateMags(input, mode='drill_run', plot=False, debug=True, overwrite=False):

    if mode == 'drill_run':
        designIDs = utils.getFromPlatePlans(input, column='designid')
    else:
        designIDs = [int(input)]

    if debug:
        for handler in log.handlers:
            handler.setLevel('DEBUG')

    nDesigns = len(designIDs)

    for nn, designID in enumerate(designIDs):

        if mode == 'drill_run':
            log.info('Creating plateMags for designID={0:d} ({1}/{2})'
                     .format(designID, nn + 1, nDesigns))

        mangaScience = utils.getMangaSciencePath(designID)

        plateMags = PlateMags(mangaScience)
        plateMags.write()

        if plot:
            plateMags.plot(overwrite=overwrite)

    return
