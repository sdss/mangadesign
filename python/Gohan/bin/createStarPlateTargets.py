#!/usr/bin/env python
# encoding: utf-8
"""
createStarPlateTargets.py

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

from Gohan import log, config
from Gohan import StarPlateTargets, StandardPlateTargets
from Gohan.utils import utils
import numpy as np

import os
import argparse


nBundles = sum(config['IFUs'].values())


def createStarPlateTargets(plateids, overwrite=False):
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
    log.info('Appended {0} targets to {1}'.format(nAppended,
                                                  starPlateTargetsPath))

    standardPlateTargetsPath, nAppended = standardPlateTargets.write()
    log.info('{0} saved'.format(os.path.basename(standardPlateTargetsPath)))
    log.info('Appended {0} targets to {1}'.format(nAppended,
                                                  standardPlateTargetsPath))

    return starPlateTargets, standardPlateTargets


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Adds targets to starPlateTargets.')
    parser.add_argument('--overwrite', '-o', action='store_true',
                        help='forces overwrite of existing records.')
    parser.add_argument('--plateid', '-p', action='store_true',
                        help='if present, the input values will be considered '
                        'plateids instead of plate runs.')
    parser.add_argument('--all', '-a', action='store_true',
                        help='runs all APOGEE2-MaNGA plate runs.')
    parser.add_argument('plateRuns', metavar='plateRuns/plateids',
                        type=str, nargs='*', help='the plate run or plateids '
                        'to process.')

    args = parser.parse_args()

    if len(args.plateRuns) == 0 and not args.all:
        parser.error('plateRuns/plateids must be specified '
                     'unless all=True.')

    if args.all:
        args.plateid = False

    if args.plateid:
        plateids = [int(plateid) for plateid in args.plateRuns]
        createStarPlateTargets(plateids, overwrite=args.overwrite)
    else:
        if not args.all:
            plateRuns = args.plateRuns
        else:
            plateRuns = utils.getStellarLibraryRuns()

        for plateRun in plateRuns:
            log.info('Plate run: ' + plateRun)
            plates = map(int,
                         utils.getFromPlatePlans(plateRun, column='plateid'))
            createStarPlateTargets(plates, overwrite=args.overwrite)
