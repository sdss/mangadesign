#!/usr/bin/env python
# encoding: utf-8
#
# create_starplate_targets.py
#
# Created by José Sánchez-Gallego on 5 Jul 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os

from Gohan import log
from Gohan.StarPlateTargets import StarPlateTargets
from Gohan.StandardPlateTargets import StandardPlateTargets

import numpy as np


def create_starplate_targets(plateids, overwrite=False):
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
