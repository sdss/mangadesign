#!/usr/bin/env python
# encoding: utf-8
"""
checkPlateTargets.py

Created by José Sánchez-Gallego on 21 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    21 Mar 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Gohan.helpers.utils import getPlateTargetsPath, getDesignID
from Gohan import exceptions
from sdss.utilities.yanny import yanny
from astropy import table
import numpy as np
import os
import argparse


def checkPlateTargets(plateTargetsPath):
    """Checks a plateTargets file and searchs for duplicate targets.

    Parameters
    ----------
    plateTargetsPath : path
        The path of the plateTargets file to check.

    Returns
    -------
    result : `astropy.table.Table`
        A table with the plates that contain duplicate targets.

    """

    plateTargets = yanny(plateTargetsPath, np=True)['PLTTRGT']

    results = table.Table(None, names=['plateid', 'designid',
                                       'duplicate_targets',
                                       'neverobserve',
                                       'duplicate_plateids'],
                          dtype=[int, int, int, int, 'S100'])

    plateids = sorted(np.unique(plateTargets['plateid']))

    for nn, plateid in enumerate(plateids):

        nDuplicates = 0
        duplicatePlateIDs = []
        designID = getDesignID(plateid)
        neverobserve = plateTargets[plateTargets['plateid'] ==
                                    plateid]['neverobserve'][0]

        for target in plateTargets[plateTargets['plateid'] == plateid]:
            for plateid2 in plateids[nn+1:]:
                plateid2Targets = plateTargets[plateTargets['plateid'] ==
                                               plateid2]

                if target['mangaid'] in plateid2Targets['mangaid'].tolist():
                    nDuplicates += 1
                    duplicatePlateIDs.append(str(plateid2))

        duplicatePlateIDs = np.unique(duplicatePlateIDs)

        if nDuplicates > 0:
            results.add_row((plateid, designID, nDuplicates, neverobserve,
                             ', '.join(duplicatePlateIDs)))

    return results


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Checks a plateTargets for duplicate targets.')
    parser.add_argument('catalogid', type=int,
                        help='the catalogid of the plateTargets to check.')

    args = parser.parse_args()

    plateTargetsPath = getPlateTargetsPath(args.catalogid)

    if not os.path.exists(plateTargetsPath):
        raise exceptions.GohanError('not plateTargets file found with path '
                                    '{0}'.format(plateTargetsPath))

    result = checkPlateTargets(plateTargetsPath)
    result.pprint()
