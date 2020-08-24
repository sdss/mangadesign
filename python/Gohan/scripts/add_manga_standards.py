#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-29
# @Filename: add_manga_standards.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os

from Gohan import log
from Gohan.StandardPlateTargets import StandardPlateTargets
from Gohan.utils import getAllMaNGAPlates


def add_manga_standards():
    """Adds the standards from galaxy plates to standardPlateTargets."""

    plates = getAllMaNGAPlates()
    standards = StandardPlateTargets()

    for plate in plates:
        log.info('Adding targets for plate_id={0}'.format(plate))
        try:
            standards.addTargets(plate)
        except BaseException as ee:
            log.error(f'Cannot process plate {plate}: {ee}')
            continue

    standardPlateTargetsPath, nAppended = standards.write()
    log.info('{0} saved'.format(os.path.basename(standardPlateTargetsPath)))
    log.info('Appended {0} targets to {1}'.format(nAppended, standardPlateTargetsPath))


if __name__ == '__main__':
    add_manga_standards()
