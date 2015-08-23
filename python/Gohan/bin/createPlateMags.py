#!/usr/bin/env python
# encoding: utf-8
"""
createPlateMags.py

Created by José Sánchez-Gallego on 25 Aug 2014.
Licensed under a 3-clause BSD license.

Revision history:
    25 Aug 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import argparse
from Gohan.utils import getFromPlatePlans, getMangaSciencePath
from Gohan.PlateMags import PlateMags
from Gohan import config, readPath, log
import os
import shutil as sh


def copyPreImaging(designID, preImage):

    path = os.path.join(
        readPath(config['preImaging']),
        'D00{0:s}XX'.format('{0:04d}'.format(designID)[0:2]),
        '{0:04d}'.format(designID))

    if not os.path.exists(path):
        os.makedirs(path)

    sh.copy(preImage, os.path.join(path, preImage))

    return


def createPlateMags(input, mode='drillRun', plot=False, debug=False,
                    force=False):

    if mode == 'drillRun':
        designIDs = getFromPlatePlans(input, column='designid')
    else:
        designIDs = [int(input)]

    if debug:
        for handler in log.handlers:
            handler.setLevel('DEBUG')

    nDesigns = len(designIDs)

    for nn, designID in enumerate(designIDs):

        if mode == 'drillRun':
            log.info('Creating plateMags for designID={0:d} ({1}/{2})'
                     .format(designID, nn+1, nDesigns))

        mangaScience = getMangaSciencePath(designID)

        plateMags = PlateMags(mangaScience)
        plateMags.write()
        plateMags.plot(overwrite=force)

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Creates the plateMags files'
                                     ' for a given drill run or deisgnid')
    parser.add_argument('input', type=str, help='the input value ' +
                        '(drill run or designid)')
    parser.add_argument('-d', '--designid', dest='designid',
                        action='store_true', default=False,
                        help='if set, the input will be considered a designid')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
                        default=False, help='if set, a plot of the ' +
                        'preimaging is generated')
    parser.add_argument('-f', '--forces', dest='force', action='store_true',
                        default=False, help='forces overwrite files')
    parser.add_argument('-b', '--debug', dest='debug', action='store_true',
                        default=False, help='turns on debug mode')

    args = parser.parse_args()

    if args.designid is True:
        mode = 'designid'
    else:
        mode = 'drillRun'

    createPlateMags(args.input, mode=mode, plot=args.plot, debug=args.debug,
                    force=args.force)
