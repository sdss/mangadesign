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
from Gohan.helpers.utils import getPlates, getMangaScience
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


def createPlateMags(input, mode='drillRun',
                    plot=True, toRepo=False, copyPreimaging=False,
                    deletePreimaging=False, debug=False):

    if mode == 'drillRun':
        designIDs = getPlates(input, column='designid')
    else:
        designIDs = [int(input)]

    if debug:
        for handler in log.handlers:
            handler.setLevel('DEBUG')

    for designID in designIDs:

        mangaScience = getMangaScience(designID)

        plateMags = PlateMags(mangaScience)
        plateMags.write(toRepo=toRepo)
        plateMags.plot(toRepo=toRepo)

        if copyPreimaging:
            for plateMagsIFU in plateMags:
                copyPreImaging(designID, plateMagsIFU.preImage)
                if deletePreimaging:
                    os.remove(plateMagsIFU.preImage)
                    os.remove(plateMagsIFU.preimageIRG)

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Creates the plateMags files'
                                     ' for a given drill run or deisgnid')
    parser.add_argument('input', type=str, help='the input value ' +
                        '(drill run or designid)')
    parser.add_argument('-d', '--designid', dest='designid',
                        action='store_true', default=False,
                        help='if set, the input will be considered a designid')
    parser.add_argument('-r', '--repo', dest='toRepo', action='store_true',
                        default=False, help='if set, files will be copied ' +
                        'to mangacore')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
                        default=False, help='if set, a plot of the ' +
                        'preimaging is generated')
    parser.add_argument('-c', '--copy', dest='copy', action='store_true',
                        default=False, help='copies the preimaging data to' +
                        ' config.preImaging')
    parser.add_argument('-m', '--move', dest='move', action='store_true',
                        default=False, help='moves the preimaging data to' +
                        ' config.preImaging')
    parser.add_argument('-b', '--debug', dest='debug', action='store_true',
                        default=False, help='turns on debug mode')

    args = parser.parse_args()

    if args.designid is True:
        mode = 'designid'
    else:
        mode = 'drillRun'

    if args.copy:
        copyPreimaging = True
        deletePreimaging = False
    if args.move:
        copyPreimaging = True
        deletePreimaging = True

    createPlateMags(args.input, mode=mode,
                    plot=args.plot, toRepo=args.toRepo,
                    copyPreimaging=copyPreimaging,
                    deletePreimaging=deletePreimaging,
                    debug=args.debug)
