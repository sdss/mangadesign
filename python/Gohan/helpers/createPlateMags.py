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


def createPlateMags(drillRun, plot=True, toRepo=False):

    designIDs = getPlates(drillRun, column='designid')

    for designID in designIDs:

        mangaScience = getMangaScience(designID)

        plateMags = PlateMags(mangaScience)
        plateMags.write(toRepo=toRepo)
        plateMags.plot(toRepo=toRepo)

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Creates the plateMags files'
                                     ' for a given drill run')
    parser.add_argument('drillRun', type=str, help='the drill run to use')
    parser.add_argument('-r', '--repo', dest='toRepo', action='store_true',
                        default=False, help='if set, files will be copied ' +
                        'to mangacore')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
                        default=False, help='if set, a plot of the preimagin' +
                        'is generated')

    args = parser.parse_args()

    createPlateMags(args.drillRun, plot=args.plot, toRepo=args.toRepo)
