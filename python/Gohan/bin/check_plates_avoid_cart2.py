#!/usr/bin/env python
# encoding: utf-8
#
# check_plates_avoid_cart2.py
#
# Created by José Sánchez-Gallego on 1 Jun 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import astropy.table as table
import os
import progressbar

from Totoro import fromPlateID
from Totoro.utils.utils import avoid_cart_2, get_closest_holes

from Gohan.utils.utils import getAllMaNGAPlates, getStellarLibraryRuns, getFromPlatePlans

platelistDir = os.path.realpath(os.environ['PLATELIST_DIR'])
platePlans = os.path.join(platelistDir, 'platePlans.par')

plates = getAllMaNGAPlates()

mastarPlates = []
for plateRun in getStellarLibraryRuns():
    mastarPlates += getFromPlatePlans(plateRun, column='plateid')

plates += mastarPlates

print('Analysing {0} plates'.format(len(plates)))

tt = table.Table(None, names=['plate', 'avoid_cart_2', 'min_separation',
                              'holeType_1', 'holeType_2'],
                 dtype=[int, bool, float, 'S20', 'S20'])

bar = progressbar.ProgressBar()

for ii in bar(range(len(plates))):

    plate = plates[ii]

    try:
        totoro_plate = fromPlateID(plate)
        if totoro_plate.isComplete:
            continue

        avoid = avoid_cart_2(plate)
        hole_data = get_closest_holes(plate)

        tt.add_row((plate, avoid, hole_data[0], hole_data[-2], hole_data[-1]))
    except Exception as ee:
        print('Found a problem processing plate {0}: {1}'.format(plate, ee))

tt.write('avoid_cart_2.dat', format='ascii.fixed_width', delimiter='|')
