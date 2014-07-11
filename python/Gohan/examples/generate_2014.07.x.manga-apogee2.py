#!/usr/bin/env python
# encoding: utf-8
"""
generate_2014.05.x.manga-apogee2.py

Created by José Sánchez-Gallego on 22 May 2014.
Licensed under a 3-clause BSD license.

Revision history:
    22 May 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Gohan import InputCatalogue, PlateInput, log
from astropy import table
import shutil as sh
import os


fields = table.Table.read(
    './fieldSelection/fieldlist_replaced.dat', format='ascii')
stdCatalogue = './aug14_stds/2014.07.x.manga-apogee2_STA.fits'

platePlansBlob = open('2014.07-platePlans.dat', 'w')
plateRun = '2014.07.x.manga-apogee2'

designID = 8367
for field in fields:

    tileID = int(field['tileID'])
    locationID = int(field['locationID'])

    if tileID in [5744, 5809]:
        log.important('Skipping field tileID={0}'.format(tileID))
        designID += 1
        continue

    sciInput = './targets/MaNGA_targets_extNSA_all_final_tiled_v0.5' + \
        '_{0:04d}.fits'.format(tileID)

    sciCat = InputCatalogue(tileid=None, format='file', type='SCI',
                            file=sciInput,
                            meta={'locationid': locationID, 'tileid': tileID},
                            removeSuperfluous=False, failOnCollision=False)

    sciCat.write('sciCatalogue_{0:04d}.fits'.format(tileID))

    stdCat = InputCatalogue(tileid=tileID, format='file', type='STD',
                            file=stdCatalogue, conversions={},
                            meta={'racen': sciCat.meta['racen'],
                                  'deccen': sciCat.meta['deccen'],
                                  'locationid': sciCat.meta['locationid'],
                                  'tileid': sciCat.meta['tileid']},
                            decollision=sciCat)
    stdCat.write('stdCatalogue_{0:04d}.fits'.format(tileID))

    plateInput = PlateInput(designID, plateRun,
                            [sciCat, stdCat])
    # plateInput.plotIFUs()
    plateInput.write(toRepo=True, platePlans=False)
    platePlansBlob.write(plateInput.getPlatePlans() + '\n')

    designID += 1

log.logToRepo(plateRun)
sh.copy(log.logFilename, os.path.join('./', plateRun + '.log'))
