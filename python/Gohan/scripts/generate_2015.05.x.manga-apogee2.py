#!/usr/bin/env python
# encoding: utf-8
"""
generate_2015.04.x.manga-apogee2.py

Created by José Sánchez-Gallego on 5 May 2015.
Licensed under a 3-clause BSD license.

Revision history:
    5 May 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
import shutil as sh
from astropy import table
from Gohan import log
from Gohan.core.colourPrint import _color_text
from Gohan.utils.assignIFUDesigns import plotIFUs
from Gohan import PlateInput, PlateDefinition


plateRun = '2015.05.x.manga-apogee2'

mangaTargetVersion = 'v1_2_12'
toRepo = False

rejectTargets = []
rejectTargetsStd = []

fields = table.Table.read('./fieldListMay2015.dat',
                          format='ascii.commented_header')

platePlansBlob = open('{0}-platePlans.dat'.format(plateRun), 'w')

sciTemplate = os.path.join(os.path.expandvars('$MANGASAMPLE_DIR'),
                           'mangawork/target/{0}'.format(mangaTargetVersion),
                           'MaNGA_targets_extNSA_tiled_ancillary_{0:d}.fits')
stdTemplate = './manga_stds_may15/manga_stds_{0:d}.fit'


for field in fields:

    mangaTileID = int(field['manga_tileid'])
    locationID = int(field['location_id'])
    designID = int(field['design_id'])

    raCen = float(field['RA'])
    decCen = float(field['Dec'])

    sciCat = sciTemplate.format(mangaTileID)

    stdCat = stdTemplate.format(mangaTileID)

    print(_color_text('\ndesignID: {0}, mangaTileID: {1}'
                      .format(designID, mangaTileID), 'lightred'))

    mangaScience = PlateInput(designID, 'science', catalogues=sciCat,
                              surveyMode='mangaLead', plateRun=plateRun,
                              silentOnCollision=True, sort=False,
                              raCen=raCen, decCen=decCen,
                              manga_tileid=mangaTileID, locationid=locationID,
                              rejectTargets=rejectTargets)
    mangaScience.write(toRepo=toRepo)

    mangaStandard = PlateInput(designID, 'standard', catalogues=stdCat,
                               surveyMode='mangaLead', plateRun=plateRun,
                               silentOnCollision=True, sort=False,
                               decollidePlateInputs=[mangaScience],
                               raCen=raCen, decCen=decCen,
                               manga_tileid=mangaTileID, locationid=locationID,
                               rejectTargets=rejectTargetsStd)
    mangaStandard.write(toRepo=toRepo)

    plateDefinition = PlateDefinition([mangaScience, mangaStandard])
    plateDefinition.write(toRepo=toRepo)

    platePlans = plateDefinition.getPlatePlans('2015-06-15')
    platePlansBlob.write(platePlans + '\n')

    plotIFUs([mangaScience], centre=None,
             filename='ifuPlot_{0}.pdf'.format(designID))

sh.copy(log.logFilename, os.path.join('./', plateRun + '.log'))

if toRepo:
    log.logToRepo(plateRun)
