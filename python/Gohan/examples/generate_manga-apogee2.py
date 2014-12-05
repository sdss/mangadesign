#!/usr/bin/env python
# encoding: utf-8
"""
generate_manga-apogee2.py

Created by José Sánchez-Gallego on 22 May 2014.
Licensed under a 3-clause BSD license.

Revision history:
    4 Dec 2014 J. Sánchez-Gallego
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


rejectTargets = []
rejectTargetsStd = []

fieldList = 'fieldListDec2014.dat'
plateRun = '2014.12.x.manga-apogee2'
epoch = '2014-12-29'
toRepo = True
savePlot = True
saveLog = True


fields = table.Table.read(fieldList, format='ascii.commented_header')

platePlansBlob = open('{0}-platePlans.dat'.format(plateRun), 'w')

sciTemplate = './targets/MaNGA_targets_extNSA_tiled_{0:d}.fits'
sciTemplateRA = './targets/MaNGA_targets_extNSA_tiled_{0:d}_ra.fits'
stdTemplate = './standards/manga_stds_{0:d}.fit'


for field in fields:

    mangaTileID = int(field['manga_tileid'])
    locationID = int(field['location_id'])
    designID = int(field['design_id'])

    raCen = float(field['RA'])
    decCen = float(field['Dec'])

    sciCat = sciTemplateRA.format(mangaTileID) \
        if os.path.exists(sciTemplateRA.format(mangaTileID)) else \
        sciTemplate.format(mangaTileID)
    stdCat = stdTemplate.format(mangaTileID)

    print(_color_text('mangaTileID: {0}'.format(mangaTileID), 'lightred'))

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

    platePlans = plateDefinition.getPlatePlans(epoch)
    platePlansBlob.write(platePlans + '\n')

    if savePlot:
        plotIFUs([mangaScience, mangaStandard], centre=None,
                 filename='ifuPlot_{0}.pdf'.format(designID))

if saveLog:
    sh.copy(log.logFilename, os.path.join('./', plateRun + '.log'))
    log.logToRepo(plateRun)
