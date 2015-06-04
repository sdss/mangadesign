#!/usr/bin/env python
# encoding: utf-8
"""
generate_2015.03.x.apogee2-manga.py

Created by José Sánchez-Gallego on 13 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    13 Mar 2015 J. Sánchez-Gallego
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
from Gohan import PlateInput
from Gohan.utils import selectSkies, getStellarLibraryTargets


toRepo = False

rejectTargets = list(getStellarLibraryTargets()['mangaid'])

# Targets we reject becuase the IFU won't get enough skies
# rejectTargets += ['4-21764', '4-21765', '4-21769', '4-21770']

rejectTargetsStd = []

fields = table.Table.read('./fieldList_APOGEE_May2015.dat',
                          format='ascii.commented_header')

platePlansBlob = open('2015.05.x.apogee2-manga-platePlans.dat', 'w')

plateRun = '2015.05.x.apogee2-manga'

urMinorStd = 'Standards/manga_stds_URMINOR.fit'
urMinorStdTable = table.Table.read(urMinorStd)
urMinorStdTable.add_column(
    table.Column(['5-{0:d}'.format(stdID)
                  for stdID in urMinorStdTable['STD_ID']],
                 name='MANGAID', dtype='S30'))

skyTemplate = './Skies/sky_{0}.fits'

for field in fields:

    designID = int(field['design_id'])
    fieldName = field['fieldName']

    raCen = float(field['RA'])
    decCen = float(field['Dec'])

    if fieldName == 'URMINOR':
        sciCat = table.Table.read(
            './Targets/URMINOR_MaNGAstellib_targets.fits')
        stdCat = urMinorStdTable

    print(_color_text('\ndesignID: {0}, fieldName: {1}'
                      .format(designID, fieldName), 'lightred'))

    mangaScience = PlateInput(designID, 'science', catalogues=sciCat,
                              surveyMode='apogeeLead', plateRun=plateRun,
                              silentOnCollision=True, sort=False,
                              raCen=raCen, decCen=decCen,
                              rejectTargets=rejectTargets)
    mangaScience.write(toRepo=toRepo)
    rejectTargets += mangaScience.getMangaIDs()

    mangaStandard = PlateInput(designID, 'standard',
                               catalogues=stdCat,
                               surveyMode='apogeeLead', plateRun=plateRun,
                               silentOnCollision=True, sort=True,
                               decollidePlateInputs=[mangaScience],
                               raCen=raCen, decCen=decCen,
                               rejectTargets=rejectTargetsStd)
    mangaStandard.write(toRepo=toRepo)

    skyCat = 'selectedSkies_{0}_{1}.fits'.format(fieldName, designID)
    if os.path.exists(skyCat):
        os.remove(skyCat)

    skies = skyTemplate.format(fieldName)
    selectSkies(skies, designID, fieldName, raCen, decCen)

    mangaSky = PlateInput(designID, 'sky', catalogues=skyCat,
                          surveyMode='apogeeLead', plateRun=plateRun,
                          silentOnCollision=True, sort=False,
                          raCen=raCen, decCen=decCen)
    mangaSky.write(toRepo=toRepo)

    plotIFUs([mangaScience, mangaStandard, mangaSky], centre=None,
             filename='ifuPlot_{0}.pdf'.format(designID))

sh.copy(log.logFilename, os.path.join('./', plateRun + '.log'))

if toRepo:
    log.logToRepo(plateRun)
