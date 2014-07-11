#!/usr/bin/env python
# encoding: utf-8
"""
generate2014.07.x.apogee2-manga.py

Created by José Sánchez-Gallego on 10 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    10 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.utilities.yanny import yanny
import os
import glob
from astropy import table
from astropysics import coords
from Gohan import InputCatalogue, PlateInput
import numpy as np
from Gohan import log


inputAPOGEE = os.path.join(
    os.environ['PLATELIST_DIR'], 'inputs', 'apogee', '2014.07.x.apogee2')
inputAPOGEE2 = os.path.join(
    os.environ['PLATELIST_DIR'], 'inputs', 'apogee', '2014.05.x.apogee2')


def getAPOGEEDesgins():

    files = glob.glob(os.path.join(inputAPOGEE, '*'))

    designs = []
    for file in files:
        design = int(os.path.splitext(file)[0].split('_')[-1])
        if design not in designs:
            designs.append(design)

    designs.sort()

    return designs


def getAPOGEEFiles(design):

    if design != 7977:
        files = glob.glob(
            os.path.join(inputAPOGEE, '*_{0}.par'.format(design)))
    else:
        files = glob.glob(
            os.path.join(inputAPOGEE2, '*_{0}.par'.format(design)))

    if len(files) != 3:
        raise ValueError('more than three files found for design {0}'.format(
                         design))

    files.sort()

    return files


def getCoordinates(inputSci, inputStd):

    inputSci = yanny(inputSci, np=True)
    inputStd = yanny(inputStd, np=True)

    coords = table.vstack(
        [table.Table(inputSci['APOGEEINPUT1'][['ra', 'dec']]),
         table.Table(inputStd['APOGEEINPUT2'][['ra', 'dec']])])

    return [(float(inputSci['raCen']), float(inputSci['decCen']),
             int(inputSci['locationid'])), coords]


def getValid(cat, cen):

    validStars = []
    centreCoords = np.array([coords.ICRSCoordinates(cen[0], cen[1])])
    stellCoords = np.array(
        [coords.ICRSCoordinates(
            cat['RA'][ii], cat['DEC'][ii]) for ii in range(len(cat))])
    validStars = (centreCoords - stellCoords) < coords.AngularSeparation(1.49)

    validCat = cat[validStars]
    validCat.sort('PRIORITY')

    return validCat


def generateDesigns():

    designIDs = [7977]  # + getAPOGEEDesgins()

    stellLibr = table.Table.read(
        'MaNGAstellib_targets_Jul4.tar/' +
        '2014.07.x.apogee2-manga_StellarLibrary.fits')
    sdssCat = table.Table.read(
        'STD_apogee_aug14/2014.07.x.apogee2-manga_STA.fits')
    apassCat = table.Table.read('APASS/2014.07.x.apogee2-manga_APASS.fits')

    for design in designIDs:

        inputs = getAPOGEEFiles(design)
        fieldName = inputs[0].split('_')[-3]

        pairs, collCoordinates = getCoordinates(inputs[0], inputs[2])
        racen = pairs[0]
        deccen = pairs[1]
        locationid = pairs[2]

        stellLibrTargets = getValid(stellLibr, (racen, deccen))
        inputSci = InputCatalogue(
            tileid=None, format='table', input=stellLibrTargets,
            decollision=collCoordinates,
            meta={'racen': racen, 'deccen': deccen, 'locationid': locationid,
                  'fieldname': fieldName},
            warnOnCollision=False)

        for target in inputSci.data:
            mangaid = target['mangaid']
            stellLibr.remove_rows(np.where(stellLibr['MANGAID'] == mangaid))

        newCollCatalogue = table.vstack([collCoordinates,
                                         inputSci.data[['ra', 'dec']]])

        sdssCatValid = getValid(sdssCat, (racen, deccen))
        apassCatValid = getValid(apassCat, (racen, deccen))
        cols = ['MANGAID', 'RA', 'DEC', 'MANGA_TARGET1',
                'PSFMAG', 'PRIORITY', 'EXTINCTION']
        stdCat = table.vstack([sdssCatValid[cols], apassCatValid[cols]])
        stdCat['PRIORITY'] = np.arange(1, len(stdCat)+1)

        inputStd = InputCatalogue(
            input=stdCat, format='table',
            decollision=newCollCatalogue, type='STD',
            meta={'racen': racen, 'deccen': deccen,
                  'locationid': locationid, 'fieldname': fieldName},
            warnOnCollision=False)

        if len(inputSci) < 17:
            log.important('Skipping field {0}: too few targets.'.format(
                          fieldName))
            continue
        elif len(inputStd) < 12:
            log.important('Skipping field {0}: too few standards.'.format(
                          fieldName))
            continue

        plateInput = PlateInput(design, '2014.07.x.apogee2-manga',
                                [inputSci, inputStd],
                                plateType='apogeeLeading',
                                reassignFerrules=True)
        plateInput[1].plotIFUs('plateIFUs_{0}_{1:04d}.pdf'.format(
            fieldName, design))
        plateInput.write(toRepo=True)

    return


if __name__ == '__main__':
    generateDesigns()
