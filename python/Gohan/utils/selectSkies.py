#!/usr/bin/env python
# encoding: utf-8
"""
selectSkies.py

Created by José Sánchez-Gallego on 24 Apr 2015.
Licensed under a 3-clause BSD license.

Revision history:
    24 Apr 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division, print_function

import fnmatch
import os

import numpy as np
from astropy import table
from astropy.coordinates import SkyCoord

from Gohan import config, log, readPath
from Gohan.exceptions import GohanError, GohanWarning

from .sortTargets import sortTargets
from .utils import getPlateDefinition
from .yanny import yanny


nSkies = {127: 8, 91: 6, 61: 4, 37: 2, 19: 2, 7: 1}


def getIFUSkies(data, coords, ra, dec, skyPatrolRadius=14 / 60.,
                minNeightborDist=4, limitTo=40, plot=False):
    """Selects skies near the targets defined in the list of `mangaInputs`
    PlateInput objects."""

    # Selects only skies within the FOV
    centre = SkyCoord(ra=ra, dec=dec, unit='deg')
    separationCentre = coords.separation(centre).deg

    ifuSkies = data[np.where((separationCentre < skyPatrolRadius) &
                             (data['neighbor_dist'] > minNeightborDist))]

    ifuSkyCoords = np.zeros((len(ifuSkies), 2))
    ifuSkyCoords[:, 0] = ifuSkies['ra']
    ifuSkyCoords[:, 1] = ifuSkies['dec']

    ra_com, dec_com = ifuSkyCoords.mean(axis=0)

    try:
        validTargets, assigned = sortTargets(ifuSkyCoords, centre=(ra_com, dec_com),
                                             limitTo=limitTo,
                                             radius=skyPatrolRadius, plot=plot)
    except Exception as ee:
        log.warning('Error: {0}'.format(ee), GohanWarning)
        return False

    return validTargets


def getInfoFromAPOGEE(designID):
    """Reads the APOGEE plateInput files and returns target and design
    information."""

    inputsPath = os.path.join(readPath(config['platelist']), 'inputs')

    plateDefinitionPath = getPlateDefinition(designID)

    plateDefinition = yanny(plateDefinitionPath, np=True)
    nInputs = int(plateDefinition['nInputs'])

    scienceInput = None
    stdInput = None

    for nn in range(nInputs):
        inp = plateDefinition['plateInput{0}'.format(nn + 1)]
        if fnmatch.fnmatch(inp, '*plateInput*_SCI_*.par'):
            scienceInput = os.path.join(inputsPath, inp)
        elif fnmatch.fnmatch(inp, '*plateInput*_STA_*.par'):
            stdInput = os.path.join(inputsPath, inp)
        # elif fnmatch.fnmatch(inp, '*plateInput*_SKY_*.par'):
        #     skyInput = os.path.join(inputsPath, inp)

    assert scienceInput is not None and stdInput is not None, \
        'science or standard input cannot be found in plateDefinition.'
    assert map(os.path.exists, [scienceInput, stdInput]), \
        'one or more of the plateInput paths in plateDefinition does ' + \
        'not exist.'

    sciData = yanny(scienceInput, np=True)
    stdData = yanny(stdInput, np=True)
    # skyData = yanny(skyInput, np=True)

    sciCoords = sciData['APOGEEINPUT1'][['ra', 'dec']]
    stdCoords = stdData['APOGEEINPUT2'][['ra', 'dec']]
    # skyCoords = skyData['APOGEEINPUT3'][['ra', 'dec']]

    apogeeCoords = np.concatenate((sciCoords, stdCoords))  # , skyCoords))
    apogeeCoords = np.array(apogeeCoords.tolist())

    return apogeeCoords


def decollide(aa, bb, distance=config['decollision']['targetAvoid']):
    """Decollides targets."""

    distance = 117 / 3600.
    coordsA = SkyCoord(aa[:, 0], aa[:, 1], unit='deg')
    coordsB = SkyCoord(bb[:, 0], bb[:, 1], unit='deg')

    closest, sep, sep3D = coordsA.match_to_catalog_sky(coordsB, nthneighbor=1)
    coordsA = coordsA[sep.deg > distance]

    if len(coordsA) < 3:
        return False

    closest, sep, sep3D = coordsA.match_to_catalog_sky(coordsA, nthneighbor=2)
    coordsA = coordsA[sep.deg > distance]

    validCoords = np.zeros((len(coordsA), 2))
    validCoords[:, 0] = coordsA.ra.deg
    validCoords[:, 1] = coordsA.dec.deg

    return validCoords


def selectSkies(skyCat, designID, fieldName, raCen, decCen,
                use_apogee=True, raise_error=True, limitTo=40, FOV=None):
    """Writes a list of skies for each one of the IFUs in a design."""

    allSkies = table.Table.read(skyCat)

    for col in allSkies.colnames:
        if col != col.lower():
            allSkies.rename_column(col, col.lower())

    skyCoords = SkyCoord(allSkies['ra'], allSkies['dec'], unit='deg')
    plateCentre = SkyCoord(ra=raCen, dec=decCen, unit='deg')
    separationCentre = skyCoords.separation(plateCentre).deg

    FOV = config['decollision']['FOV'] if FOV is None else FOV

    validSkies = allSkies[
        np.where((separationCentre < FOV) &
                 (separationCentre > config['decollision']['centreAvoid']))]
    validSkies = validSkies[validSkies['neighbor_dist'] <= 60]
    validSkies.sort('neighbor_dist')
    validSkies.reverse()
    skyCoords = SkyCoord(validSkies['ra'], validSkies['dec'], unit='deg')

    mangaScience = table.Table(
        yanny('./mangaScience_{0}_{1}.par'.format(fieldName, designID),
              np=True)['MANGAINPUT'])

    mangaStandard = table.Table(
        yanny('./mangaStandard_{0}_{1}.par'.format(fieldName, designID),
              np=True)['MANGAINPUT'])

    mangaTargets = table.vstack([mangaScience[['ra', 'dec', 'mangaid',
                                               'ifudesign']],
                                 mangaStandard[['ra', 'dec', 'mangaid',
                                                'ifudesign']]])

    mangaTargets.sort('ifudesign', reverse=True)

    mangaCoords = np.zeros((len(mangaTargets), 2))
    mangaCoords[:, 0] = mangaTargets['ra']
    mangaCoords[:, 1] = mangaTargets['dec']

    if use_apogee:
        apogeeCoords = getInfoFromAPOGEE(designID)
        coords = np.concatenate((mangaCoords, apogeeCoords), axis=0)
    else:
        coords = mangaCoords

    skyTable = table.Table(
        None, names=['MANGAID', 'RA', 'DEC', 'MANGA_TARGET2', 'IFUID'],
        dtype=['S20', float, float, int, int])

    mangaID = 1
    for target in mangaTargets:
        ifuDesign = int(target['ifudesign'])

        candidateSkies = getIFUSkies(validSkies, skyCoords,
                                     target['ra'], target['dec'],
                                     limitTo=limitTo, plot=False)

        if candidateSkies is False:
            err_msg = 'target {0} ({1}) has not enough sky candidates'.format(target['mangaid'],
                                                                              ifuDesign)
            return False, err_msg, target['mangaid'], ifuDesign

        decollidedSkies = decollide(candidateSkies, coords)

        if decollidedSkies is False:
            err_msg = 'target {0} ({1}) has not enough skies'.format(target['mangaid'], ifuDesign)
            if raise_error:
                raise GohanError(err_msg)
            else:
                return False, err_msg, target['mangaid'], ifuDesign

        ifuDesignSize = int(str(ifuDesign)[0:-2])

        if nSkies[ifuDesignSize] > len(decollidedSkies):
            err_msg = 'target {0} ({1}) has not enough skies ({2} when {3} are needed)'.format(
                target['mangaid'], ifuDesign, len(decollidedSkies), nSkies[ifuDesignSize])

            if raise_error:
                raise GohanError(err_msg)
            else:
                return False, err_msg, target['mangaid'], ifuDesign

        ifuSkies = decollidedSkies[:nSkies[ifuDesignSize], :]

        for ra, dec in ifuSkies:
            skyTable.add_row(('0-{0}'.format(mangaID), ra, dec, 2, ifuDesign))
            coords = np.append(coords, np.array([[ra, dec]]), axis=0)
            mangaID += 1

    if len(skyTable) < 92:
        raise GohanError('not enough targets ({0})'.format(len(skyTable)))

    # fileName = 'selectedSkies_{0}_{1}.fits'.format(fieldName, designID)
    # skyTable.write(fileName)

    return skyTable
