#!/usr/bin/env python
# encoding: utf-8
"""
autocomplete.py

Created by José Sánchez-Gallego on 10 Nov 2014.
Licensed under a 3-clause BSD license.

Revision history:
    10 Nov 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Gohan import log, readPath, config
from astropy import table
from astropy import coordinates as coo
from Gohan.exceptions import GohanUserWarning
import warnings
import numpy as np


defaultValues = {
    'ifudesign': -999,
    'ifudesignsize': -999,
    'manga_target1': 0,
    'manga_target2': 0,
    'manga_target3': 0,
    'psfmag': [0.0, 0.0, 0.0, 0.0, 0.0],
    'sourcetype': None
}


def autocomplete(targets, targettype, centre, **kwargs):
    """Autocompletes a list of targets using NSA to a number of `nTargets`."""

    if targettype == 'science':
        bundles = config['IFUs'].copy()
    elif targettype == 'standard':
        bundles = config['miniBundles'].copy()
    elif targettype == 'sky':
        bundles = config['skies'].copy()

    nBundles = np.sum(bundles.values())

    if len(targets) >= nBundles:
        log.info('the are at least as many targets as bundles. No need '
                 'to run autocomplete.')
        return targets

    nBundlesToAssign = nBundles - len(targets)
    log.info('autocompleting {0} bundles.'.format(nBundlesToAssign))

    NSATargets = getNSATargets(targets, centre)
    unassignedBundleSizes = getUnassignedBundleSizes(targets, bundles)

    for bundleSize in unassignedBundleSizes:

        target = _getOptimalTarget(targets, NSATargets, bundleSize, centre)

        if target is not None:
            targets = addTarget(targets, target, bundleSize, centre)
        else:
            warnings.warn('no valid replacement found for '
                          'ifudesignsize={0}'.format(int(bundleSize)),
                          GohanUserWarning)

        if len(targets) >= nBundles:
            break

    if len(targets) < nBundles:
        log.important('even after autocomplete, {0} '.format(
            nBundles-len(targets)) + 'bundles are still unassigned.')

    return targets


def _getOptimalTarget(targets, candidateTargets, bundleSize, centre):

    reffField = config['plateMags']['reffField'].upper()
    if reffField in candidateTargets.colnames:
        return _getOptimalTargetSize(targets, candidateTargets, bundleSize,
                                     centre)
    else:
        for target in candidateTargets:
            if _checkNSATarget(targets, target, centre):
                return target
        return None


def _getOptimalTargetSize(targets, candidateTargets, bundleSize, centre):

    if 'IFUDESIGNSIZE' not in candidateTargets.colnames:
        candidateTargets.add_column(
            table.Column(
                [-999] * len(candidateTargets),
                name='IFUDESIGNSIZE', dtype=int))

    reffField = config['plateMags']['reffField'].upper()
    candidateTargets.sort(reffField)
    candidateTargets.reverse()

    if bundleSize == 19:
        minSize = 5 * 2.5
    elif bundleSize == 37:
        minSize = 6 * 2.5
    elif bundleSize == 61:
        minSize = 6 * 2.5
    elif bundleSize == 91:
        minSize = 7 * 2.5
    elif bundleSize == 127:
        minSize = 8 * 2.5
    else:
        minSize = 0.

    for target in candidateTargets:

        if target[reffField] >= minSize:
            if _checkNSATarget(targets, target, centre):
                target['IFUDESIGNSIZE'] = bundleSize
                return target

    # If it comes to this, it means that there are no target with valid
    # size. We fall back to use whatever target that does not collide.
    log.info('no target of right size ({0})'.format(bundleSize) +
             ' can be used to autocomplete. '
             'Trying to find a suitable target.')

    for target in candidateTargets:
            if _checkNSATarget(targets, target, centre):
                return target

    return None


def getUnassignedBundleSizes(targets, bundles):

    for target in targets:
        if target['IFUDESIGNSIZE'] > 0:
            bundles[target['IFUDESIGNSIZE']] -= 1

    unassignedBundleSizes = []
    for size in bundles:
        for nn in range(bundles[size]):
            unassignedBundleSizes.append(size)

    return unassignedBundleSizes


def getNSATargets(targets, centre):

    NSACat = table.Table.read(readPath(config['catalogues']['NSA']))
    NSACat.add_column(
        table.Column(np.arange(len(NSACat)), dtype='int', name='CATIND'))

    raCen, decCen = centre

    coords = coo.SkyCoord(NSACat['RA'], NSACat['DEC'], unit='deg')
    separation = coords.separation(
        coo.SkyCoord(ra=raCen, dec=decCen, unit='deg')).deg

    NSATargets = NSACat[np.where(separation <= config['decollision']['FOV'])]

    NSATargets.add_column(
        table.Column([-999] * len(NSATargets),
                     name='IFUDESIGNSIZE', dtype=int))
    NSATargets.add_column(
        table.Column([-999] * len(NSATargets),
                     name='IFUDESIGN', dtype=int))

    return _getValidNSATargets(targets, NSATargets)


def _getValidNSATargets(targets, NSATargets):
    """Selects targets that are not in the target list already."""

    currentNSATargets = [int(row['MANGAID'].split('-')[1]) for row in targets
                         if row['MANGAID'].split('-')[0] == '1']

    validIdx = [True if target['CATIND'] not in currentNSATargets
                else False for target in NSATargets]

    validNSATargets = NSATargets[np.where(validIdx)]
    validNSATargets.sort('Z')

    mangaIDs = ['1-{0}'.format(row['CATIND']) for row in validNSATargets]

    validNSATargets.add_column(
        table.Column(mangaIDs, name='MANGAID', dtype='S20'))

    return validNSATargets


def _checkNSATarget(targets, target, centre):

    if target['MANGAID'] in targets['MANGAID']:
        return False

    FOV = config['decollision']['FOV']
    centreAvoid = config['decollision']['centreAvoid']
    targetAvoid = config['decollision']['targetAvoid']

    centralPost = coo.SkyCoord(centre[0], centre[1], unit='deg')
    targetCoords = coo.SkyCoord(target['RA'], target['DEC'], unit='deg')

    if targetCoords.separation(centralPost).deg < centreAvoid:
        return False
    if targetCoords.separation(centralPost).deg > FOV:
        return False

    for inputTarget in targets:

        inputCoords = coo.SkyCoord(inputTarget['RA'], inputTarget['DEC'],
                                   unit='deg')

        if inputCoords.separation(targetCoords).deg < targetAvoid:
            return False

    return True


def addTarget(targets, target, bundleSize, centre):

    newRow = []

    for column in targets.colnames:
        if column.upper() == 'IFUDESIGNSIZE':
            newRow.append(bundleSize)
        elif column.upper() == 'PLATERA':
            newRow.append(centre[0])
        elif column.upper() == 'PLATEDEC':
            newRow.append(centre[1])
        elif column.upper() in target.dtype.names:
            newRow.append(target[column.upper()])
        elif column.upper() == 'PRIORITY':
            newRow.append(np.max(targets['priority'])+1)
        elif column.lower() in defaultValues:
            newRow.append(defaultValues[column.lower()])
        else:
            newRow.append(-999)

    targets.add_row(newRow)

    log.important('autocomplete: added target with mangaid=' +
                  target['MANGAID'] + ' (ifudesignsize=' +
                  str(int(bundleSize)) + ')')

    return targets