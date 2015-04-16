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
from Gohan.utils.sortTargets import sortTargets
from Gohan.utils import getMaskBitFromLabel
import warnings
import numpy as np
import os


try:
    NSAPath = readPath(config['catalogues']['NSA'])

    if not os.path.exists(NSAPath):
        warnings.warn('NSA catalogue {0} not found'.format(NSAPath),
                      GohanUserWarning)
        NSACat = None
    else:
        NSACat = table.Table.read(NSAPath)
except:
    NSACat = None

qaWarningIssued = False

defaultValues = {
    'ifudesign': -999,
    'ifudesignsize': -999,
    'manga_target1': 2**getMaskBitFromLabel('MANGA_TARGET1', 'FILLER')[0],
    'manga_target2': 0,
    'manga_target3': 0,
    'psfmag': [0.0, 0.0, 0.0, 0.0, 0.0],
    'sourcetype': None
}


def autocomplete(targets, centre, **kwargs):
    """Autocompletes a list of targets using the MaNGa sample or the
    NSA catalogue."""

    bundles = config['IFUs'].copy()
    nBundles = np.sum(bundles.values())

    if len(targets) >= nBundles:
        log.info('the are at least as many targets as bundles. No need '
                 'to run autocomplete.')
        return targets

    nBundlesToAssign = nBundles - len(targets)
    log.info('autocompleting {0} bundles.'.format(nBundlesToAssign))

    MaNGATargets = getMaNGATargets(targets, centre)

    NSATargets = getNSATargets(targets, centre)
    if NSATargets is not None:
        NSATargetCoords = np.zeros((len(NSATargets), 2), np.float64)
        NSATargetCoords[:, 0] = NSATargets['RA']
        NSATargetCoords[:, 1] = NSATargets['DEC']
        sortedCoords, order = sortTargets(NSATargetCoords, centre)
        NSATargets = NSATargets[order]

    unassignedBundleSizes = getUnassignedBundleSizes(targets, bundles)

    for bundleSize in unassignedBundleSizes:

        target = None

        # Tries first the MaNGA sample catalogue.
        if MaNGATargets is not None:
            target = _getOptimalTarget(targets, MaNGATargets,
                                       bundleSize, centre, **kwargs)

        if target is None and NSATargets is not None:
            # If no targets can be found in the MaNGA sample, uses NSA.
            target = _getOptimalTarget(targets, NSATargets,
                                       bundleSize, centre, **kwargs)

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


def _getOptimalTarget(targets, candidateTargets, bundleSize, centre,
                      useReff=True, **kwargs):

    reffField = config['plateMags']['reffField'].upper()
    if reffField in candidateTargets.colnames and useReff:
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

    for target in candidateTargets:
        # First we try to use the IFUDESIGNSIZE column. Not that if all values
        # in IFUDESIGNSIZE are -999, all targets will be skipped.
        if (target['IFUDESIGNSIZE'] >= bundleSize and
                _checkNSATarget(targets, target, centre)):
            target['IFUDESIGNSIZE'] = bundleSize
            return target

    # If the previous process fails, we try using reff and a minimum galaxy
    # size depending on the bundle size.

    reffField = config['plateMags']['reffField'].upper()
    candidateTargets.sort(reffField)
    candidateTargets.reverse()

    if bundleSize == 19:
        minSize = 2 * 2.5
    elif bundleSize == 37:
        minSize = 3 * 2.5
    elif bundleSize == 61:
        minSize = 4 * 2.5
    elif bundleSize == 91:
        minSize = 5 * 2.5
    elif bundleSize == 127:
        minSize = 6 * 2.5
    else:
        minSize = 0.

    for target in candidateTargets:
        if (target[reffField] >= minSize and
                _checkNSATarget(targets, target, centre)):
            target['IFUDESIGNSIZE'] = bundleSize
            return target

    # If it comes down to this, it means that there are no target with valid
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

    if NSACat is None:
        return None

    # Adds the CATID column to the NSA catalogue for comparison with MaNGA
    # targets
    if 'CATIND' not in NSACat.colnames:
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


def getMaNGATargets(targets, centre):
    """Gets non-selected targets from the general MaNGA sample catalogue."""

    MaNGAPath = readPath(config['catalogues']['science'])

    if not os.path.exists(MaNGAPath):
        warnings.warn('MaNGA catalogue {0} not found'.format(MaNGAPath),
                      GohanUserWarning)
        return None

    MaNGACat = table.Table.read(MaNGAPath)

    raCen, decCen = centre

    coords = coo.SkyCoord(MaNGACat['RA'], MaNGACat['DEC'], unit='deg')
    separation = coords.separation(
        coo.SkyCoord(ra=raCen, dec=decCen, unit='deg')).deg

    MaNGATargets = MaNGACat[np.where(separation <=
                                     config['decollision']['FOV'])]

    MaNGATargets.add_column(table.Column([-999] * len(MaNGATargets),
                                         name='IFUDESIGN', dtype=int))

    mangaIDs = [target['MANGAID'] for target in targets]

    removeIdx = []
    for nn in range(len(MaNGATargets)):
        if MaNGATargets[nn]['MANGAID'] in mangaIDs:
            removeIdx.append(nn)
    MaNGATargets.remove_rows(removeIdx)

    MaNGATargets.sort('Z')

    if len(MaNGATargets) > 0:
        return MaNGATargets
    else:
        return None


def _checkNSATarget(targets, target, centre):

    # Gets mangaids with bad photometry
    badPhotometry = getBadPhotometry()

    if target['MANGAID'] in badPhotometry:
        log.debug('target {0} rejected because it has bad photometry'
                  .format(target['MANGAID']))
        return False

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


def getBadPhotometry():
    """If available, reads the QA_adjusted file linked to the science catalogue
    and returns a list of mangaids with bad photometry."""

    global qaWarningIssued

    scienceCat = readPath(config['catalogues']['science'])

    qaAdjusted = os.path.join(
        os.path.dirname(scienceCat), 'nsa_v1_0_0_QA_adjusted.dat')

    if not os.path.exists(qaAdjusted):
        if not qaWarningIssued:
            # Issues this warning only once
            warnings.warn('nsa_v1_0_0_QA_adjusted.dat not found',
                          GohanUserWarning)
            qaWarningIssued = True
        return []

    qaTable = table.Table.read(qaAdjusted, format='ascii.commented_header')

    # Adds the catalogid and returns the targets with bad photometry
    badPhotometry = qaTable[qaTable['bad_phot'] == 1]

    return ['1-{0}'.format(str(catind)) for catind in badPhotometry['catind']]


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
    row = targets[targets['MANGAID'] == target['MANGAID']][0]

    log.important('autocomplete: added target with mangaid={0}'
                  ' (ifudesignsize={1}, manga_target1={2})'
                  .format(target['MANGAID'], str(int(bundleSize)),
                          row['MANGA_TARGET1']))

    return targets
