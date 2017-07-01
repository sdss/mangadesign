#!/usr/bin/env python
# encoding: utf-8
"""
assignIFUDesigns.py

Created by José Sánchez-Gallego on 10 Nov 2014.
Licensed under a 3-clause BSD license.

Revision history:
    10 Nov 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from astropy import coordinates as coo
from astropy import table

from Gohan import config, log, readPath
from Gohan.exceptions import GohanError
from Gohan.utils import yanny


try:
    cartmap = yanny.yanny(
        readPath(config['plateInputs']['cartmap']), np=True)['CMAP']
except:
    cartmap = None


def ifuDesign2Size(ifudesign):
    row = cartmap[cartmap['ifuDesign'] == ifudesign][0]
    if len(row) == 0:
        raise GohanError('Invalid ifudesign {0}'.format(ifudesign))
    return row['ifusize']


def assignIFUDesigns(targets, centre, targettype='science',
                     failOnIFUDesignSize=False, reassignAll=False,
                     plot=False, plotFilename=None, **kwargs):
    """Checks a list of targets and assigns, if needed, ifudesignsizes and
    ifudesigns in an efficient way."""

    assert cartmap is not None
    assert len(centre) == 2

    if targettype == 'science':
        bundleSizes = config['IFUs'].copy()
    elif targettype == 'standard':
        bundleSizes = config['miniBundles'].copy()
    elif targettype == 'sky':
        bundleSizes = config['skies'].copy()

    for size in targets['ifudesignsize']:
        if size > 0:
            bundleSizes[size] -= 1

    ifuDesignSizeAssigned = False
    for target in targets:
        if target['ifudesignsize'] < 0:
            if failOnIFUDesignSize:
                raise GohanError('at least one target found with '
                                 'ifudesignsize < 0')
            for size in bundleSizes:
                if bundleSizes[size] > 0:
                    target['ifudesignsize'] = size
                    target['ifudesign'] = -999
                    ifuDesignSizeAssigned = True
                    bundleSizes[size] -= 1
                    log.debug(
                        'mangaid={0} has been automatically '.format(
                            target['mangaid']) +
                        'assigned an ifudesignsize={0}'.format(
                            target['ifudesignsize']))
                    break

    if np.sum(bundleSizes.values()) > 0:
        raise GohanError('some bundles have not been assigned.')
    elif np.sum(bundleSizes.values()) < 0:
        raise GohanError('more bundles assigned than available.')

    if ifuDesignSizeAssigned:
        log.info('some ifudesignsizes have been assigned automatically.')

    # Checks if one or more ifuDesigns are missing.
    missingDesign = False
    for ifuDesign in targets['ifudesign']:
        if ifuDesign < 0.:
            missingDesign = True
            break

    # If at least one ifuDesign is undefined, reassign ferrules. If
    # reassignFerrules is True, reassign all the ferrules even if they
    # are defined in the input catalogue.
    if missingDesign:
        log.debug('one or all ifudesigns are missing. Reassigning IFUs.')
        targets = assignFerrulesAnchorBlock(targets, centre,
                                            reassignAll=reassignAll)
    elif reassignAll:
        log.info('reassigning all IFUs.')
        targets = assignFerrulesAnchorBlock(targets, centre, reassignAll=True)
    else:
        log.debug('all ifudesigns are correctly assigned')

    if plot:
        plotIFUs(targets, centre, filename=plotFilename, **kwargs)

    return targets


def assignFerrulesAnchorBlock(targets, centre, reassignAll=False):

    struct1 = targets.copy()
    raCen, decCen = centre

    # Gets the targets for which IFUs need to be assigned.
    if not reassignAll:
        struct1ToKeep = struct1[struct1['ifudesign'] > 0]
        struct1ToAssign = struct1[struct1['ifudesign'] <= 0]
    else:
        struct1ToKeep = None
        struct1ToAssign = struct1.copy()
        struct1ToAssign['ifudesign'] = -999

    assignOrder, blockPositionsCoo = _getAssignOrder(struct1ToAssign,
                                                     raCen, decCen)

    # Gets already assigned ifudesigns
    usedIFUDesigns = []
    if struct1ToKeep is not None:
        for row in struct1ToKeep:
            usedIFUDesigns.append(int(row['ifudesign']))

    for idx in assignOrder:

        row = struct1ToAssign[idx]

        ra = row['ra']
        dec = row['dec']
        targetCoo = coo.SkyCoord(ra=ra, dec=dec, unit='deg')

        ifuSize = row['ifudesignsize']
        ifuDesign = _assignFerrule(
            targetCoo, ifuSize, blockPositionsCoo, usedIFUDesigns)

        usedIFUDesigns.append(ifuDesign)
        struct1ToAssign[idx]['ifudesign'] = ifuDesign
        log.debug('mangaid={0} assigned ifudesign={1}'.format(
            struct1ToAssign[idx]['mangaid'],
            struct1ToAssign[idx]['ifudesign']))

    if struct1ToKeep is None:
        return struct1ToAssign
    else:
        if len(struct1ToKeep) == 0:
            return struct1ToAssign
        elif len(struct1ToAssign) == 0:
            return struct1ToKeep
        else:
            return table.vstack((struct1ToKeep, struct1ToAssign))


def _getAssignOrder(struct, raCen, decCen):
    """Calculates the distances of the targets to the anchor blocks.

    Returns the indices of the input structure sorted by distance to
    the anchor blocks (ordered from farther to closer). The metric used
    is the square root of the sum of the squares of the distance of the
    target to each block. The (RA, Dec) position of the blocks is also
    returned in `time` format.

    """

    if len(struct) == 0:
        return struct

    distances = np.zeros(len(struct), dtype=float)

    blockPositions = getRADecAnchorBlocks(raCen, decCen)
    blockPositionsCoo = coo.SkyCoord(
        ra=list(zip(*blockPositions)[0]),
        dec=list(zip(*blockPositions)[1]),
        unit='deg')

    for nn, row in enumerate(struct):
        ra = row['ra']
        dec = row['dec']

        targetCoo = coo.SkyCoord(ra=ra, dec=dec, unit='deg')
        separations = targetCoo.separation(blockPositionsCoo)

        distances[nn] = np.sum(separations.degree**2)

    return np.argsort(distances)[::-1], blockPositionsCoo


def _assignFerrule(targetCoo, ifuSize,
                   blockPositionsCoo, usedIFUDesigns):
    """
    Returns the optimum ifuDesign for a target based on the target
    position, the desired IFU size and the number of available ifuDesigns.

    """

    ifus = config['IFUs'].copy()
    ifus.update(config['miniBundles'])

    for usedIFUDesign in usedIFUDesigns:
        usedIFUSize = ifuDesign2Size(usedIFUDesign)
        ifus[usedIFUSize] -= 1

    if ifus[ifuSize] <= 0:
        raise ValueError('trying to assign too many '
                         'IFUs of size {0}.'.format(ifuSize))

    anchorBlockOrder = np.argsort(
        targetCoo.separation(blockPositionsCoo).degree)

    sortedBlocks = [config['anchorBlocks']['anchors'][anchorIdx]
                    for anchorIdx in anchorBlockOrder]
    ifuDesign = 0
    for anchorName in sortedBlocks:
        ferrules = config['anchorBlocks'][anchorName]
        for ferrule in ferrules:
            cartMapRow = cartmap[cartmap['frlPlug'] == ferrule][0]
            ferruleSize = cartMapRow['ifusize']
            ferruleDesign = cartMapRow['ifuDesign']
            if ferruleDesign not in usedIFUDesigns and \
                    ferruleSize == ifuSize:
                ifuDesign = ferruleDesign
                break
            else:
                continue
        if ifuDesign != 0:
            break

    return ifuDesign


def getRADecAnchorBlocks(ra, dec):
    """Calculates RA and Dec of the anchor blocks.

    Returns a list of tuples, each tuple containing RA and Dec.
    The order of blocks in the returned list is [NE, SE, NW, SE]

    """

    raAngToAnchor = 1.5 / 2. / np.cos(dec * np.pi / 180.)
    decAngToAnchor = 1.5 / 2

    return [
        (ra + raAngToAnchor, dec + decAngToAnchor),
        (ra + raAngToAnchor, dec - decAngToAnchor),
        (ra - raAngToAnchor, dec + decAngToAnchor),
        (ra - raAngToAnchor, dec - decAngToAnchor)]


def plotIFUs(structures, centre=None, filename='ifuPlot.pdf', **kwargs):
    """Plots the position of the IFUs on the plate.

    This method plots the position of the IFUs (and, thus,
    of the targets) on the plate with information about
    MaNGAID and IFUDESIGN.

    Parameters
    ----------
    struct1 : `Table` or list of Tables or PlateInput instances
        The targets to plot.
    centre : tuple
        The coordinates of the centre of the plot.
    filename : str, optional
        The plot filename.

    """

    from Gohan import PlateInput  # To avoid circular import

    skies = []

    if isinstance(structures, table.Table):
        if centre is None:
            raise GohanError('centre is required')
        tables = [structures]
    else:
        if isinstance(structures[0], PlateInput):
            if centre is None:
                centre = (structures[0].raCen, structures[0].decCen)
            tables = []
            for structure in structures:
                if structure.targettype == 'sky':
                    skies.append(structure.mangaInput)
                else:
                    tables.append(structure.mangaInput)

    assert len(tables) > 0

    struct1 = tables[0]['mangaid', 'ra', 'dec', 'ifudesign', 'ifudesignsize']
    for tt in tables[1:]:
        tmpData = tt['mangaid', 'ra', 'dec', 'ifudesign', 'ifudesignsize']
        struct1 = table.vstack([struct1, tmpData], join_type='exact')

    if len(skies) > 0:
        skyStruct1 = table.vstack([sky['ra', 'dec'] for sky in skies],
                                  join_type='exact')
    else:
        skyStruct1 = None

    plt.cla()
    plt.clf()

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)

    raCen, decCen = centre

    plate = Ellipse((raCen, decCen),
                    height=3.,
                    width=3 / np.cos(decCen * np.pi / 180.),
                    linewidth=2,
                    edgecolor='k', facecolor='None')
    ax.add_patch(plate)

    for row in struct1:

        ra = row['ra']
        dec = row['dec']
        ifuSize = row['ifudesignsize']
        ifuDesign = row['ifudesign']
        mangaID = row['mangaid']

        cartMapInfo = cartmap[cartmap['ifuDesign'] == ifuDesign][0]
        ferrule = cartMapInfo['frlPlug']

        anchorName = ''
        for key in config['anchorBlocks']['anchors']:
            if ferrule in config['anchorBlocks'][key]:
                anchorName = key
                break

        closestAnchor = _isClosestAnchor(
            ra, dec, raCen, decCen, anchorName)

        if closestAnchor:
            color = 'b'
        else:
            color = 'r'

        ax.scatter(ra, dec, s=ifuSize / 2, marker='o',
                   edgecolor=color, facecolor='None', zorder=10)

        ax.text(ra, dec - 0.06,
                '{0} ({1})'.format(ifuDesign, anchorName),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=8, color=color, zorder=10)

        ax.text(ra, dec + 0.05, r'{0}'.format(mangaID),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=6, color=color, zorder=10)

    if skyStruct1 is not None:
        ax.scatter(skyStruct1['ra'], skyStruct1['dec'], marker='x', s=6,
                   edgecolor='0.8', zorder=0)

    ax.set_xlim(raCen + 1.6 / np.cos(decCen * np.pi / 180.),
                raCen - 1.6 / np.cos(decCen * np.pi / 180.))
    ax.set_ylim(decCen - 1.6, decCen + 1.6)

    ax.set_xlabel(r'$\alpha_{2000}$')
    ax.set_ylabel(r'$\delta_{2000}$')

    plt.savefig(filename)

    plt.close('all')

    log.info('IFU plot for saved to {0}'.format(filename))


def _isClosestAnchor(ra, dec, raCen, decCen, anchorName):
    """Returns True if an ifuDesign is from the closest anchor."""

    blockPositions = getRADecAnchorBlocks(raCen, decCen)
    blockPositionsCoo = coo.SkyCoord(
        ra=list(zip(*blockPositions)[0]),
        dec=list(zip(*blockPositions)[1]),
        unit='deg')

    targetCoo = coo.SkyCoord(
        ra=ra, dec=dec,
        unit='deg')

    closestAnchorIdx = np.argmin(
        targetCoo.separation(blockPositionsCoo).degree)
    closestAnchor = config['anchorBlocks']['anchors'][closestAnchorIdx]

    if closestAnchor == anchorName:
        return True
    else:
        return False
