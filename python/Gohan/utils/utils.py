#!/usr/bin/env python
# encoding: utf-8
"""
utils.py

Created by José Sánchez-Gallego on 11 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    11 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
from Gohan import readPath, config, log
import glob
from Gohan.utils import yanny
from Gohan.exceptions import GohanError, GohanUserWarning
from numbers import Number
from collections import OrderedDict
import numpy as np
from astropy import table
import warnings


# __all__ = [
#     'getMaskBitFromLabel',
#     'getPlateTargetsPath',
#     'getPlateListDir',
#     'getPlates',
#     'getPlateDir',
#     'getPlateHolesSortedPath',
#     'getMaNGAIDs',
#     'sortByCatID',
#     'getParsingFunction',
#     'getSampleCatalogue',
#     'getPlateInputData',
#     'getPlateHolesSortedData',
#     'getPointing',
#     'getMangaScience',
#     'getPlateDefinition',
#     'getTargetFix',
#     'getDesignID',
#     'getCatalogueRow',
#     'getCataloguePath',
#     'getPlateTargetsTemplate'
# ]


# Dictionary to cache catalogues after being read
cachedCatalogues = {}


try:
    sdssMaskBitsFile = readPath(config['sdssMaskBits'])
    sdssMaskBits = table.Table(
        yanny.yanny(sdssMaskBitsFile, np=True)['MASKBITS'])
except:
    sdssMaskBits = None


def getMaskBitFromLabel(flag, label):
    """Returns the maskbit and description associated to a label for a
    certain flag. E.g., getMaskBitFromLabel('MANGA_TARGET1', 'FILLER')."""

    flag = flag.upper()
    label = label.upper()

    if sdssMaskBits is None:
        raise GohanError('sdssMaskBits could not be found. Use gohan.yaml '
                         'to define sdssMaskBits.')

    assert flag in sdssMaskBits['flag'], \
        'could not find {0} in sdssMaskBits flags'.format(flag)

    flaggedBits = sdssMaskBits[sdssMaskBits['flag'] == flag]

    assert label in flaggedBits['label'], \
        'could not find {0} in {1} labels'.format(label, flag)

    row = flaggedBits[flaggedBits['label'] == label][0]

    return row['bit'], row['description']


def getPlateTargetsPath(catID):

    plateTargetsPath = os.path.join(
        readPath(config['mangacore']), 'platedesign/platetargets/',
        'plateTargets-{0}.par'.format(catID))

    return plateTargetsPath


def getPlateListDir():

    plateListDir = readPath(config['platelist'])
    if 'trunk' in glob.glob(os.path.join(plateListDir, '*')):
        plateListDir = os.path.join(plateListDir, 'trunk')

    return plateListDir


plateListDir = getPlateListDir()
platePlansFile = os.path.join(plateListDir, 'platePlans.par')

if os.path.exists(platePlansFile):
    platePlans = yanny.yanny(platePlansFile, np=True)['PLATEPLANS']
else:
    platePlans = None


def getFromPlatePlans(plateRun, column='plateid'):
    """Returns a column (defaults to plateid) for all the plates in a certain
    plate run."""

    if platePlans is None:
        raise GohanError('platePlans file not found')

    return sorted(
        platePlans[np.where(platePlans['platerun'] == plateRun)][column])


def getPlateDir(plateid):

    if not isinstance(plateid, (Number, np.int32)):
        raise GohanError('plate must be a number')

    plateid = str(int(plateid)).zfill(4)

    platePath = os.path.join(getPlateListDir(),
                             'plates/00{0}XX/00{1}'.format(plateid[0:2],
                                                           plateid))

    if not os.path.exists(platePath):
        raise GohanError(
            'no plate path found for plateid {0}'.format(plateid))

    return platePath


def getPlateHolesSortedPath(plateid):

    if not isinstance(plateid, (Number, np.int32)):
        raise GohanError('plate must be a number')

    plateDir = getPlateDir(plateid)

    plateid = str(int(plateid)).zfill(6)
    plateHolesSortedPath = os.path.join(
        plateDir, 'plateHolesSorted-{0}.par'.format(plateid))

    if not os.path.exists(plateHolesSortedPath):
        raise GohanError(
            'no plateHolesSorted found for plate {0}'.format(plateid))

    return plateHolesSortedPath


def getMaNGAIDs(plateid):

    ynPlateHolesSorted = yanny.yanny(
        getPlateHolesSortedPath(plateid), np=True)['STRUCT1']

    return [ynPlateHolesSorted['mangaid'][ii]
            for ii in range(len(ynPlateHolesSorted))
            if ynPlateHolesSorted['holetype'][ii] == 'MANGA' and
            ynPlateHolesSorted['targettype'][ii] == 'science']


def sortByCatID(mangaids):

    catIDdict = OrderedDict()

    for mangaid in mangaids:
        catID = int(mangaid.split('-')[0])
        if catID in catIDdict.keys():
            catIDdict[catID].append(mangaid)
        else:
            catIDdict[catID] = [mangaid]

    return catIDdict


def getParsingFunction(catID):

    if catID == 1:
        from .catalogueParsers import parseCatalogID_1
        return parseCatalogID_1
    elif catID == 12:
        from .catalogueParsers import parseCatalogID_12
        return parseCatalogID_12
    elif catID == 50:
        from .catalogueParsers import parseCatalogID_50
        return parseCatalogID_50
    else:
        return None


def getSampleCatalogue(catID):

    catalogIDs = os.path.join(
        readPath(config['mangacore']),
        'platedesign/catalog_ids.dat')

    if not os.path.exists(catalogIDs):
        raise GohanError('no catalog_ids.dat found')

    catNames = table.Table.read(catalogIDs, format='ascii.no_header',
                                names=['catID', 'file'])

    catFile = catNames[catNames['catID'] == catID]['file']

    if len(catFile) == 0:
        raise GohanError('no catalogID={0} found in catalog_ids.dat'
                         .format(catID))

    catPath = os.path.join(
        readPath(config['cataloguePath']),
        '{0}-{1}'.format(catID, catFile[0]))

    if os.path.exists(catPath):
        return catPath
    else:
        if os.path.exists(catPath + '.gz'):
            return catPath + '.gz'
        else:
            raise GohanError('catalogID={0} could not be found in path {1}'
                             .format(catID, catPath))


def getPlateInputData(mangaid, plateHolesSorted):

    inputs = []
    for key in plateHolesSorted.keys():
        if 'plateinput' in key:
            inputs.append(plateHolesSorted[key])

    for input in inputs:

        plateInputPath = os.path.join(
            readPath(config['platelist']), 'inputs', input)

        plateInput = yanny.yanny(plateInputPath, np=True)

        if 'MANGAINPUT' in plateInput.keys():
            structName = 'MANGAINPUT'
        elif 'STRUCT1' in plateInput.keys():
            structName = 'STRUCT1'
        else:
            # warnings.warn('cannot identify structure name for plateInput '
            #               '{0}'.format(plateInputPath), GohanUserWarning)
            continue

        if 'manga_tileid' in plateInput:
            manga_tileid = plateInput['manga_tileid']
        elif 'tileid' in plateInput:
            manga_tileid = plateInput['tileid']
        else:
            raise GohanError('no maga_tileid or tileid field found in '
                             'plateInput {0}'.format(plateInputPath))

        tbStructure = table.Table(
            plateInput[structName],
            names=[nn.lower() for nn in plateInput[structName].dtype.names])

        for row in tbStructure:

            if 'mangaid' not in row.colnames:
                raise GohanError('no mangaid field found in plateInput '
                                 '{0}'.format(plateInputPath))

            if row['mangaid'] == mangaid:
                return {'MANGAINPUT': row, 'racen': plateInput['racen'],
                        'deccen': plateInput['deccen'],
                        'designid': plateInput['designid'],
                        'manga_tileid': manga_tileid}

    raise GohanError('it has not been possible to find a plateInput for '
                     'mangaid={0}'.format(mangaid))


def getPlateHolesSortedData(plateid):
    """Returns a formatted version of the plateHolesSorted data for a
    plateid."""

    plateHolesSortedPath = getPlateHolesSortedPath(plateid)

    plateHolesSorted = yanny.yanny(plateHolesSortedPath, np=True)
    plateHolesSortedStruct = table.Table(plateHolesSorted.pop('STRUCT1'))

    for colname in plateHolesSortedStruct.colnames:
        if colname != colname.lower():
            plateHolesSortedStruct.rename_column(colname, colname.lower())

    plateHolesSorted.pop('symbols')
    for key in plateHolesSorted:
        if key != key.lower():
            plateHolesSorted[key.lower()] = plateHolesSorted.pop(key)

    # Strips mangaids
    for row in plateHolesSortedStruct:
        row['mangaid'] = row['mangaid'].strip()

    return (plateHolesSortedStruct, plateHolesSorted)


def getPointing(plateInputData):

    fields = ['ra', 'dec', 'target_ra', 'target_dec', 'ifu_ra', 'ifu_dec']
    ra = plateInputData['ra']
    dec = plateInputData['dec']

    for field in fields:
        if field not in plateInputData.colnames:
            log.debug('some pointing info missing for mangaid={0}. Using '
                      'only ra and dec.'.format(plateInputData['mangaid']))
            return {'target_ra': ra, 'target_dec': dec,
                    'ifu_ra': ra, 'ifu_dec': dec}

    return {'target_ra': plateInputData['target_ra'],
            'target_dec': plateInputData['target_dec'],
            'ifu_ra': plateInputData['ifu_ra'],
            'ifu_dec': plateInputData['ifu_dec']}


def getMangaScience(input, format='designid'):
    """Returns the path of the mangaScience data for a design or plate."""

    if format == 'plateid':
        designID = getDesignID(input)
    else:
        designID = input

    plateDefinition = getPlateDefinition(designID)

    plateDefYanny = yanny.yanny(plateDefinition)

    for key in plateDefYanny.keys():
        if 'plateInput' in key:
            if 'mangaScience' in plateDefYanny[key]:
                return os.path.join(readPath(config['platelist']), 'inputs',
                                    plateDefYanny[key])

    raise ValueError('no plateInput file found')


def getPlateDefinition(designID):
    """Returns the path of the plateDefinition for a designid."""

    path = os.path.join(
        readPath(config['platelist']),
        'definitions/00{0}XX/plateDefinition-{1:06d}.par'
        .format('{0:04d}'.format(designID)[0:2], designID))

    if os.path.exists(path):
        return path
    else:
        raise ValueError('plateDefinition does not exist.')


def getTargetFix(plateID):

    path = os.path.join(
        readPath(config['mangacore']), 'platedesign/targetfix/',
        '{0:04d}XX'.format(int(plateID/100)),
        'targetfix-{0}.par'.format(plateID))

    if not os.path.exists(path):
        return None
    else:
        return path


def getDesignID(plateid):
    """Returns the designid for a plateid."""

    if platePlans is None:
        raise GohanError('platePlans file not found')

    try:
        return platePlans[platePlans['plateid'] == plateid]['designid'][0]
    except:
        raise GohanError('plateid={0} not found'.format(plateid))


def getCatalogueRow(mangaid, catalogue=None):
    """Recovers the row in the parent catalogue matching the mangaid.
    An `astropy.table.Table` instance of the catalogue can be passed."""

    # List of catalogids for catalogues in which the targetid in mangaid is
    # the index (zero-indexed) in the parent catalogue.
    indexedCatalogues = [1, 8, 12]

    catalogID, targetID = mangaid.strip().split('-')

    if catalogue is None:
        if catalogID in cachedCatalogues:
            # If catalogue is cached
            catalogue = cachedCatalogues[catalogID]

        else:
            # Reads catalogue and caches is

            cataloguePath = getCataloguePath(catalogID)

            if cataloguePath is None:
                return None

            catalogue = table.Table.read(cataloguePath)

            # Changes all columns to lowercase
            for col in catalogue.colnames:
                if col != col.lower():
                    catalogue.rename_column(col, col.lower())

            cachedCatalogues[catalogID] = catalogue

    if int(catalogID) in indexedCatalogues:
        return catalogue[int(targetID)]
    else:
        if 'MANGAID' not in catalogue.colnames:
            warnings.warn('MANGAID column not be found', GohanUserWarning)
            return None

        row = catalogue[catalogue['MANGAID'] == mangaid]
        if len(row) == 1:
            return row
        else:
            warnings.warn('no row found in catalogue for mangaid={0}'
                          .format(mangaid), GohanUserWarning)
            return None


def getCataloguePath(catalogid):
    """Returns the path of a catalogue for a certain catalogid."""

    catalogid = int(catalogid)

    catalogidsPath = os.path.join(readPath(config['mangacore']),
                                  'platedesign/catalog_ids.dat')
    if not os.path.exists(catalogidsPath):
        warnings.warn('catalog_ids.dat could not be found', GohanUserWarning)
        return None

    catalogids = open(catalogidsPath, 'r').read().splitlines()

    cataloguePath = None
    for catalogidRow in catalogids:
        if int(catalogidRow.strip().split()[0]) == catalogid:
            cataloguePath = catalogidRow
            break

    if cataloguePath is None:
        warnings.warn('no entry in catalog_ids.dat for catalogid={0}'
                      .format(catalogid), GohanUserWarning)
        return None

    cataloguePath = os.path.join(
        readPath(config['catalogues']['catalogueDir']),
        '-'.join(cataloguePath.strip().split()))

    if not os.path.exists(cataloguePath):
        warnings.warn('no catalogue found in {0}'.format(cataloguePath),
                      GohanUserWarning)
        return None

    return cataloguePath


def getPlateTargetsTemplate(catalogid):
    """Returns the path to the plateTargets template for a given catalogid."""

    if catalogid in [1, 12]:
        return readPath('+plateTargets/plateTargets-1.template')
    elif catalogid >= 30:
        return readPath('+plateTargets/plateTargets-Ancillary.template')
    else:
        warnings.warn('no template found for catalogid={0}'.format(catalogid),
                      GohanUserWarning)
        return None
