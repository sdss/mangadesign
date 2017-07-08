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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import glob
import os
import warnings

from collections import OrderedDict
from numbers import Number

import numpy as np

from astropy import table
from astropy import time
from astropy.io import fits

from Gohan import readPath, config, log
from Gohan.exceptions import GohanError, GohanUserWarning
from Gohan.utils import yanny


# Dictionary to cache catalogues and plateTargets after being read
cachedCatalogues = {}
cachedPlateTargets = {}


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
        readPath(config['mangacore']), 'platedesign/platetargets/')

    if catID == 'MaStar':
        plateTargetsPath = os.path.join(plateTargetsPath,
                                        'starPlateTargets.par')
    elif catID == 'standards':
        plateTargetsPath = os.path.join(plateTargetsPath,
                                        'standardPlateTargets.par')
    else:
        plateTargetsPath = os.path.join(plateTargetsPath,
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
    platePlans = table.Table(yanny.yanny(platePlansFile, np=True)['PLATEPLANS'])
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

    plateidDir = str(int(plateid)).zfill(6)
    plateidPre = plateidDir[0:4] + 'XX'

    platePath = os.path.join(getPlateListDir(),
                             'plates/{0}/{1}'.format(plateidPre, plateidDir))

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
    plateHolesSorted = dict((kk.lower(), vv) for kk, vv in plateHolesSorted.items())

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


def getMangaSciencePath(input, format='designid'):
    """Returns the path of the mangaScience data for a design or plate."""

    return getPlateInputPath(input, mode='science', format=format)


def getPlateInputPath(input, mode='science', format='designid'):
    """Returns the path of the plateInput data for a design or plate.

    Parameters
    ----------
    input : int
        Either the plateid or designid of the design to be used.
    format : str
        Either `'designid'` or `'plateid'`. Defines the type of `input`.
    mode : str
        One of `'science', 'standard', 'sky'`.

    """

    if format == 'plateid':
        designID = getDesignID(input)
    else:
        designID = input

    plateDefinition = getPlateDefinition(designID)

    plateDefYanny = yanny.yanny(plateDefinition)

    mangaInput = 'manga' + mode.capitalize()

    for key in plateDefYanny.keys():
        if 'plateInput' in key:
            if ((mangaInput in plateDefYanny[key]) or
                    ('plateInput' in plateDefYanny[key] and
                     designID in [7878, 7879, 7888])):
                return os.path.join(readPath(config['platelist']), 'inputs',
                                    plateDefYanny[key])

    raise ValueError('no plateInput file found')


def getPlateDefinition(designID):
    """Returns the path of the plateDefinition for a designid."""

    path = os.path.join(
        readPath(config['platelist']),
        'definitions/{0}/plateDefinition-{1:06d}.par'
        .format('{0:06d}'.format(designID)[0:4] + 'XX', designID))

    if os.path.exists(path):
        return path
    else:
        raise ValueError('plateDefinition does not exist.')


def getTargetFix(plateID, alwaysReturnPath=False):

    path = os.path.join(
        readPath(config['mangacore']), 'platedesign/targetfix/',
        '{0:04d}XX'.format(int(plateID / 100)),
        'targetfix-{0}.par'.format(plateID))

    if alwaysReturnPath:
        return path

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
    indexedCatalogues = [1, 12, 8]

    catalogID, targetID = mangaid.strip().split('-')

    if catalogue is None:
        if catalogID in cachedCatalogues:
            # If catalogue is cached
            catalogue = cachedCatalogues[catalogID]

        else:
            # Reads catalogue and caches it
            cataloguePath = getCataloguePath(catalogID)

            if cataloguePath is None:
                return None

            fData = fits.open(cataloguePath)
            catalogue = table.Table(fData[1].data)

            # Changes all columns to lowercase
            for col in catalogue.colnames:
                if col != col.lower():
                    catalogue.rename_column(col, col.lower())

            cachedCatalogues[catalogID] = catalogue

    if int(catalogID) in indexedCatalogues:
        return catalogue[int(targetID)]
    else:
        if 'mangaid' not in catalogue.colnames:
            warnings.warn('mangaid column not be found', GohanUserWarning)
            return None

        row = catalogue[catalogue['mangaid'] == mangaid]
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

    if catalogid == 1:
        return readPath('+plateTargets/plateTargets-1.template')
    elif catalogid == 12:
        return readPath('+plateTargets/plateTargets-12.template')
    elif catalogid == 50:
        return readPath('+plateTargets/plateTargets-50.template')
    elif catalogid >= 30 or catalogid == 18:
        return readPath('+plateTargets/plateTargets-Ancillary.template')
    else:
        warnings.warn('no template found for catalogid={0}'.format(catalogid),
                      GohanUserWarning)
        return None


def getAPOGEEPlateTemperature(date):
    """Returns the drilling temperature for a plate.

    Calculated the optimal temperature for which a plate must be drilled if it
    is planned to be observed at a certain date, to minimise the possibility
    that the scale of the plate falls out of the limits.

    Uses the APOGEE formula from get_drilltemp.pro.

    Parameters
    ----------
    date : float or `astropy.time.Time` instance
        Either an integer with the year fraction (e.g., 0.37 for May 15th) or
        an `astropy.time.Time` object with the date at which the plate will be
        observed

    Returns
    -------
    temperature : float
        The optimal temperature for which the plate must be drilled.

    """

    # Max/min reference temperatures for each month in Fahrenheit. Last element
    # is the same as the first one for wrapping.
    maxtempsf = np.array([38.7, 40.3, 46.7, 55.5, 66.3, 74.0, 71.8, 68.9, 65.7,
                          57.2, 46.7, 40.6, 38.7])
    mintempsf = np.array([21.5, 20.6, 24.6, 31.0, 39.4, 46.9, 50.3, 48.6, 44.9,
                          36.2, 27.4, 22.7, 21.5])

    # Converts to Celsius
    maxtemps = (maxtempsf - 32.) * 5. / 9.
    mintemps = (mintempsf - 32.) * 5. / 9.

    # Gets the month fraction.
    if isinstance(date, time.Time):
        month = 12. * (date.byear - int(date.byear))
    else:
        month = 12. * date

    # Just makes sure the month fraction is between 0 and 12
    month = month % 12.

    # Interpolates min/max temp
    interpMaxTemp = np.interp(month, np.arange(13), maxtemps)
    interpMinTemp = np.interp(month, np.arange(13), mintemps)

    temperature = (interpMaxTemp - interpMinTemp) / 3. + interpMinTemp

    return temperature


def getStellarLibraryTargets(designid=None, reject_neverobserved=True, reject_saturated=True):
    """Returns a plateTargets-like table with stellar library targets.

    Parameters:
        designid (list or None):
            If a list, only the targets for those designids will be returned.
            If None, targets for all stellar library plates are returned.
        reject_neverobserved (bool):
            If True, targets in plates listed in ``stellibtargets/etc/neverobserved.txt``
            are not included.
        reject_saturated (bool):
            targets that meet the condition that
            ``(psfmag[:, 3] < 12.7 ) and (psfmag[:, 3] > 0) and
              (abs(ifu_dec - object_dec) < 0.1 / 3600.)``
            are not returned.

    Returns:
        targets (astropy.table.Table object):
            A Astropy table, with the same format as the ``PLTTRGT`` structure in
            starPlateTargets, but filtered to include only the desired targets.

    """

    starPlateTargetsPath = getPlateTargetsPath('MaStar')
    starPlateTargets = table.Table(yanny.yanny(starPlateTargetsPath, np=True)['PLTTRGT'])

    if reject_neverobserved:
        neverobserved = getMaStarNeverobserved()
        for plateid in neverobserved:
            starPlateTargets.remove_rows(np.where(starPlateTargets['plateid'] == plateid))

    if reject_saturated:
        idx_remove = np.where((starPlateTargets['psfmag'][:, 3] < 12.7) &
                              (starPlateTargets['psfmag'][:, 3] > 0) &
                              (np.abs((starPlateTargets['ifu_dec'] -
                                       starPlateTargets['object_dec'])) < 0.1 / 3600))
        starPlateTargets.remove_rows(idx_remove)

    if designid is None:
        return starPlateTargets
    else:
        return starPlateTargets[starPlateTargets['designid'] == designid]


def getMaStarNeverobserved():
    """Returns a list of plates that will never be observed by MaStar."""

    if ('STELLIBTARGETS_DIR' not in os.environ or
            not os.path.exists(os.environ['STELLIBTARGETS_DIR'])):
        raise ValueError('$STELLIBTARGETS_DIR is undefined or missing.')

    neverobserved_path = os.path.join(os.environ['STELLIBTARGETS_DIR'],
                                      'etc/neverobserved.txt')

    if not os.path.exists(neverobserved_path):
        raise ValueError('cannot find file {0}'.format(neverobserved_path))

    neverobserved = table.Table.read(neverobserved_path,
                                     format='ascii.commented_header')

    return neverobserved['Plate'].tolist()


def getStellarLibraryDesigns():
    """Returns, from apodb, all the designids for stellar library plates."""

    apogeeMangaRuns = getStellarLibraryRuns()

    apogeeMangaDesigns = []
    for run in apogeeMangaRuns:
        apogeeMangaDesigns += (platePlans[platePlans['platerun'] == run]
                               ['designid'].tolist())

    return apogeeMangaDesigns


def getStellarLibraryRuns():
    """Returns all APOGEE2-MaNGA plate runs."""

    plateRuns = np.unique(platePlans['platerun'])

    apogeeMangaRuns = [plateRun for plateRun in plateRuns
                       if 'apogee2-manga' in plateRun.strip().lower() and
                       plateRun.strip().lower() != '2014.04.f.apogee2-manga']

    return apogeeMangaRuns


def getRequiredPlateTargetsColumns():
    """Returns a list with the mandatory plateTargets columns."""

    path = os.path.join(os.path.dirname(getPlateTargetsPath(1)),
                        'requiredColumns.dat')

    assert os.path.exists(path), 'requiredColumns.dat not found in mangacore'

    return open(path, 'r').read().splitlines()


def getLastLocationID():
    """Returns the value of the last locationid used for plate design.

    This is a convenience function that simply returns the value of
    Gohan/etc/lastLocationID.dat

    """

    file = readPath('+etc/lastLocationID.dat')

    assert os.path.exists(file), 'lastLocationID.dat not found'

    lines = open(file, 'r').read().splitlines()

    for line in lines:
        lineStripped = line.strip()
        if lineStripped != '' and lineStripped[0] != '#':
            return int(line)

    raise ValueError('the format of lastLocationID.dat does not seem valid.')


def getAllMaNGAPlateRuns():
    """Returns a list of all plate runs."""

    plateRuns = np.unique(platePlans['platerun'])

    mangaRuns = [plateRun for plateRun in plateRuns if 'manga' in plateRun.lower()]

    mangaLeadRuns = []
    for mangaRun in mangaRuns:
        if mangaRun in ['2014.02.f.manga', '2014.04.e.manga-apogee2']:
            continue
        year = int(mangaRun.split('.')[0])
        if year < 2014:
            continue
        surveys = mangaRun.split('.')[-1].split('-')
        if len(surveys) == 1:
            mangaLeadRuns.append(mangaRun)
        else:
            if surveys[0] == 'manga':
                mangaLeadRuns.append(mangaRun)

    return mangaLeadRuns


def getAllMaNGAPlates():
    """Returns a list of all MaNGA plates. Rejects some old plates."""

    runs = getAllMaNGAPlateRuns()

    plates = []
    for run in runs:
        runPlates = platePlans[platePlans['platerun'] == run]['plateid']
        plates += runPlates.tolist()

    return plates


def getPlateID(designID):
    """Returns the plateid associated with a designid."""

    plateID = platePlans[platePlans['designid'] == designID]['plateid']

    if len(plateID) == 0:
        return None

    return plateID[0]


def getPlateTargetsRow(mangaid, plateid=None):
    """Returns the row in plateTarget for a certain target.

    Note that if plateid is None, more than one row may be returned."""

    mangaid = mangaid.strip()

    catalogid = int(mangaid.split('-')[0])
    if catalogid in cachedPlateTargets:
        plateTargets = cachedPlateTargets[catalogid]
    else:
        plateTargets = yanny.yanny(getPlateTargetsPath(catalogid), np=True)
        plateTargets = table.Table(plateTargets['PLTTRGT'])
        cachedPlateTargets[catalogid] = plateTargets

    if plateid is not None:
        idx = np.where((plateTargets['mangaid'] == mangaid) &
                       (plateTargets['plateid'] == plateid))
    else:
        idx = np.where(plateTargets['mangaid'] == mangaid)

    row = plateTargets[idx]

    if len(row) == 0:
        return None
    else:
        return row


def getNSArow(nsaid):
    """Returns the NSA row for a certain nsaid."""

    if 1 in cachedCatalogues:
        cat = cachedCatalogues[1]
    else:
        ff = fits.open(getCataloguePath(1))
        cat = table.Table(ff[1].data)
        cachedCatalogues[1] = cat

    row = cat[cat['NSAID'] == nsaid]

    if len(row) == 0:
        return None
    else:
        return row


def getAllDrilledTargets():
    """Returns a stacked table with the plateTargets for all targets."""

    allPlateTargets = glob.glob(
        os.path.join(
            readPath(config['mangacore']), 'platedesign/platetargets/',
            'plateTargets-*.par'))

    requiredColumns = getRequiredPlateTargetsColumns()
    tables = [table.Table(yanny.yanny(pT)['PLTTRGT'])[requiredColumns]
              for pT in allPlateTargets]

    return table.vstack(tables)


def get_manga_targets_path():
    """Returns the full MaNGA_targets path from the configuration file."""

    path = os.path.expandvars(os.path.join(config['targets']['path'], config['targets']['version'],
                                           config['targets']['mangaTargets']))

    return path


def calculate_special_targets(plate_input_path, input_path):
    """Determines how many special MaStar targets survived decollision.

    Parameters:
        plate_input_path (str):
            The path to the generated plateInput file/
        input_path (str):
            The path to the original FITS table containing all the potential
            targets. The table should contain a ``CATALOG`` column. Targets for
            which ``CATALOG='Special'`` will be compared to determined whether
            they are also include in ``plate_input_path``.

    Returns:
        output:
            Either ``None`` if ``input_path`` does not contain a ``CATALOG``
            column, or a tuple of two elements, the first being the number of
            special targets in ``plate_input_path``, and the second the total
            number of original targets in ``input_path``.

    """

    assert os.path.exists(plate_input_path), 'cannot find file {0}'.format(plate_input_path)
    assert os.path.exists(input_path), 'cannot find file {0}'.format(input_path)

    plate_input = table.Table(yanny.yanny(plate_input_path, np=True)['MANGAINPUT'])
    input_data = table.Table.read(input_path)

    if 'CATALOG' not in input_data.colnames:
        return None

    # Strips some columns for good measure
    input_data['MANGAID'] = list(map(lambda xx: xx.strip(), input_data['MANGAID']))
    input_data['CATALOG'] = list(map(lambda xx: xx.strip(), input_data['CATALOG']))

    special_targets = input_data[input_data['CATALOG'] == 'Special']
    n_special = len(special_targets)

    if n_special == 0:
        return (0, 0)

    n_survive = 0

    for m_id in special_targets['MANGAID']:
        if m_id in plate_input['mangaid']:
            n_survive += 1

    return (n_survive, n_special)


def print_special_summary(plate_data_path):
    """Prints a table summarising the MaStar special targets that survived decollision.

    This function assumes that the MaStar input files are in a directory called
    ``input`` in the same path as ``plate_data_path`` and that the target
    files conform to the MaStar standards.

    Parameters:
        plate_data_path (str):
            The path to the APOGEE plate data file.

    """

    fields = table.Table.read(plate_data_path, format='ascii.commented_header')

    special_table = table.Table(None,
                                names=['design_id', 'field_name', 'special_targets', 'allocated'],
                                dtype=[int, 'U20', int, int])

    for field in fields:

        designID = int(field['DesignID'])
        fieldName = field['FieldName']

        plate_input_path = 'mangaScience_{0}_{1}.par'.format(fieldName, designID)
        input_path = os.path.join(os.path.dirname(plate_data_path),
                                  'inputs/{0}_target.fits'.format(fieldName))

        special = calculate_special_targets(plate_input_path, input_path)

        if special is None:
            special_table.add_row((designID, fieldName, 0, 0))
        else:
            special_table.add_row((designID, fieldName, special[1], special[0]))

    print()
    special_table.pprint()

    return special_table


def platerun_is_mastar(platerun):
    """Returns True if the platerun is a MaStar one."""

    manga_runs = getAllMaNGAPlateRuns()

    if platerun in manga_runs:
        return False

    surveys = platerun.split('.')[3].split('-')

    if len(surveys) == 1:
        raise ValueError('platerun {0} is not a MaNGA run '
                         'but is not co-observed.'.format(platerun))
    else:
        if 'manga' in surveys and 'apogee' in platerun:
            return True
        else:
            raise ValueError('cannot determine type for platerun {0}'.format(platerun))
