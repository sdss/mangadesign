#!/usr/bin/env python
# encoding: utf-8
"""
PlateTargets.py

Created by José Sánchez-Gallego on 26 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    26 Mar 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

from Gohan.exceptions import GohanPlateTargetsError, GohanPlateTargetsWarning
from Gohan.utils import utils
from Gohan.utils import yanny
from Gohan import log, readPath, config

from astropy import table
import numpy as np

import six
import os


def _toLowerCase(structure):
    """Renames all the columns in an `astropy.table.Table` to lowercase."""

    for colname in structure.colnames:
        if colname != colname.lower():
            structure.rename_column(colname, colname.lower())

    return structure


# Fields from mangaScience to exclude in ancillary plateTargets.
excludeFields = ['priority', 'sourcetype', 'manga_tileids']


# Field names conversion
conversions = {
    'all': {
        'nsa_redshift': 'z',
        'nsa_mstar': 'mass',
        'nsa_petro_th50': 'petroth50'
    }
}


neverobserve = list(map(int, open(readPath(config['plateTargets']['neverobserve']),
                                  'r').read().splitlines()[1:]))

try:
    mangaTargetsExtNSA = table.Table.read(utils.get_manga_targets_path())

    mangaTargetsExtNSA['MANGAID'] = list(map(lambda xx: xx.strip(), mangaTargetsExtNSA['MANGAID']))
    mangaTargetsExtNSA = _toLowerCase(mangaTargetsExtNSA)
except:
    raise GohanPlateTargetsError('failed loading targeting catalogue')


class PlateTargets(object):

    def __init__(self, input, **kwargs):
        """A class to handle plateTargets files.

        Parameters
        ----------
        input : integer or string
            Either the catalogid or the path to the plateTargets file to use.

        """

        self.template = False
        self._nAppended = 0

        if isinstance(input, six.string_types):
            self.path = input
            if not os.path.exists(self.path):
                raise GohanPlateTargetsError(
                    'path {0} cannot be found'.format(self.path))

        else:
            catalogid = input
            self.path = utils.getPlateTargetsPath(catalogid)
            if not os.path.exists(self.path):
                self.path = utils.getPlateTargetsTemplate(catalogid)
                if self.path is None:
                    raise GohanPlateTargetsError(
                        'neither the plateTargets nor the template '
                        'for catalogid={0} can be found'.format(catalogid))
                log.warning('using template for catalogid={0}'.format(catalogid),
                            GohanPlateTargetsWarning)
                self.template = True

        data = yanny.yanny(self.path, np=True)

        self.comments = self._getComments(data)

        if self.template:
            self.catalogid = catalogid
        else:
            typeid = int(data['typeid'])
            self.catalogid = typeid

        self.ancillary = True if self.catalogid >= 30 else False

        self.structure = table.Table(data['PLTTRGT'])

        if self.template:
            self.structure.remove_row(0)

        log.debug('using targeting catalogue {0}'.format(utils.get_manga_targets_path()))

    def _getComments(self, data):
        """Gets the comments in the Yanny file."""

        header = []
        for row in data._contents.split('\n'):
            if row.strip() == '':
                header.append('')
            elif row.strip()[0] == '#':
                header.append(row)
            else:
                break

        # Strips empty lines
        header = '\n'.join(list(map(str, header))).strip() + '\n\n'

        return header

    def write(self, path=None, useCatID=True):
        """Writes the current instance to a Yanny file.

        Parameters
        ----------
        path : string or None, optional
            The path where the plateTargets file will be written.
            If `path=None`, the usual path will be used.

        Returns
        -------
        result : tuple
            A tuple containing, in the first element, the path where the
            file was writen, and in the second element the number of new
            targets added to the plateTargets file.
        useCatID : bool
            If False, not typeid will be written in the Yanny file. This is
            mainly to be used when StarPlateTargets super calls this method.

        """

        if len(self.structure) == 0:
            raise GohanPlateTargetsError('there are no targets added to the '
                                         'plateTargets file.')

        if path is None:
            if self.template:
                path = utils.getPlateTargetsPath(self.catalogid)
            else:
                path = self.path

        if os.path.exists(path):
            os.remove(path)

        # Writes the yanny file.
        if useCatID:
            yanny.write_ndarray_to_yanny(path, self.structure.as_array(),
                                         structname='PLTTRGT',
                                         hdr={'typeid': self.catalogid})
        else:
            yanny.write_ndarray_to_yanny(path, self.structure.as_array(),
                                         structname='PLTTRGT')

        # Now adds the comments on top of the file.
        lines = [line for line in open(path, 'r').read().splitlines()
                 if len(line) == 0 or line[0] != '#']
        lines = self.comments + '\n'.join(lines) + '\n'
        unit = open(path, 'w')
        unit.write(lines)
        unit.close()

        if useCatID:
            log.debug('plateTargets-{0} saved to {1}'
                      .format(self.catalogid, path))
            log.debug('{0} targets appended to plateTargets-{1}'
                      .format(self._nAppended, self.catalogid))

        # If this was a template, it stops being it after writing the file
        # to disk.
        if self.template:
            self.path = path
            self.template = False

        nAppended = self._nAppended
        self._nAppended = 0

        return (self.path, nAppended)

    def addTargets(self, mangaids=None, plateid=None, mangaScience=None,
                   overwrite=False):
        """Adds a new target to the current instance.

        The method requires a `mangaid` and either a `plateid` or a
        `mangaScience`. In `plateid` is defined but not `mangaScience`,
        the mangaScience will be determined from the plateHolesSorted file.
        If only `mangaScience` is defined, the plate-related information will
        be filled with nulls. A warning will be issued.

        This method gathers information from mangaScience, catalogue data and
        plateHolesSorted. It then calls either the routine for the main sample,
        or for ancillary targets. The returned values are then processed for
        targetfixes and added to the main structure.

        Parameters
        ----------
        mangaids : string, list of strings or None
            The mangaid or list of mangaids to add. If None, all the mangaids
            in the plate.
        plateid : integer or None, optional
            The plateid of the plate on which the new target(s) have been
            drilled. If `mangaids` is a list, it is assumed that all the
            targets correspond to the same `plateid`.
        mangaScience : string or None, optional
            The path of the mangaScience for the new target. As with `plateid`,
            it is assumed that all the `mangaids` are included in that
            mangaScience file.

        Returns
        -------
        result : `astropy.table.Table`
            The new rows added to the `PlateTargets` instance, as well
            as returned.

        """

        addedIndices = []

        if plateid is None and mangaScience is None:
            raise GohanPlateTargetsError(
                'either plateid or mangaScience needs to be defined')

        # If plateid is not defined, issues a warning.
        if plateid is None:
            log.warning('no plateid information provided. The plate '
                        'information will be filled out with nulls.',
                        GohanPlateTargetsWarning)

            plateHolesSortedData = plateHolesSortedPairs = None

        # If plateid is defined we get the plateHolesSorted data.
        if plateid is not None:

            if mangaScience is None:
                # If mangaScience is not defined, we use
                mangaScience = utils.getMangaSciencePath(plateid, format='plateid')

            plateHolesSortedData, plateHolesSortedPairs = utils.getPlateHolesSortedData(plateid)

        # Checks mangaScience path
        assert os.path.exists(mangaScience), ('path {0} does not exist'.format(mangaScience))

        # Gets data from mangaScience structure and pairs
        mangaScienceData, mangaSciencePairs = self._readMangaScience(mangaScience)
        mangaScienceData = _toLowerCase(mangaScienceData)

        if mangaids is None:
            mangaids = mangaScienceData['mangaid']

        # Makes sure that all the mangaids have the catalogid of the current
        # instance of PlateTargets.
        catalogids = np.array([int(mangaid.split('-')[0]) for mangaid in mangaids])

        if self.catalogid not in ['MaSTAR', 'standard']:
            assert np.all(catalogids == self.catalogid), \
                'one or more of the mangaids has catalogid != {0}'.format(
                    self.catalogid)

        # Defines mangaids either from input or from the mangaScience data.
        if mangaids is not None:
            mangaids = [mangaid.strip() for mangaid in np.atleast_1d(mangaids)]
        else:
            mangaids = [mangaid.strip()
                        for mangaid in mangaScienceData['mangaid']]

        # Gets the common fields for all plateTargets.
        commonData = self._getCommonData(
            mangaids, mangaScienceData, mangaSciencePairs,
            plateHolesSortedData, plateHolesSortedPairs)

        if self.catalogid in ['MaSTAR', 'standard']:
            return commonData, mangaids, plateid

        # Creates a list of specific fields for this plateTargets file
        specificColumns = [col for col in self.structure.colnames
                           if col not in utils.getRequiredPlateTargetsColumns()
                           ]

        # Now we get the specific fields for this catalogid.
        if self.catalogid in [1, 12]:
            # We call the specific method for the main sample.
            specificData = self._getMainSampleData(mangaids, specificColumns)
        else:
            # Completes specific columns using information in mangaScience
            # specificData = self._getAncillaryData(mangaids, mangaScienceData,
            #                                       commonData)
            specificData = {}

        # Finally we add the data target by target.
        for mangaid in mangaids:

            # We combine both dictionaries
            targetData = commonData[mangaid]
            if mangaid in specificData:
                targetData.update(specificData[mangaid])

            # Checks if the targets already exists in plateTargets.
            existing = False

            if plateid is not None:

                # Checks if the tuple (mangaid, plateid) already exists in
                # plateTargets.
                plateTargetRow = self.structure[
                    (self.structure['mangaid'] == mangaid) &
                    (self.structure['plateid'] == plateid)]

                if len(plateTargetRow) > 0:
                    # If it exists, checks if overwrite is True
                    if overwrite:
                        existing = True
                        log.warning('replacing target mangaid={0} in plateid={1}'
                                    .format(mangaid, plateid), GohanPlateTargetsWarning)
                    else:
                        # If overwrite is False, skips this target.
                        log.debug('skipping mangaid={0} because it is already '
                                  'in plateTargets.'.format(mangaid))
                        continue

            # Cleans up values
            targetData = self._cleanupTargetData(targetData)

            # Applies target fixes
            targetData = self._applyTargetFix(targetData)

            # Check that redshift is well populated
            if targetData['z'] == -999. and 'nsa_z' in targetData:
                targetData['z'] = targetData['nsa_z']

            if targetData['iauname'] == -999:
                targetData['iauname'] = '-999'

            # Adds the new targets
            if not existing:
                self.structure.add_row(targetData)
                addedIndices.append(len(self.structure) - 1)
            else:
                # If the target already exists, replaces it values
                idx = np.where((self.structure['mangaid'] == mangaid) &
                               (self.structure['plateid'] == plateid))

                for field in targetData:
                    self.structure[field][idx] = targetData[field]
                addedIndices.append(idx[0][0])

            self._nAppended += 1

            if self.catalogid == 'MaSTAR':
                fileName = 'starPlateTargets'
            elif self.catalogid == 'standard':
                fileName = 'standardPlateTargets'
            else:
                fileName = 'plateTargets-{0}'.format(self.catalogid)

            log.debug('mangaid={0} added to {1}'.format(mangaid, fileName))

        return self.structure[addedIndices]

    def _getCommonData(self, mangaids, mangaScienceData, mangaSciencePairs,
                       plateHolesSortedData, plateHolesSortedPairs):
        """Returns a dictionary with the required plateTargets columns.

        Creates a dictionary with keys each one of the `mangaids`. For each
        element, the value is another dictionary containing the mandatory
        columns, common to all plateTargets files.

        """

        result = {}

        requiredColumns = utils.getRequiredPlateTargetsColumns()

        designid = int(plateHolesSortedPairs['designid'])

        if self.catalogid == 12:
            targets12 = readPath(config['plateTargets']['targets12'])
            nsaV1bCat = _toLowerCase(table.Table.read(targets12))

        for mangaid in mangaids:
            result[mangaid] = {}

            # We get the appropriate row in mangaScience
            assert mangaid in mangaScienceData['mangaid'], \
                'mangaid={0} not found in mangaScience file'.format(mangaid)
            mangaScienceRow = mangaScienceData[
                mangaScienceData['mangaid'] == mangaid]

            # Idem with plateHolesSorted, if defined
            if plateHolesSortedData is not None:
                assert mangaid in plateHolesSortedData['mangaid'], \
                    'mangaid={0} not found in plateHolesSorted file'.format(
                        mangaid)
                plateHolesSortedRow = plateHolesSortedData[
                    plateHolesSortedData['mangaid'] == mangaid]
            else:
                plateHolesSortedRow = None

            # Gets the row for the mangaid from te targeting catalogue
            if self.catalogid not in ['MaSTAR', 'standard']:
                targetRow = mangaTargetsExtNSA[
                    mangaTargetsExtNSA['mangaid'] == mangaid.strip()]
                if len(targetRow) == 0:
                    log.warning('mangaid={0} not found in targeting catalogue'.format(mangaid),
                                GohanPlateTargetsWarning)
                    targetRow = None
            else:
                targetRow = None

            for column in requiredColumns:

                # Handles the neverobserve column. Note that if plateid is None
                # (plateHolesSorted not defined), neverobserve defaults to 0.
                if column == 'neverobserve':
                    result[mangaid]['neverobserve'] = \
                        1 if designid in neverobserve else 0
                    continue

                # We get the coordinate information
                if column in ['ifu_ra', 'ifu_dec', 'object_ra', 'object_dec',
                              'catalog_ra', 'catalog_dec']:
                    result[mangaid][column] = \
                        self._getCoordinate(column, targetRow, mangaScienceRow,
                                            plateHolesSortedRow)
                    continue

                if column == 'epoch':
                    result[mangaid][column] = plateHolesSortedRow[column]
                    continue
                elif column == 'z':
                    if targetRow is not None:
                        result[mangaid][column] = targetRow['nsa_z'][0]
                    else:
                        result[mangaid][column] = -999.
                    continue

                # Now uses plateHolesSorted, the targeting catalogue, and
                # mangaScience, in that order, to complete the information.
                if (plateHolesSortedPairs is not None and
                        column in plateHolesSortedPairs):
                    result[mangaid][column] = plateHolesSortedPairs[column]
                elif (plateHolesSortedRow is not None and
                        column in plateHolesSortedRow.colnames):
                    result[mangaid][column] = plateHolesSortedRow[column]
                elif column in mangaSciencePairs:
                    result[mangaid][column] = mangaSciencePairs[column]
                elif targetRow is not None and column in targetRow.colnames:
                    result[mangaid][column] = targetRow[column][0]
                elif column in mangaScienceRow.colnames:
                    result[mangaid][column] = mangaScienceRow[column][0]
                elif (self.catalogid == 12 and column in nsaV1bCat.colnames and
                        mangaid in nsaV1bCat['mangaid']):
                    result[mangaid][column] = nsaV1bCat[
                        nsaV1bCat['mangaid'] == mangaid][column]

                else:
                    if 'manga_target' in column:
                        # If manga_targetX is not defined, we set it to 0.
                        result[mangaid][column] = 0
                    else:
                        # If other value is not found, we set it to -999 but
                        # raise a warning. Skips columns that we know will
                        # systematically fail for stellar library targets
                        if (self.catalogid not in ['MaSTAR', 'standard'] or
                            column not in ['iauname', 'ifutargetsize',
                                           'ifudesignwrongsize']):
                            log.warning('mangaid={0}: no value found for mandatory field {1}. '
                                        'Setting it to -999'.format(mangaid, column),
                                        GohanPlateTargetsWarning)
                            result[mangaid][column] = -999

        return result

    def _getCoordinate(self, column, targetRow, mangaScienceRow,
                       plateHolesSortedRow):
        """Returns the appropriate coordinate for a column and target."""

        coord_type, coord_name = column.split('_')

        if coord_type == 'ifu':
            return plateHolesSortedRow['target_' + coord_name][0]

        # Handles MaStar a bit differently
        if self.catalogid == 'MaSTAR':
            if coord_type == 'object':
                if column in mangaScienceRow.colnames:
                    return mangaScienceRow[column]
                else:
                    mastar_obj = ('obj' + coord_name)
                    if mastar_obj in mangaScienceRow.colnames:
                        log.debug('using {} to populate {}'.format(mastar_obj, column))
                        return mangaScienceRow[mastar_obj]

        if targetRow is not None:
            return targetRow[column]
        else:
            return mangaScienceRow[coord_name]

    def _getMainSampleData(self, mangaids, fields):
        """Returns the specific columns for the MaNGA main sample.

        Returns a dictionary similar to _getCommonData but with the columns
        specific to plateTargets-1.par. It also handles plateTargets-12.par.

        """

        # Reads the NSA v1b to v1 match
        nsaV1bToV1 = table.Table.read(readPath(config['plateTargets']['NSA_v1b_v1']),
                                      format='ascii.fixed_width',
                                      delimiter='|')

        result = {}
        for mangaid in mangaids:

            result[mangaid] = {}
            mangaidDict = result[mangaid]

            targetID = int(mangaid.split('-')[1])

            nsaCatPath = utils.getCataloguePath(self.catalogid)

            # We want to get the right row from NSA v1_0_1. If catalogid=1, we
            # just use utils.getCatalogueRow with the mangaid. If catalogid=12
            # we create a mock mangaid with the format 1-{indexID_101} taking
            # NSAID_v101 from the list of conversions in nsaV1bToV1.

            if self.catalogid == 1:
                nsa101_mangaid = mangaid
                indexID_101 = targetID
            elif self.catalogid == 12:
                nsaV1bToV1_row = nsaV1bToV1[nsaV1bToV1['MaNGAID'] == mangaid]
                indexID_101 = nsaV1bToV1_row['catid'][0]
                nsa101_mangaid = '1-{0}'.format(indexID_101)

            nsaRow = utils.getCatalogueRow(nsa101_mangaid)
            targetRow = mangaTargetsExtNSA[
                mangaTargetsExtNSA['mangaid'] == nsa101_mangaid.strip()]

            # If this is from catalogid=12, let's make a quick check and make
            # sure we have selected the right column.
            if self.catalogid == 12:
                assert nsaRow['nsaid'] == nsaV1bToV1_row['NSAID_v1_0_0'][0], \
                    ('row selected from NSA v1_0_1 does not match the '
                     'value expected for NSAID ({0} != {1})'.format(
                         nsaRow['nsaid'], nsaV1bToV1_row['NSAID_v1_0_1'][0]))

            if nsaRow is None:
                raise GohanPlateTargetsError(
                    'mangaid={0} cannot be found in catalogue {1}'
                    .format(mangaid, nsaCatPath))

            for field in fields:

                if field == 'nsa_version':
                    mangaidDict[field] = self._getNSAVersion()

                elif field == 'nsa_nsaid_v1b':
                    # For commissioning targets, we also include the nsaid
                    # in the nsa_v1b catalogue.
                    if self.catalogid == 12:
                        mangaidDict[field] = utils.getCatalogueRow(mangaid)['nsaid']
                    else:
                        mangaidDict[field] = -999

                elif field in config['plateTargets']['mangaTargetFields']:
                    if len(targetRow) == 0:
                        mangaidDict[field] = -999.
                    else:
                        mangaidDict[field] = targetRow[field][0]

                elif 'nsa_' in field:
                    nsaField = field[4:].lower()
                    mangaidDict[field] = nsaRow[nsaField]

                else:
                    raise GohanPlateTargetsError(
                        'unexpected field {0} when compiling data for '
                        'plateTargets-{1}. Could not find field {2} in '
                        'NSA catalogue'
                        .format(field, self.catalogid, field))

        return result

    def _getAncillaryData(self, mangaids, mangaScienceData, commonData):
        """Returns a dictionary with all the columns in mangaScience.

        This method returns a dictionary with the data in the mangaScience
        structure that are not in the common fields. It is mainly used for
        ancillary plateTargets that do not have a set list of specific columns.

        """

        result = {}
        for mangaid in mangaids:

            result[mangaid] = {}
            mangaidDict = result[mangaid]

            row = mangaScienceData[mangaScienceData['mangaid'] == mangaid]

            for field in mangaScienceData.colnames:

                if field in commonData[mangaid]:
                    # If the column is in the common fields, skips it because
                    # we have already compiled this information
                    continue

                if field not in self.structure.colnames:
                    # If the fields is new, adds the column to the structure
                    # only if the structure is empty. Once the structure
                    # columns are set, they won't change.
                    if len(self.structure) == 0:
                        newCol = table.Column(
                            None, name=field,
                            dtype=mangaScienceData[field].descr[1],
                            shape=mangaScienceData[field].descr[2])
                        self.structure.add_column(newCol)
                    else:
                        continue

                mangaidDict[field] = row[field]

        return result

    def _readMangaScience(self, mangaSciencePath):
        """Reads a mangaScience file and returns the structure and keywords."""

        assert os.path.exists(mangaSciencePath)

        mangaScienceData = yanny.yanny(mangaSciencePath, np=True)

        if 'MANGAINPUT' in mangaScienceData.keys():
            structName = 'MANGAINPUT'
        elif 'STRUCT1' in mangaScienceData.keys():
            structName = 'STRUCT1'
        else:
            raise GohanPlateTargetsError(
                'cannot identify structure name for mangaScience {0}'
                .format(mangaSciencePath))

        if 'manga_tileid' in mangaScienceData:
            pass
        elif 'tileid' in mangaScienceData:
            mangaScienceData['manga_tileid'] = mangaScienceData.pop('tileid')
        else:
            raise GohanPlateTargetsError(
                'no maga_tileid or tileid field found in mangaScience {0}'
                .format(mangaSciencePath))

        structure = _toLowerCase(
            table.Table(mangaScienceData.pop(structName)))
        mangaScienceData.pop('symbols')

        # Strips mangaids in structure. Important for mangaid comparisons
        for row in structure:
            row['mangaid'] = row['mangaid'].strip()

        # return the mangaScience structure and all the keywords.
        return (structure, mangaScienceData)

    def _getNSAVersion(self):
        """Returns the version of the NSA catalogue used."""

        if self.catalogid == 1:
            return 'v1_0_1'
        elif self.catalogid == 12:
            return 'v1_0_1'
        else:
            return '-999'

    def _cleanupTargetData(self, targetData):
        """Does some cleaning on the targetData dictionary."""

        for key in targetData:

            # Strips strings
            if isinstance(targetData[key], six.string_types):
                targetData[key] = targetData[key].strip()

            # Checks if a value can be converted to float
            try:
                float(targetData[key])
                targetData[key] = float(targetData[key])
            except:
                pass

            try:
                # Checks if value(s) are NaN
                isNaN = np.isnan(targetData[key])
                if np.any(isNaN):
                    if np.isscalar(isNaN):
                        targetData[key] = -999.
                    else:
                        targetData[key][isNaN] = -999.

                # Does the same for infs
                isInf = np.isinf(targetData[key])
                if np.any(isInf):
                    if np.isscalar(isInf):
                        targetData[key] = -999.
                    else:
                        targetData[key][isInf] = -999.
            except TypeError:
                pass

        return targetData

    def _applyTargetFix(self, targetData):
        """Reads the targetfix file for a target and updates the
        target data."""

        plateid = int(targetData['plateid'])

        if plateid is None or plateid == -999:
            return targetData

        targetFixPath = utils.getTargetFix(plateid)
        if targetFixPath is None:
            return targetData

        # Reads targetfix and makes sure all mangaids are stripped
        targetFix = yanny.yanny(targetFixPath, np=True)['OPTARFIX']
        for row in targetFix:
            row['mangaid'] = row['mangaid'].strip()

        # Determines the rows in targetfix referring to the mangaid of
        # targetData
        mangaid = targetData['mangaid'][0]
        targetFix_mangaid = targetFix[targetFix['mangaid'] == mangaid]

        # If targetfixes are found for our mangaid, proceeds to update
        # targetData.
        if len(targetFix_mangaid) > 0:
            log.important('Applying target fix to mangaid={0}'.format(mangaid))
            for row in targetFix_mangaid:
                keyword = row['keyword']
                if keyword in targetData:
                    oldValue = targetData[keyword]

                    newValue = row['value'].strip()

                    # Case that the value is an array in IDL format
                    if newValue[0] == '{':
                        newValue = np.fromstring(
                            newValue[1:-1], dtype=float, sep=' ')

                    targetData[keyword] = newValue
                    log.debug('mangaid={0}: {1}={2} -> {3}'.format(
                              mangaid, keyword, oldValue, newValue))

        return targetData
