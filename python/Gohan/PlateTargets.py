#!/usr/bin/env python
# encoding: utf-8
"""
function.py

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
import warnings
import os


# Fields to exclude
excludeFields = ['priority', 'sourcetype', 'manga_tileids']


# Field names conversion
conversions = {
    'all': {
        'nsa_redshift': 'z',
        'nsa_mstar': 'mass',
        'nsa_petro_th50': 'petroth50'
    }
}


neverobserve = map(
    int,
    open(readPath(config['plateTargets']['neverobserve']),
         'r').read().splitlines()[1:])


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

        if isinstance(input, basestring):
            self.path = basestring
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
                warnings.warn(
                    'using template for catalogid={0}'.format(catalogid),
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
        header = '\n'.join(map(str, header)).strip() + '\n\n'

        return header

    def write(self, path=None):
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

        """

        if path is None:
            if self.template:
                path = utils.getPlateTargetsPath(self.catalogid)
            else:
                path = self.path

        if os.path.exists(path):
            os.remove(path)

        # Writes the yanny file.
        yanny.write_ndarray_to_yanny(path, self.structure.as_array(),
                                     structname='PLTTRGT',
                                     hdr={'typeid': self.catalogid})

        # Now adds the comments on top of the file.
        lines = [line for line in open(path, 'r').read().splitlines()
                 if len(line) == 0 or line[0] != '#']
        lines = self.comments + '\n'.join(lines) + '\n'
        unit = open(path, 'w')
        unit.write(lines)
        unit.close()

        log.debug('plateTargets-{0} saved to {1}'.format(self.catalogid, path))
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

    def addTargets(self, mangaids, plateid=None, mangaInput=None,
                   overwrite=False):
        """Adds a new target to the current instance.

        The method requires a `mangaid` and either a `plateid` or a
        `mangaInput`. In `plateid` is defined but not `mangaInput`,
        the mangaInput will be determined from the plateHolesSorted file.
        If only `mangaInput` is defined, the plate-related information will be
        filled. with nulls A warning will be issued.

        It uses a combination of mangaScience, catalogue data and
        plateHolesSorted (in that order) to obtain all possible information
        about the target and fill out all the plateTargets fields.

        Parameters
        ----------
        mangaids : string or list of strings
            The mangaid or list of mangaids to add.
        plateid : integer or None, optional
            The plateid of the plate on which the new target(s) have been
            drilled. If `mangaids` is a list, it is assumed that all the
            targets correspond to the same `plateid`.
        mangaInput : string or None, optional
            The path of the mangaInput for the new target. As with `plateid`,
            it is assumed that all the `mangaids` are included in that
            mangaInput file.

        Returns
        -------
        result : `astropy.table.Table`
            The new rows added to the `PlateTargets` instance.

        """

        mangaids = np.atleast_1d(mangaids)

        if mangaInput is not None:
            assert os.path.exists(mangaInput), ('path {0} does not exist'
                                                .format(mangaInput))

        if plateid is None:
            warnings.warn(
                'no plateid information provided. The plate information will '
                'be filled out with nulls.', GohanPlateTargetsWarning)

        if plateid is None and mangaInput is None:
            raise GohanPlateTargetsError(
                'either plateid or mangaInput needs to be defined')

        if plateid is not None and mangaInput is None:
            mangaInput = utils.getMangaScience(plateid)

        # Gets data from mangaInput structure and keywords
        mangaInputStructure, mangaInputKyw = self._readMangaInput(mangaInput)
        mangaInputStructure = self._toLowerCase(mangaInputStructure)

        designid = int(mangaInputKyw['designid'])

        # Reads the plateHolesSorted file
        if plateid is not None:
            plateHolesSortedStructure, plateHolesSortedKyw = \
                utils.getPlateHolesSortedData(plateid)
        else:
            plateHolesSortedStructure = plateHolesSortedKyw = None

        # Defines the fields to fill out
        # In principle, the fields are the columns of the current instance.
        fields = self.structure.colnames

        # If this is a template of an ancillary catalogue, we add all the
        # fields from the mangaInput structure.
        if self.template and self.ancillary:
            for ii in range(len(mangaInputStructure.columns)):
                column = mangaInputStructure.columns[ii].copy()
                colName = column.name
                if colName not in fields and colName not in excludeFields:
                    # newColumn = table.Column(None, name=colName,
                    #                          dtype=column.dtype)
                    self.structure.add_column(column[0:0])
                    fields.append(column.name)

        for mangaid in mangaids:

            # Defines if the targets already exists in plateTargets.
            existing = False

            if plateid is not None:
                if (mangaid in self.structure['mangaid'] and
                        plateid in self.structure['plateid']):
                    if overwrite:
                        existing = True
                        warnings.warn(
                            'replacing target mangaid={0} in plateid={1}'
                            .format(mangaid, plateid),
                            GohanPlateTargetsWarning)
                    else:
                        continue

            # Selects the row for mangaid from the mangaScience structure
            if mangaid not in mangaInputStructure['mangaid']:
                raise GohanPlateTargetsError(
                    'mangaid={0} not found in mangaInput'.format(mangaid))

            mangaInputRow = mangaInputStructure[
                mangaInputStructure['mangaid'] == mangaid][0]

            # Selects the row for mangaid from the plateHolesSorted structure
            # if plateHolesSorted is defined.
            if plateHolesSortedStructure is not None:

                if mangaid not in plateHolesSortedStructure['mangaid']:
                    raise GohanPlateTargetsError(
                        'mangaid={0} not found in plateHolesSorted'
                        .format(mangaid))

                plateHolesSortedRow = plateHolesSortedStructure[
                    plateHolesSortedStructure['mangaid'] == mangaid][0]

            catalogueRow = utils.getCatalogueRow(mangaid)
            if catalogueRow is None:
                warnings.warn('no catalogue data found for mangaid={0}'
                              .format(mangaid), GohanPlateTargetsWarning)

            targetData = {}
            for field in fields:

                if field == 'mangaid':
                    targetData[field] = mangaid
                    continue

                if field == 'neverobserve':
                    targetData[field] = 1 if designid in neverobserve else 0
                    continue

                if field == 'nsa_version':
                    targetData[field] = self._getNSAVersion()
                    continue

                # Some particular NSA cases
                if 'nsa' in field:

                    if (field == 'nsa_inclination' and
                            catalogueRow is not None and
                            'ba90' in catalogueRow.colnames):
                        targetData[field] = np.round(
                            np.rad2deg(np.arccos(catalogueRow['ba90'])), 4)
                        continue

                    elif field == 'nsa_id' or field == 'nsa_id100':
                        targetData[field] = mangaInputRow['nsaid']
                        continue

                # Defines the fieldname to search in mangaInput, plateHoles,
                # and catalogue data.

                searchField = None

                if field in conversions['all']:
                    searchField = conversions['all'][field]

                if self.catalogid in conversions:
                    if field in conversions[self.catalogid]:
                        searchField = conversions[self.catalogid][field]

                if searchField is None:
                    if 'nsa_' in field:
                        searchField = field[4:]
                    else:
                        searchField = field

                # Fills values using the mangaInput structure
                if searchField in mangaInputRow.colnames:
                    targetData[field] = mangaInputRow[searchField]
                    continue

                # Fills values using the mangaInput keywords
                if searchField in mangaInputKyw:
                    targetData[field] = mangaInputKyw[searchField]
                    continue

                # Fills values using the catalogue row
                if catalogueRow is not None:
                    if searchField in catalogueRow.colnames:
                        targetData[field] = catalogueRow[searchField]
                        continue

                # Fills values using plateHolesSorted and keywords
                if plateHolesSortedStructure:
                    if searchField in plateHolesSortedRow.colnames:
                        targetData[field] = plateHolesSortedRow[searchField]
                        continue

                    if searchField in plateHolesSortedKyw:
                        targetData[field] = plateHolesSortedKyw[searchField]
                        continue

                # If not found, fills values with -999
                targetData[field] = -999

                # TODO: test with plateTargets-12

            # Cleans up values
            targetData = self._cleanupTargetData(targetData)

            # Applies target fixes
            targetData = self._applyTargetFix(targetData)

            # Adds the new targets
            if not existing:
                self.structure.add_row(targetData)
            else:
                # If the target already exists, replaces it values
                row = self.structure[
                    self.structure['mangaid'] == mangaid &
                    self.structure['plateid'] == plateid]
                for field in targetData:
                    row[field] = targetData[field]

            self._nAppended += 1
            log.debug(
                'mangaid={0} added to plateTargets-{1}'.format(mangaid,
                                                               self.catalogid))

        return self.structure[-len(mangaids):]

    def _readMangaInput(self, mangaInput):
        """Reads a mangaInput file and returns the structure and keywords."""

        mangaInputData = yanny.yanny(mangaInput, np=True)

        if 'MANGAINPUT' in mangaInputData.keys():
            structName = 'MANGAINPUT'
        elif 'STRUCT1' in mangaInputData.keys():
            structName = 'STRUCT1'
        else:
            raise GohanPlateTargetsError(
                'cannot identify structure name for mangaInput {0}'
                .format(mangaInput))

        if 'manga_tileid' in mangaInputData:
            pass
        elif 'tileid' in mangaInputData:
            mangaInputData['manga_tileid'] = mangaInputData.pop('tileid')
        else:
            raise GohanPlateTargetsError(
                'no maga_tileid or tileid field found in mangaInput {0}'
                .format(mangaInput))

        structure = table.Table(mangaInputData.pop(structName))
        mangaInputData.pop('symbols')

        # return the MANGAINPUT structure and all the keywords.
        return (structure, mangaInputData)

    def _toLowerCase(self, structure):
        """Renames all the columns in an `astropy.table.Table` to lowercase."""

        for colname in structure.colnames:
            if colname != colname.lower():
                structure.rename_column(colname, colname.lowercase())

        return structure

    def _getNSAVersion(self):
        """Returns the version of the NSA catalogue used."""

        if self.catalogid == 1:
            return 'v1_0_0'
        elif self.catalogid == 12:
            return 'v1b_0_0_v2'
        else:
            return '-999'

    def _cleanupTargetData(self, targetData):
        """Does some cleaning on the targetData dictionary."""

        for key in targetData:
            if isinstance(targetData[key], basestring):
                targetData[key] = targetData[key].strip()

        return targetData

    def _applyTargetFix(self, targetData):
        """Reads the targetfix file for a target and updates the
        target data."""

        plateid = targetData['plateid']

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
        mangaid = targetData['mangaid']
        targetFix_mangaid = targetFix[targetFix['mangaid'] == mangaid]

        # If targetfixes are found for our mangaid, proceeds to update
        # targetData.
        if len(targetFix_mangaid) > 0:
            log.important('Applying target fix to mangaid={0}'.format(mangaid))
            for row in targetFix_mangaid:
                keyword = row['keyword']
                if keyword in targetData:
                    oldValue = targetData[keyword]
                    newValue = row['value']
                    targetData[keyword] = newValue
                    log.debug('mangaid={0}: {1}={2} -> {3}'.format(
                              mangaid, keyword, oldValue, newValue))

        return targetData
