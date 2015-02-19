#!/usr/bin/env python
# encoding: utf-8
"""
PlateInput.py

Created by José Sánchez-Gallego on 6 Feb 2014.
Licensed under a 3-clause BSD license.

Major revision history:
    3 Feb 2014 J. Sánchez-Gallego
      Initial version.
    13 Feb 2014 J. Sánchez-Gallego
      Modified to use InputCatalogue class as input method.
    28 Feb 2014 J. Sánchez-Gallego
      Added documentation.
    22 May 2014 J. Sánchez-Gallego
      Partially rewritten to work with the changes to InputCatalogue.
    9 Nov 2014 J. Sánchez-Gallego
      Extense rewrite. Now it does not work with an InputCatalogue instance.

"""

import numpy as np
import os
import shutil as sh
import glob
import warnings
from collections import OrderedDict

from astropy import coordinates as coo
from astropy import table

from Gohan import config, log, readPath
from Gohan import exceptions
from Gohan.utils import yanny, sortTargets
from Gohan.utils import assignIFUDesigns
from Gohan.utils import autocomplete


def reformatAstropyColumn(inputTable, columnName, newFormat):
    """Changes the format of a column for an `astropy.table.Column`.

    Returns a new table in which the desired column has been cast to a
    new dtype format.

    Parameters
    ----------
    inputTable : `astropy.table.Table`
        The input table whose column format wants to be changed.

    columnName : string or list of strings
        The name or names (as a list) of the columns to modify.

    newFormat : dtype or list of dtypes
        The new dtype format of the column. If a list, it must have the same
        size as `columnName`.

    Returns
    -------
    outputTable : `astropy.table.Table`
        A new `astropy.table.Table`, identical to `inputTable` but with the
        desired column(s) cast to a new format.

    Example
    -------
      >> table1 = table.Table([[1,2,3],[4,5,6]], names=['columnA', 'columnB'])
      >> table2 = reformatAstropyColumn(table1, 'columnB', 'S5')

    """

    if isinstance(columnName, (list, tuple, np.ndarray)):
        assert isinstance(newFormat, (list, tuple, np.ndarray))
        assert len(columnName) == len(newFormat)
    else:
        columnName = [columnName]
        newFormat = [newFormat]

    dtypes = [inputTable[col].dtype for col in inputTable.colnames]
    colIndices = [inputTable.colnames.index(colName) for colName in columnName]

    newDtype = dtypes
    for ii, colIndex in enumerate(colIndices):
        newDtype[colIndex] = newFormat[ii]

    return table.Table(inputTable, dtype=newDtype)


class PlateInput(object):
    """A class to construct plateInput files.

    To be written.

    Parameters
    ----------
    designid : int
        The designID for the plate.
    targettype : str
        Either `'science'`, `'standard'` or `'sky'` depending on the type of
        targets the plateInput file will contain.
    plateRun : str, optional
        A string with the plateRun (e.g., '2014.02.x.manga').
    catalogues : file or list of files, optional
        The input catalogue data. If None, the default value based on the
        targettype, will be used. It can be a list of catalogues, and they
        will be used in the input order. Catalogues must be in FITS format.
    surveyMode : str, optional
        Either `'mangaLead'` or `'apogeeLead'`. Defaults to `'mangaLead'`
    raCen, decCen : float, optional
        The coordinates of the centre of the field. Mandatory if the plateInput
        is MaNGA lead. If APOGEE lead, the coordinates are taken from the
        APOGEE plateInput files.
    locationid : int, optional
        The locationid of the field. If `surveyMode='mangaLead'`, `locationid`
        has to be defined. If `'apogeeLead'`, the value is taken from the
        APOGEE plateInput files.
    manga_tileid : int, optional
        Same as `locationid`, but if `surveyMode='apogeeLead'`,
        `manga_tileid=-999`
    fieldName : str, optional
        The `fieldName` value to be used. Defaults to `none`.
    decollide : bool, optional
        If True (the default), decollides the input targets against themselves
        and, if `surveyMode='apogeeLead'`, against APOGEE targets.
    decollidePlateInputs : list, optional
        A list of `PlateInput` instances that have higher priority than then
        current PlateInput. Their targets will be used for decollision.
    sort : bool, optional
        If True, sorts the targets so that they are evenly distributed on the
        field. Default is False. The number of targets that will be selected to
        be sorted is defined in the `defaults` file.
    rejectTargets : list, optional
        A list of mangaids of targets to be rejected.
    plotIFUs : bool, optional
        If True, a plot with the position of each IFU on the plate is saved.
    silentOnCollision : bool, optional
        If False, does not raise a warning on collision. Default is False.
    autocomplete : bool, optional
        If True (the deafult) and `targettype='science'`, autocompletes
        unallocated bundles using NSA targets.

    """

    def __init__(self, designid, targettype, plateRun=None, catalogues=None,
                 surveyMode='mangaLead', raCen=None, decCen=None,
                 locationid=None, manga_tileid=None, fieldName=None,
                 rejectTargets=[], plotIFUs=False, **kwargs):

        assert isinstance(rejectTargets, (list, tuple, np.ndarray))
        assert surveyMode in ['mangaLead', 'apogeeLead', 'mangaOnly']

        self.designid = designid
        self.targettype = targettype
        self.plateRun = plateRun
        self.surveyMode = surveyMode
        self.locationid = locationid
        self.fieldName = fieldName
        self.manga_tileid = manga_tileid

        log.info('creating {0} plateInput for design {1}'.format(
                 targettype, designid))
        log.debug('surveyMode={0}'.format(surveyMode))

        if surveyMode == 'apogeeLead':

            self.manga_tileid = manga_tileid \
                if manga_tileid is not None else -999

            if plateRun is None:
                raise exceptions.GohanPlateInputError(
                    'if surveyMode=apogeeLead, plateRun needs to be defined.')

            targets = self._getAPOGEELeadTargets(
                catalogues, raCen=raCen, decCen=decCen,
                rejectTargets=rejectTargets, **kwargs)

        else:
            if raCen is None or decCen is None:
                raise exceptions.GohanPlateInputError(
                    'if surveyMode=mangaLed or mangaOnly, raCen and decCen '
                    'need to be defined.')

            assert self.manga_tileid is not None
            assert self.locationid is not None

            self.raCen, self.decCen = raCen, decCen

            if fieldName is None:
                self.fieldName = 'MJ{0:.5f}{1:+.5f}'.format(
                    self.raCen, self.decCen)

            targets = self._getMaNGALeadTargets(
                catalogues, raCen, decCen,
                rejectTargets=rejectTargets, **kwargs)

        targets = self._tidyUpTargets(targets)

        if self.targettype != 'sky':
            log.debug('reassigning IFUs ... ')
            targets = assignIFUDesigns(targets, (self.raCen, self.decCen),
                                       targettype=self.targettype,
                                       plot=plotIFUs,
                                       filename='ifuPlot_{0}.pdf'.format(
                                           self.designid))

        self.mangaInput = targets

        # Tweaks the default fill_values
        for column in self.mangaInput.colnames:
            dd = self.mangaInput[column].dtype.type
            if dd == np.string_:
                self.mangaInput[column].fill_value = 'none'
            else:
                self.mangaInput[column].fill_value = -999

    def getTargetCoords(self):
        """Returns a Nx2 array with the coordinates of the targets in this
        plateInput file."""

        coords = np.zeros((len(self.mangaInput), 2), np.float64)
        coords[:, 0] = self.mangaInput['ra']
        coords[:, 1] = self.mangaInput['dec']

        return coords

    def getMangaIDs(self):
        """Returns the mangaids of the targets in this plateInput file."""

        mangaids = self.mangaInput['mangaid']

        return [mangaid.strip() for mangaid in mangaids]

    def _selectTargets(self, catalogue):
        """Selects targets from a catalogue that are within the FOV."""

        data = table.Table.read(catalogue, format='fits')
        coords = coo.SkyCoord(data['RA'], data['DEC'], unit='deg')
        separation = coords.separation(
            coo.SkyCoord(ra=self.raCen, dec=self.decCen, unit='deg')).deg

        return data[np.where(separation <= config['decollision']['FOV'])]

    def _getInfoFromAPOGEE(self):
        """Reads the APOGEE plateInput files and returns target and design
        information."""

        platelist = readPath(config['platelist'])
        apogeeInputs = os.path.join(
            platelist, 'inputs', 'apogee', self.plateRun)

        scienceInput = glob.glob(
            os.path.join(
                apogeeInputs,
                'plateInput_*_SCI_{0:d}.par'.format(self.designid)))

        stdInput = glob.glob(
            os.path.join(
                apogeeInputs,
                'plateInput_*_STA_{0:d}.par'.format(self.designid)))

        if len(scienceInput) != 1 and len(stdInput) != 1:
            raise exceptions.GohanPlateInputError(
                'no APOGEE plateInput found.')

        sciData = yanny.yanny(scienceInput[0], np=True)
        stdData = yanny.yanny(stdInput[0], np=True)

        locationid = int(sciData['locationid'])
        raCen = float(sciData['raCen'])
        decCen = float(sciData['decCen'])

        sciCoords = sciData['APOGEEINPUT1'][['ra', 'dec']]
        stdCoords = stdData['APOGEEINPUT2'][['ra', 'dec']]

        apogeeCoords = np.concatenate((sciCoords, stdCoords))
        apogeeCoords = apogeeCoords.view(np.float).reshape(
            apogeeCoords.shape + (-1,))

        fieldParts = []
        fieldSplits = scienceInput[0].split('/')[-1].split('_')
        fieldParts = [fieldSplits[ii]
                      for ii in range(1, len(fieldSplits))
                      if '.par' not in fieldSplits[ii]]
        self.fieldName = '_'.join(fieldParts[0:-1])

        return raCen, decCen, locationid, apogeeCoords

    def _getAPOGEELeadTargets(self, catalogues, raCen=None, decCen=None,
                              rejectTargets=[], decollide=True, sort=False,
                              **kwargs):
        """Creates an APOGEE led plateInput file."""

        if catalogues is None:
            if self.targettype == 'science':
                self.catalogues = readPath(config['catalogues']['stelLib'])
            elif self.targettype == 'standard':
                self.catalogues = [readPath(config['catalogues']['standard']),
                                   readPath(config['catalogues']['APASS'])]
        else:
            if isinstance(catalogues, basestring):
                catalogues = [catalogues]
            self.catalogues = [readPath(cat) for cat in catalogues]

        log.debug('catalogue paths are {0}'.format(str(self.catalogues)))

        assert all(map(os.path.exists, self.catalogues))

        self.raCen, self.decCen, \
            self.locationid, apogeeCoords = self._getInfoFromAPOGEE()

        log.debug('raCen={0:.4f}, decCen={1:.4f}'.format(
                  self.raCen, self.decCen))
        log.debug('locationid={0}'.format(self.locationid))

        if self.targettype == 'science':
            nBundlesToAllocate = np.sum(config['IFUs'].values())
        elif self.targettype == 'standard':
            nBundlesToAllocate = np.sum(config['miniBundles'].values())
        elif self.targettype == 'sky':
            nBundlesToAllocate = np.sum(config['skies'].values())

        targetCats = []
        coords = apogeeCoords
        nAllocated = 0
        for catalogue in self.catalogues:
            targetsInField = self._selectTargets(catalogue)
            targetsInField = self._rejectTargets(targetsInField, rejectTargets)
            log.debug('{0} targets selected from catalogue {1}'
                      .format(len(targetsInField), catalogue))

            if len(targetsInField) == 0:
                continue

            if decollide:
                targetsInCat = self._decollide(targetsInField, coords=coords,
                                               **kwargs)
                log.debug('{0} targets remaining after decollision.'
                          .format(len(targetsInCat)))
                targetCats.append(targetsInCat)

            nAllocated += len(targetsInCat)

            if nAllocated >= nBundlesToAllocate:
                break
            else:
                if len(targetsInCat) > 0:
                    catCoords = np.zeros((len(targetsInCat), 2), np.float64)
                    catCoords[:, 0] = targetsInCat['RA']
                    catCoords[:, 1] = targetsInCat['DEC']
                    coords = np.concatenate((coords, catCoords), axis=0)

        if nAllocated < nBundlesToAllocate:
            raise exceptions.GohanPlateInputError(
                'not enough targets.')

        targets = self._combineTargetCatalogues(targetCats)

        if sort:
            targetCoords = np.zeros((len(targets), 2), np.float64)
            targetCoords[:, 0] = targets['RA']
            targetCoords[:, 1] = targets['DEC']
            log.debug('sorting targets')
            newCoords, order = sortTargets(
                targetCoords, (self.raCen, self.decCen), plot=True,
                limitTo=nBundlesToAllocate,
                filename=('sortedTargets_{0:d}_{1}.pdf'
                          .format(self.designid, self.targettype[0:3].upper()))
                )
            targets = targets[order]
        else:

            targets = targets[0:nBundlesToAllocate]

        for col in targets.colnames:
            if col.lower() == 'priority':
                targets[col] = np.arange(len(targets), dtype=int) + 1
                break

        return targets

    def _getMaNGALeadTargets(self, catalogues, raCen=None, decCen=None,
                             rejectTargets=[], decollide=True, sort=False,
                             **kwargs):
        """Creates a MaNGA led plateInput file."""

        if catalogues is None:
            if self.targettype == 'science':
                self.catalogues = readPath(config['catalogues']['science'])
            elif self.targettype == 'standard':
                self.catalogues = readPath(config['catalogues']['standard'])
        else:
            if isinstance(catalogues, basestring):
                catalogues = [catalogues]
            self.catalogues = [readPath(cat) for cat in catalogues]

        log.debug('catalogue paths are {0}'.format(str(self.catalogues)))

        assert all(map(os.path.exists, self.catalogues))

        log.debug('raCen={0:.4f}, decCen={1:.4f}'.format(
                  self.raCen, self.decCen))
        log.debug('locationid={0}'.format(self.locationid))

        if self.targettype == 'science':
            nBundlesToAllocate = np.sum(config['IFUs'].values())
        elif self.targettype == 'standard':
            nBundlesToAllocate = np.sum(config['miniBundles'].values())
        elif self.targettype == 'sky':
            nBundlesToAllocate = np.sum(config['skies'].values())

        targetCats = []
        coords = None
        nAllocated = 0
        for catalogue in self.catalogues:
            targetsInField = self._selectTargets(catalogue)
            targetsInField = self._rejectTargets(targetsInField, rejectTargets)
            log.debug('{0} targets selected from catalogue {1}'
                      .format(len(targetsInField), catalogue))

            if len(targetsInField) == 0:
                continue

            if decollide:
                targetsInCat = self._decollide(targetsInField, coords=coords,
                                               **kwargs)
                log.debug('{0} targets remaining after decollision.'
                          .format(len(targetsInCat)))
                targetCats.append(targetsInCat)

            nAllocated += len(targetsInCat)

            if nAllocated >= nBundlesToAllocate:
                break
            else:
                if len(targetsInCat) > 0:
                    catCoords = np.zeros((len(targetsInCat), 2), np.float64)
                    catCoords[:, 0] = targetsInCat['RA']
                    catCoords[:, 1] = targetsInCat['DEC']

                    if coords is None:
                        coords = catCoords
                    else:
                        coords = np.concatenate((coords, catCoords), axis=0)

        targets = self._combineTargetCatalogues(targetCats)

        if nAllocated < nBundlesToAllocate:
            if self.targettype != 'science':
                pass
            else:
                targets = autocomplete(targets, 'science',
                                       (self.raCen, self.decCen), **kwargs)
            if len(targets) < nBundlesToAllocate:
                raise exceptions.GohanPlateInputError(
                    'not enough targets.')

        if sort:
            targetCoords = np.zeros((len(targets), 2), np.float64)
            targetCoords[:, 0] = targets['RA']
            targetCoords[:, 1] = targets['DEC']
            log.debug('sorting targets')
            newCoords, order = sortTargets(
                targetCoords, (self.raCen, self.decCen), plot=True,
                limitTo=nBundlesToAllocate,
                filename=('sortedTargets_{0:d}_{1}.pdf'
                          .format(self.designid, self.targettype[0:3].upper()))
                )
            targets = targets[order]
        else:

            targets = targets[0:nBundlesToAllocate]

        for col in targets.colnames:
            if col.lower() == 'priority':
                targets[col] = np.arange(len(targets), dtype=int) + 1
                break

        return targets

    @staticmethod
    def _rejectTargets(targets, mangaids):
        """Rejects targets if they match a list of mangaids."""

        idxToReject = []
        for ii, target in enumerate(targets):
            if target['MANGAID'].strip() in mangaids:
                idxToReject.append(ii)

        targets.remove_rows(idxToReject)

        if len(idxToReject) > 0:
            log.debug('removed {0} targets from catalogue because they '
                      'where in a manual rejection list'
                      .format(len(idxToReject)))

        return targets

    def _combineTargetCatalogues(self, targetCats):
        """Creates a single table from different catalogues with potentially
        different fields."""

        if len(targetCats) == 1:
            return targetCats[0]

        targets = targetCats[0]
        for targetCat in targetCats[1:]:
            targets = table.vstack([targets, targetCat])

        log.debug('vstacked {0} catalogues'.format(len(targetCats)))

        return targets

    def _decollide(self, targets, coords=None, decollidePlateInputs=[],
                   **kwargs):
        """Rejects targets that collided either with coords or with targets
        in `decollidePlateInputs`."""

        nTargets = len(targets)
        targets = self.internalDecollision(targets, **kwargs)
        log.debug('{0} targets rejected because internal collisions.'
                  .format(nTargets-len(targets)))
        nTargets = len(targets)

        if coords is None and len(decollidePlateInputs) == 0:
            return targets
        else:
            log.debug('Decolliding against other catalogues.')
            if len(decollidePlateInputs) > 0:
                plateInputCoords = np.concatenate(
                    [plate.getTargetCoords()
                     for plate in decollidePlateInputs], axis=0)
                targets = self.decollideCoords(targets, plateInputCoords,
                                               **kwargs)
            if coords is not None:
                targets = self.decollideCoords(targets, coords, **kwargs)
            log.debug('{0} targets rejected from collision with other '
                      'catalogues.'.format(nTargets-len(targets)))

        return targets

    def internalDecollision(self, targets, **kwargs):
        """Decollides targets against themselves. `targets` must be an
        `astropy.table.Table` with `RA` and `DEC` and `MANGAID` fields."""

        silent = kwargs.get('silentOnCollision', False)

        centreAvoid = config['decollision']['centreAvoid']
        FOV = config['decollision']['FOV']
        targetAvoid = config['decollision']['targetAvoid']

        centralPost = coo.SkyCoord(self.raCen, self.decCen, unit='deg')
        coords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        separationToPost = centralPost.separation(coords).deg

        centralPostCollisions = np.where(separationToPost < centreAvoid)[0]
        for ii in centralPostCollisions:
            self.logCollision(
                'mangaid={0} '.format(targets[ii]['MANGAID']) +
                'rejected: collides with central post', silent=silent)

        outOfField = np.where(separationToPost > FOV)[0]
        for ii in centralPostCollisions:
            self.logCollision(
                'mangaid={0} '.format(targets[ii]['MANGAID']) +
                'rejected: outside FOV', silent=silent)
        targets.remove_rows(
            np.concatenate((centralPostCollisions, outOfField)))

        collisions = []
        coords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        for ii in range(len(coords)-1, -1, -1):
            separations = coords[ii].separation(coords[0:ii]).deg

            if any((separations < targetAvoid)):
                collisions.append(ii)

        collisions = np.sort(np.unique(collisions))
        for jj in collisions:
            self.logCollision('mangaid={0} rejected: internal collision'
                              .format(targets[jj]['MANGAID'].strip()),
                              silent=silent)

        targets.remove_rows(collisions)

        return targets

    def decollideCoords(self, targets, coords, **kwargs):
        """Decollides a list of targets against a list of coordinates.
        `coords` must be a np.ndarray of Nx2 dimensions. `targets` must be an
        `astropy.table.Table` with `RA` and `DEC` and `MANGAID` fields."""

        silent = kwargs.get('silentOnCollision', False)

        targetAvoid = config['decollision']['targetAvoid']

        skyCoords = coo.SkyCoord(coords, unit='deg')
        targetCoords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        targetsToReject = []
        for ii in range(len(targetCoords)):

            separations = targetCoords[ii].separation(skyCoords).deg

            if np.any(separations < targetAvoid):
                targetsToReject.append(ii)
                collisions = np.where(separations < targetAvoid)[0]
                self.logCollision(
                    'mangaid={0} '.format(targets['MANGAID'][ii]) +
                    'rejected: collides with ' +
                    'RA={0:.5f}, Dec={0:.5f}'.format(
                        skyCoords[collisions[0]].ra.deg,
                        skyCoords[collisions[0]].dec.deg), silent=silent)

        targets.remove_rows(targetsToReject)

        return targets

    def logCollision(self, message, silent=False):
        """Raises a warning due to collision."""
        if not silent:
            warnings.warn(message, exceptions.GohanCollisionWarning)
        else:
            log.debug(message)

    def _tidyUpTargets(self, targets):
        """Makes sure all mandatory columns are present, adds necessary but
        autofillable columns and tidies up the target list."""

        # Makes column names low case.
        for colname in targets.colnames:
            if colname.lower() != colname:
                targets.rename_column(colname, colname.lower())

        mandatoryColumns = ['mangaid', 'ra', 'dec']

        fillableColumns = ['sourcetype', 'manga_target1', 'manga_target2',
                           'manga_target3', 'priority']
        if self.targettype != 'sky':
            fillableColumns = ['psfmag', 'ifudesign'] + fillableColumns

        if (self.targettype == 'science' and
                self.surveyMode in ['mangaLead', 'mangaOnly']):
            mandatoryColumns.append('ifudesignsize')
        elif self.targettype == 'sky':
            pass
        else:
            fillableColumns.insert(2, 'ifudesignsize')

        for colname in mandatoryColumns:
            if colname not in targets.colnames:
                raise exceptions.GohanPlateInputError(
                    'mandatory column {0} not present'.format(colname))

        for column in fillableColumns:
            if column in targets.colnames:
                continue

            elif column == 'psfmag':
                data = [[0.0, 0.0, 0.0, 0.0, 0.0]
                        for ii in range(len(targets))]

            elif column == 'sourcetype':

                if self.targettype == 'science':
                    sourcetype = 'SCI'
                elif self.targettype == 'standard':
                    sourcetype = 'STD'
                else:
                    sourcetype = 'SKY'

                data = [sourcetype] * len(targets)

            elif column in ['ifudesignsize', 'ifudesign']:
                data = [-999] * len(targets)

            elif column in ['manga_target1', 'manga_target2', 'manga_target3']:
                data = [0] * len(targets)

            elif column == 'priority':
                data = np.arange(len(targets), dtype=int) + 1

            else:
                raise exceptions.GohanPlateInputError(
                    'failed when trying to autocomplete column ' + column)

            targets.add_column(table.Column(data, column))
            log.debug('Automatically added column {0}.'.format(column))

        if self.targettype == 'sky':
            # If these are skies, ends here.
            return self.reorder(targets, mandatoryColumns + fillableColumns)

        # Makes sure that the ife_ra, ifu_dec, target_ra and target_dec
        # columns exist.
        for coordType in ['ifu', 'target']:
            for col in ['ra', 'dec']:
                if not '{0}_{1}'.format(coordType, col) in targets.colnames:
                    targets.add_column(
                        table.Column(
                            [-999.] * len(targets),
                            '{0}_{1}'.format(coordType, col),
                            dtype=float))
                    log.debug('Adding {0}_{1} column'.format(coordType, col))

        # Replace ra, dec from ifu_ra, ifu_dec if they are defined.
        for col in ['ra', 'dec']:
            for target in targets:

                targetCol = 'target_{0}'.format(col)
                ifuCol = 'ifu_{0}'.format(col)

                if target[targetCol] < -100:
                    target[targetCol] = target[col]

                if target[ifuCol] < -100:
                    target[ifuCol] = target[col]
                elif target[ifuCol] > -100 and target[ifuCol] != target[col]:
                    target[col] = target[ifuCol]
                    log.important('setting {0} from ifu_{0} for target'.format(
                        col) + ' with mangaid={0}'.format(target['mangaid']))

        targets = reformatAstropyColumn(
            targets, ['manga_target1', 'manga_target2', 'manga_target3'],
            [int, int, int])

        return self.reorder(targets, mandatoryColumns + fillableColumns)

    @staticmethod
    def reorder(tt, orderedColumns):
        """Reorders a table by putting the columns in `orderedColumns`
        first."""

        newOrder = []
        for col in orderedColumns:
            if col in tt.colnames:
                newOrder.append(col)
        for col in tt.colnames:
            if col not in newOrder:
                newOrder.append(col)
        tt = tt[newOrder]
        log.debug('Column order rearranged')

        return tt

    def write(self, toRepo=False):
        """Writes the plateInput file to disk or to the repo."""

        filename = self.getDefaultFilename()

        if os.path.exists(filename):
            os.remove(filename)

        template = yanny.yanny(
            readPath('+etc/mangaInput_Default.par'), np=True)

        template['locationid'] = self.locationid
        template['locationid'] = self.locationid
        template['racen'] = self.raCen
        template['deccen'] = self.decCen
        template['designid'] = self.designid
        template['targettype'] = self.targettype
        template['manga_tileid'] = self.manga_tileid
        template['fieldname'] = self.fieldName

        enums = OrderedDict()
        hdr = OrderedDict()
        structname = config['plateInputs']['mangaInputStructure']

        template['symbols'] = template.dtype_to_struct(
            self.mangaInput.dtype,
            structname=structname,
            enums=enums)

        template[structname.upper()] = self.mangaInput.filled()

        for key in hdr:
            template[key] = hdr[key]

        template.write(filename)

        log.info('{0} saved.'.format(filename))

        if toRepo:
            self._toRepo(filename)

        return filename

    def getDefaultFilename(self):

        filename = 'manga' + self.targettype.capitalize()

        if self.surveyMode == 'mangaOnly' or self.surveyMode == 'mangaLead':
            filename += '_{0:d}_{1:d}.par'.format(self.locationid,
                                                  self.designid)
        else:
            filename += '_{0:s}_{1:d}.par'.format(self.fieldName,
                                                  self.designid)

        return filename

    def _toRepo(self, filename):
        """Copies a file to platelist."""

        inputPath = os.path.join(
            readPath(config['platelist']), 'inputs/manga',
            self.plateRun)

        if not os.path.exists(inputPath):
            os.makedirs(inputPath)

        sh.copy(filename, inputPath)

        log.info(filename + ' copied $PLATELIST/inputs.')
