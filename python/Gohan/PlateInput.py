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

import copy
import fnmatch
import os
import shutil as sh
from collections import OrderedDict

import numpy as np
import six
from astropy import coordinates as coo
from astropy import table

from Gohan import config, exceptions, log, readPath
from Gohan.utils import (assignIFUDesigns, getCatalogueRow, getMaskBitFromLabel,
                         getPlateDefinition, sortTargets, yanny)


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
        targets the plateInput file will contain. If `targettype='sky'`, a
        keyword `mangaInputs=` is required and must be equal to a list of all
        the `PlateInput` objects for the design (normally science and
        standard). These inputs are used to determine the position of the
        targets and standards to constrain the sky catalogue.
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
    decollideExternal : bool, optional
        If False, skips the decollisions with other targets.
    decollidePlateInputs : list, optional
        A list of `PlateInput` instances that have higher priority than then
        current PlateInput. Their targets will be used for decollision.
    FOV : float, optional
        The radius around ``(raCen, decCen)``, in degrees, for which we can allocate
        targets.
    sort : bool, optional
        If True, sorts the targets so that they are evenly distributed on the
        field. Default is False. The number of targets that will be selected to
        be sorted is defined in the `defaults` file.
    plotSorted : bool, optional
        If `sort=True` and `plotSorted=True` (default), a plot of the selected
        sorted targets vs all targets is created.
    rejectTargets : list, optional
        A list of mangaids of targets to be rejected.
    plotIFUs : bool, optional
        If True, a S with the position of each IFU on the plate is saved.
    silentOnCollision : bool, optional
        If False, does not raise a warning on collision. Default is False.
    autocomplete : bool, optional
        If True (the deafult) and `targettype='science'`, autocompletes
        unallocated bundles using NSA targets.
    fillFromCatalogue : bool, optional
        If True (default) and some of the targets contain fields with -999
        values, tries to replace them with values from the parent catalogue
        for that target.
    exclude_ifudesigns: list
        A list of ifudesigns that will not be allocated.

    """

    def __init__(self, designid, targettype, plateRun=None, catalogues=None,
                 surveyMode='mangaLead', raCen=None, decCen=None,
                 locationid=None, manga_tileid=None, fieldName=None,
                 rejectTargets=[], plotIFUs=False,
                 fillFromCatalogue=True, exclude_ifudesigns=[], **kwargs):

        assert isinstance(rejectTargets, (list, tuple, np.ndarray))
        assert surveyMode in [None, 'mangaLead', 'apogeeLead', 'mangaOnly', 'MaStar']

        self.designid = designid
        self.targettype = targettype
        self.plateRun = plateRun
        self.surveyMode = surveyMode
        self.locationid = locationid
        self.fieldName = fieldName
        self.manga_tileid = manga_tileid

        self.mangaInput = None

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
                rejectTargets=rejectTargets,
                exclude_ifudesigns=exclude_ifudesigns, **kwargs)

        else:
            if raCen is None or decCen is None:
                raise exceptions.GohanPlateInputError(
                    'if surveyMode=mangaLed or mangaOnly, raCen and decCen '
                    'need to be defined.')

            autocomplete = kwargs.pop('autocomplete', False)

            if surveyMode is None:
                # Right now this survey mode is used for
                self.manga_tileid = self.manga_tileid or -999
                self.locationid = self.locationid or -999
                check_min_targets = False
                autocomplete = False

            else:
                assert self.manga_tileid is not None
                assert self.locationid is not None
                check_min_targets = True

            self.raCen, self.decCen = raCen, decCen

            if fieldName is None:
                self.fieldName = 'MJ{0:.5f}{1:+.5f}'.format(
                    self.raCen, self.decCen)

            targets = self._getMaNGALeadTargets(
                catalogues, raCen, decCen,
                rejectTargets=rejectTargets,
                check_min_targets=check_min_targets,
                autocomplete=autocomplete,
                **kwargs)

        if targets is None:
            return

        targets = self._tidyUpTargets(targets,
                                      fillFromCatalogue=fillFromCatalogue)

        if self.targettype != 'sky':
            log.debug('reassigning IFUs ... ')
            check_min_targets = True if surveyMode is not None else False
            targets = assignIFUDesigns(targets, (self.raCen, self.decCen),
                                       targettype=self.targettype,
                                       plot=plotIFUs,
                                       plotFilename='ifuPlot_{0}.pdf'.format(self.designid),
                                       check_min_targets=check_min_targets,
                                       exclude_ifudesigns=exclude_ifudesigns)

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

        if self.mangaInput is None:
            return []

        mangaids = self.mangaInput['mangaid']

        return [mangaid.strip() for mangaid in mangaids]

    def copy(self):
        """Returns a copy of self."""

        return copy.copy(self)

    def _selectTargets(self, data):
        """Selects targets from a catalogue that are within the FOV."""

        nData = len(data)

        coords = coo.SkyCoord(data['RA'], data['DEC'], unit='deg')
        separation = coords.separation(
            coo.SkyCoord(ra=self.raCen, dec=self.decCen, unit='deg')).deg

        data = data[np.where(separation <= config['decollision']['FOV'])]

        if len(data) < nData:
            log.warning('{0} targets rejected because '
                        'they are out of the FOV of the design.'.format(nData - len(data)),
                        exceptions.GohanUserWarning)

        return data

    def _formatCatalogue(self, catalogue):
        """Makes sure that the catalogue format is uniform. Reads it if
        `catalogue` is a path"""

        if isinstance(catalogue, six.string_types):
            data = table.Table.read(catalogue, format='fits')
        elif isinstance(catalogue, table.Table):
            data = catalogue
        else:
            raise exceptions.GohanPlateInputError('catalogue has wrong format')

        for col in data.colnames:
            if col != col.upper():
                data.rename_column(col, col.upper())
            if data[col].dtype.type == np.string_:
                data[col] = list(map(lambda xx: xx.strip(), data[col]))

        if 'MANGAID' not in data.colnames and self.targettype == 'sky':
            mangaidCol = table.Column(['0-{0:d}'.format(targetid)
                                       for targetid in range(len(data))],
                                      name='MANGAID', dtype='S20')
            data.add_column(mangaidCol, 0)

        # If the data does not have RA/DEC columns, we add them by copying the
        # IFU_RA/DEC ones.
        for col in ['RA', 'DEC']:
            if col not in data.colnames:
                if 'IFU_' + col in data.colnames:
                    ifuCol = table.Column(data=data['IFU_' + col].data,
                                          name=col)
                    data.add_column(ifuCol)
                elif 'IFU' + col in data.colnames:
                    ifuCol = table.Column(data=data['IFU' + col].data,
                                          name=col)
                    data.add_column(ifuCol)
                else:
                    raise exceptions.GohanPlateInputError(
                        'input catalogue does not contain a column {0} or '
                        'IFU_{0} or IFU{0}'.format(col))

        return data

    def _selectSkies(self, data, skyPatrolRadius=16 / 60.,
                     minNeightborDist=4, **kwargs):
        """Selects skies near the targets defined in the list of `mangaInputs`
        PlateInput objects."""

        for col in data.colnames:
            if col != col.upper():
                data.rename_column(col, col.upper())

        try:
            mangaInputs = kwargs.get('mangaInputs')
        except:
            raise exceptions.GohanPlateInputError(
                'if targettype=sky, you must call PlateInput with a keyword '
                'mangaInputs.')

        # Selects only skies within the FOV
        catCoords = coo.SkyCoord(data['RA'], data['DEC'], unit='deg')
        centre = coo.SkyCoord(ra=self.raCen, dec=self.decCen, unit='deg')
        separationCentre = catCoords.separation(centre).deg

        data = data[np.where(separationCentre <= config['decollision']['FOV'])]
        data = data[data['NEIGHBOR_DIST'] <= 60]
        data.sort('NEIGHBOR_DIST')
        data.reverse()

        # Creates a list of skies within the fibre patrol radius of each target
        # in mangaInputs.
        validSkies = []
        maxSkies = config['plateInputs']['maxSkiesPerTarget']
        for mangaInput in mangaInputs:

            targetsRADec = mangaInput.getTargetCoords()
            mangaIDs = mangaInput.getMangaIDs()

            for ii, targetRADec in enumerate(targetsRADec):

                ra, dec = targetRADec

                catCoords = coo.SkyCoord(data['RA'], data['DEC'], unit='deg')
                targetCoords = coo.SkyCoord(ra, dec, unit='deg')
                sep = catCoords.separation(targetCoords).deg

                indices = np.where(
                    (sep <= skyPatrolRadius) &
                    (data['NEIGHBOR_DIST'] > minNeightborDist))[0]
                valid = data[indices]

                # Selects a group of the skies with largest neighbour distance
                maxSkies = len(valid) if len(valid) < maxSkies else maxSkies
                validSkies.append(valid[0:maxSkies])

                if maxSkies < 10:
                    log.warning('only {0} skies around target {1}. '
                                'Try reducing minNeightborDist.'.format(maxSkies, mangaIDs[ii]),
                                exceptions.GohanUserWarning)

                data.remove_rows(indices[0:maxSkies])

        return table.vstack(validSkies)

    def _parseCatalogues(self, catalogues, leadSurvey='manga'):
        """Returns a list of astropy tables with the catalogues to be used."""

        assert leadSurvey in ['apogee', 'manga'], 'leadSurvey must be apogee or manga'

        if not isinstance(catalogues, six.string_types) and len(catalogues) == 0:
            catalogues = None

        if catalogues is None:

            if self.targettype == 'science':
                if leadSurvey == 'apogee':
                    catalogues = [readPath(config['catalogues']['stelLib'])]
                else:
                    catalogues = [readPath(config['catalogues']['science'])]
                    log.warning('using parent science catalogue', exceptions.GohanUserWarning)

            elif self.targettype == 'standard':
                if leadSurvey == 'apogee':
                    catalogues = [readPath(config['catalogues']['standard']),
                                  readPath(config['catalogues']['APASS'])]
                else:
                    catalogues = readPath(config['catalogues']['standard'])
                    log.warning('using parent standard catalogue', exceptions.GohanUserWarning)

            else:
                raise exceptions.GohanPlateInputError(
                    'no default catalogue for skies.')

        else:
            if isinstance(catalogues, (six.string_types, table.Table)):
                catalogues = [catalogues]
            for ii in range(len(catalogues)):
                if isinstance(catalogues[ii], six.string_types):
                    catalogues[ii] = readPath(catalogues[ii])

        cataloguePaths = [cat for cat in catalogues
                          if isinstance(cat, six.string_types)]
        catalogueTables = [cat for cat in catalogues
                           if isinstance(cat, table.Table)]

        if len(cataloguePaths) > 0:
            log.debug('catalogue paths are {0}'.format(str(cataloguePaths)))
        if len(catalogueTables):
            log.debug('catalogue tables: {0}'.format(len(catalogueTables)))

        for ii in range(len(catalogues)):
            if isinstance(catalogues[ii], six.string_types):
                assert os.path.exists(catalogues[ii]), \
                    'catalogue {0} does not exist'.format(catalogues[ii])
                catalogues[ii] = table.Table.read(catalogues[ii],
                                                  format='fits')

        return catalogues

    def _getInfoFromAPOGEE(self):
        """Reads the APOGEE plateInput files and returns target and design
        information."""

        inputsPath = os.path.join(readPath(config['platelist']), 'inputs')

        plateDefinitionPath = getPlateDefinition(self.designid)
        if not os.path.exists(plateDefinitionPath):
            raise ValueError('plateDefinition {0} not found'.format(plateDefinitionPath))

        plateDefinition = yanny.yanny(plateDefinitionPath, np=True)
        nInputs = int(plateDefinition['nInputs'])

        scienceInput = None
        stdInput = None

        for nn in range(nInputs):
            inp = plateDefinition['plateInput{0}'.format(nn + 1)]
            if fnmatch.fnmatch(inp, '*plateInput*_SCI_*.par'):
                scienceInput = os.path.join(inputsPath, inp)
            elif fnmatch.fnmatch(inp, '*plateInput*_STA_*.par'):
                stdInput = os.path.join(inputsPath, inp)

        assert scienceInput is not None and stdInput is not None, \
            'science or standard input cannot be found in plateDefinition.'
        assert all(map(os.path.exists, [scienceInput, stdInput])), \
            'one or more of the plateInput paths in plateDefinition does ' + \
            'not exist.'

        if not os.path.exists(scienceInput) or not os.path.exists(stdInput):
            raise ValueError('cannot find APOGEE science or standard input file.')

        sciData = yanny.yanny(scienceInput, np=True)
        stdData = yanny.yanny(stdInput, np=True)

        locationid = int(sciData['locationid'])
        raCen = float(sciData['raCen'])
        decCen = float(sciData['decCen'])

        sciCoords = sciData['APOGEEINPUT1'][['ra', 'dec']]
        stdCoords = stdData['APOGEEINPUT2'][['ra', 'dec']]

        apogeeCoords = np.concatenate((sciCoords, stdCoords))
        apogeeCoords = np.array(apogeeCoords.tolist())

        fieldParts = []
        fieldSplits = scienceInput.split('/')[-1].split('_')
        fieldParts = [fieldSplits[ii]
                      for ii in range(1, len(fieldSplits))
                      if '.par' not in fieldSplits[ii]]
        self.fieldName = '_'.join(fieldParts[0:-1])

        return raCen, decCen, locationid, apogeeCoords

    def _getAPOGEELeadTargets(self, catalogues, raCen=None, decCen=None,
                              rejectTargets=[], decollide=True,
                              decollideExternal=True, sort=False,
                              exclude_ifudesigns=[], **kwargs):
        """Creates an APOGEE led plateInput file."""

        self.catalogues = self._parseCatalogues(catalogues,
                                                leadSurvey='apogee')

        self.raCen, self.decCen, self.locationid, \
            apogeeCoords = self._getInfoFromAPOGEE()

        log.debug('raCen={0:.4f}, decCen={1:.4f}'.format(
                  self.raCen, self.decCen))
        log.debug('locationid={0}'.format(self.locationid))

        n_exclude = len(exclude_ifudesigns)

        if self.targettype == 'science':
            nBundlesToAllocate = sum(config['IFUs'].values()) - n_exclude
        elif self.targettype == 'standard':
            nBundlesToAllocate = sum(config['miniBundles'].values()) - n_exclude
        elif self.targettype == 'sky':
            nBundlesToAllocate = sum(config['skies'].values()) - n_exclude

        targetCats = []

        coords = apogeeCoords

        nAllocated = 0

        for catalogue in self.catalogues:

            catalogue = self._formatCatalogue(catalogue)

            if self.targettype != 'sky':
                targetsInField = self._selectTargets(catalogue)
                targetsInField = self._rejectTargets(targetsInField,
                                                     rejectTargets)
            else:
                # Sky selection is different, so we call a special method.
                targetsInField = catalogue
                # self._selectSkies(catalogue, **kwargs)

            if isinstance(catalogue, six.string_types):
                log.debug('{0} targets selected from catalogue {1}'
                          .format(len(targetsInField), catalogue))
            else:
                log.debug('{0} targets selected from table'
                          .format(len(targetsInField)))

            if len(targetsInField) == 0:
                continue

            # We don't decollide for skies
            if self.targettype == 'sky':
                log.warning('Skipping decollision because this is a sky '
                            'catalogue.', exceptions.GohanUserWarning)
                targetCats.append(targetsInField)
                nAllocated += len(targetsInField)
                continue

            if decollide:
                targetsInCat = self._decollide(
                    targetsInField, coords=coords if decollideExternal else None,
                    targetAvoid=config['decollision']['targetAvoidAPOGEE'], **kwargs)
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
            raise exceptions.GohanPlateInputError('not enough targets.')

        targets = self._combineTargetCatalogues(targetCats)
        if targets is None:
            return None

        if self.targettype == 'sky':
            limitTargets = 1e4  # Maximum number of skies per field
        else:
            limitTargets = nBundlesToAllocate

        if sort:
            targetCoords = np.zeros((len(targets), 2), np.float64)
            targetCoords[:, 0] = targets['RA']
            targetCoords[:, 1] = targets['DEC']
            log.debug('sorting targets')
            newCoords, order = sortTargets(
                targetCoords, (self.raCen, self.decCen),
                plot=kwargs.get('plotSorted', True),
                limitTo=int(limitTargets),
                filename=('sortedTargets_{0:d}_{1}.png'.format(self.designid,
                                                               self.targettype[0:3].upper())))
            targets = targets[order]
        else:
            targets = targets[0:int(limitTargets)]

        for col in targets.colnames:
            if col.lower() == 'priority':
                targets[col] = np.arange(len(targets), dtype=int) + 1
                break

        return targets

    def _getMaNGALeadTargets(self, catalogues, raCen=None, decCen=None,
                             rejectTargets=[], decollide=True, sort=False,
                             autocomplete=True, check_min_targets=True, **kwargs):
        """Creates a MaNGA led plateInput file."""

        self.catalogues = self._parseCatalogues(catalogues, leadSurvey='manga')

        log.debug('raCen={0:.4f}, decCen={1:.4f}'.format(
                  self.raCen, self.decCen))
        log.debug('locationid={0}'.format(self.locationid))

        if self.targettype == 'science':
            nBundlesToAllocate = np.sum(list(config['IFUs'].values()))
        elif self.targettype == 'standard':
            nBundlesToAllocate = np.sum(list(config['miniBundles'].values()))
        elif self.targettype == 'sky':
            nBundlesToAllocate = np.sum(list(config['skies'].values()))

        targetCats = []
        coords = None
        nAllocated = 0
        for catalogue in self.catalogues:
            catalogue = self._formatCatalogue(catalogue)
            targetsInField = self._selectTargets(catalogue)
            targetsInField = self._rejectTargets(targetsInField, rejectTargets)

            if isinstance(catalogue, six.string_types):
                log.debug('{0} targets selected from catalogue {1}'
                          .format(len(targetsInField), catalogue))
            else:
                log.debug('{0} targets selected from table'
                          .format(len(targetsInField)))

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
        if targets is None:
            return None

        if nAllocated < nBundlesToAllocate:
            if self.targettype == 'science' and autocomplete:
                # Imports autocomplete here to avoid loading the NSA catalogue
                # if not needed.
                from Gohan.utils.autocomplete import autocomplete
                # Autocompletes targets
                targets = autocomplete(targets, (self.raCen, self.decCen),
                                       **kwargs)
            if check_min_targets and len(targets) < nBundlesToAllocate:
                raise exceptions.GohanPlateInputError(
                    'not enough targets.')

        if sort:
            targetCoords = np.zeros((len(targets), 2), np.float64)
            targetCoords[:, 0] = targets['RA']
            targetCoords[:, 1] = targets['DEC']
            log.debug('sorting targets')
            newCoords, order = sortTargets(
                targetCoords, (self.raCen, self.decCen), plot=False,
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

        if len(targetCats) == 0:
            return None

        if len(targetCats) == 1:
            return targetCats[0]

        targets = targetCats[0]
        for targetCat in targetCats[1:]:
            targets = table.vstack([targets, targetCat])

        log.debug('vstacked {0} catalogues'.format(len(targetCats)))

        return targets

    def _decollide(self, targets, coords=None, targetAvoid=None,
                   decollidePlateInputs=[], **kwargs):
        """Rejects targets that collided either with coords or with targets
        in `decollidePlateInputs`."""

        nTargets = len(targets)

        if coords is not None or len(decollidePlateInputs) == 0:

            log.debug('Decolliding against other catalogues.')

            if self.targettype == 'sky':
                # If this is a sky plateInput, we must have defined a
                # mangaInputs keyword, and those plateInputs are the ones
                # against which we want to decollide.
                mangaInputs = kwargs.get('mangaInputs', [])
                for mangaInput in mangaInputs:
                    if mangaInput not in decollidePlateInputs:
                        decollidePlateInputs.append(mangaInput)
                        log.debug('added mangaInput to decollidePlateInputs')

            if len(decollidePlateInputs) > 0:
                plateInputCoords = np.concatenate(
                    [plate.getTargetCoords()
                     for plate in decollidePlateInputs], axis=0)
                targets = self.decollideCoords(targets, plateInputCoords,
                                               **kwargs)
            if coords is not None:
                targets = self.decollideCoords(targets, coords,
                                               targetAvoid=targetAvoid, **kwargs)
            log.debug('{0} targets rejected from collision with other '
                      'catalogues.'.format(nTargets - len(targets)))

        nTargets = len(targets)
        targets = self.internalDecollision(targets, **kwargs)
        log.debug('{0} targets rejected because internal collisions.'
                  .format(nTargets - len(targets)))

        return targets

    def internalDecollision(self, targets, FOV=None, **kwargs):
        """Decollides targets against themselves. `targets` must be an
        `astropy.table.Table` with `RA` and `DEC` and `MANGAID` fields."""

        silent = kwargs.get('silentOnCollision', False)

        centreAvoid = config['decollision']['centreAvoid']
        FOV = config['decollision']['FOV'] if FOV is None else FOV
        targetAvoid = config['decollision']['targetAvoid']

        centralPost = coo.SkyCoord(self.raCen, self.decCen, unit='deg')
        coords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        separationToPost = centralPost.separation(coords).deg

        centralPostCollisions = np.where(separationToPost < centreAvoid)[0]
        for ii in centralPostCollisions:
            self.logCollision(
                'mangaid={0} '.format(targets[ii]['MANGAID'].strip()) +
                'rejected: collides with central post', silent=silent)

        outOfField = np.where(separationToPost > FOV)[0]
        for ii in outOfField:
            self.logCollision(
                'mangaid={0} '.format(targets[ii]['MANGAID'].strip()) +
                'rejected: outside FOV', silent=silent)
        targets.remove_rows(
            np.concatenate((centralPostCollisions, outOfField)))

        collisions = []
        collided = []
        coords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        for ii in range(len(coords) - 1, -1, -1):
            separations = coords[ii].separation(coords[0:ii]).deg

            if any(separations < targetAvoid):
                collisions.append(ii)
                collided.append(
                    (separations < targetAvoid).tolist().index(True))

        collisionPairs = zip(collisions, collided)
        for ii, jj in collisionPairs:
            self.logCollision(
                'mangaid={0} rejected: internal collision with {1}. '
                'Separation={2:.2f} arcsec'
                .format(targets[ii]['MANGAID'].strip(),
                        targets[jj]['MANGAID'].strip(),
                        coords[ii].separation(coords[jj]).deg * 3600,
                        silent=silent))

        if len(collisions) > 0:
            targets.remove_rows(collisions)

        return targets

    def decollideCoords(self, targets, coords, targetAvoid=None, **kwargs):
        """Decollides a list of targets against a list of coordinates.
        `coords` must be a np.ndarray of Nx2 dimensions. `targets` must be an
        `astropy.table.Table` with `RA` and `DEC` and `MANGAID` fields."""

        silent = kwargs.get('silentOnCollision', False)

        targetAvoid = config['decollision']['targetAvoid'] if targetAvoid is None else targetAvoid

        skyCoords = coo.SkyCoord(coords, unit='deg')
        targetCoords = coo.SkyCoord(targets['RA'], targets['DEC'], unit='deg')

        targetsToReject = []
        for ii in range(len(targetCoords)):

            separations = targetCoords[ii].separation(skyCoords).deg

            if np.any(separations < targetAvoid):
                targetsToReject.append(ii)
                idx_min_sep = np.argmin(separations)
                self.logCollision(
                    'mangaid={0} '.format(targets['MANGAID'][ii].strip()) +
                    'rejected: collides with ' +
                    'RA={0:.5f}, Dec={1:.5f}. Separation={2:.1f} arcsec.'.format(
                        skyCoords[idx_min_sep].ra.deg,
                        skyCoords[idx_min_sep].dec.deg,
                        separations[idx_min_sep] * 3600.), silent=silent)

        targets.remove_rows(targetsToReject)

        return targets

    def logCollision(self, message, silent=False):
        """Raises a warning due to collision."""
        if not silent:
            log.warning(message, exceptions.GohanCollisionWarning)
        else:
            log.debug(message)

    def _tidyUpTargets(self, targets, fillFromCatalogue=True):
        """Makes sure all mandatory columns are present, adds necessary but
        autofillable columns and tidies up the target list."""

        # Makes column names low case and strips strings.
        for colname in targets.colnames:
            if targets[colname].dtype.kind in ['U', 'S']:
                targets[colname] = list(map(lambda xx: xx.strip(), targets[colname]))
            if colname.lower() != colname:
                targets.rename_column(colname, colname.lower())

        # Defines columns that are mandatory and others that can be filled
        # automatically
        mandatoryColumns = ['mangaid', 'ra', 'dec']

        fillableColumns = ['sourcetype', 'manga_target1', 'manga_target2',
                           'manga_target3', 'priority']

        # Adds some more columns depending on the targettype
        if self.targettype != 'sky':
            fillableColumns = ['psfmag', 'ifudesign'] + fillableColumns

        if (self.targettype == 'science' and
                self.surveyMode in ['mangaLead', 'mangaOnly']):
            mandatoryColumns.append('ifudesignsize')
        elif self.targettype == 'sky':
            pass
        else:
            fillableColumns.insert(2, 'ifudesignsize')

        # Checks that all mandatory columns are present
        for colname in mandatoryColumns:
            if colname not in targets.colnames:
                raise exceptions.GohanPlateInputError(
                    'mandatory column {0} not present'.format(colname))

        # Fills out auto fillable columns
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

            elif column == 'manga_target2' and self.targettype == 'sky':
                data = np.zeros(len(targets), dtype=int) + 2

            elif column in ['manga_target1', 'manga_target2', 'manga_target3']:
                data = [0] * len(targets)

            elif column == 'priority':
                data = np.arange(len(targets), dtype=int) + 1

            else:
                raise exceptions.GohanPlateInputError(
                    'failed when trying to autocomplete column ' + column)

            targets.add_column(table.Column(data, column))
            log.debug('Automatically added column {0}.'.format(column))

        # If these are skies, ends here.
        if self.targettype == 'sky':
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

        # Reformats manga_targetX to integer
        targets = reformatAstropyColumn(
            targets, ['manga_target1', 'manga_target2', 'manga_target3'],
            [int, int, int])

        # Makes sure that the maskbits are correctly defined
        for target in targets:
            if (target['manga_target1'] == 0 and
                    target['manga_target2'] == 0 and
                    target['manga_target3'] == 0):
                # If manga_target1 = manga_target3 = 0, this is a filler target
                fillerBit = 2**getMaskBitFromLabel('MANGA_TARGET1',
                                                   'FILLER')[0]
                target['manga_target1'] = fillerBit

        # Checks that the proper motions, if present, contain all the needed
        # columns
        colNames = target.colnames
        if (('pmra' in colNames and 'pmdec' not in colNames) or
                ('pmdec' in colNames and 'pmra' not in colNames)):
            raise exceptions.GohanPlateInputError(
                'for proper motions, both pmra and pmdec must be '
                'present')
        if ('pmra' in colNames and 'pmdec' in colNames and
                'epoch' not in colNames):
            raise exceptions.GohanPlateInputError(
                'pmra and pmdec are present but missing epoch column.')

        # Checks ancillary targets for incomplete data
        for target in targets:
            if target['manga_target3'] > 0:
                target = self.fillFromCatalogue(target)

        # Returns a reordered table
        return self.reorder(targets, mandatoryColumns + fillableColumns)

    def fillFromCatalogue(self, target):
        """Checks a target for fields with -999 values and tries to complete
        them from the parent catalogue."""

        mangaid = target['mangaid']

        # First checks if the target needs to be filled with extra data
        fill = False
        for col in target.colnames:
            if (np.isscalar(target[col]) and
                    (target[col] == -999 or
                     (isinstance(target[col], six.string_types) and
                      '-999' in target[col]))):
                fill = True
                break
            elif np.any(np.array(target[col] == -999)):
                fill = True
                break

        if not fill:
            return target

        log.info('filling out incomplete information for mangaid={0}'
                 .format(mangaid))

        # We assume that this method is called after renaming the columns
        # to lowercase.
        catalogueData = table.Table(getCatalogueRow(mangaid))

        # Lowercases the catalogue columns
        for col in catalogueData.colnames:
            if col != col.lower():
                catalogueData.rename_column(col, col.lower())

        # Wherever necessary, replaces the information in the target
        # with the one in the catalogue
        for col in target.colnames:
            if col not in catalogueData.colnames:
                continue
            if (np.isscalar(target[col]) and
                    (target[col] == -999 or
                     (isinstance(target[col], six.string_types) and
                      '-999' in target[col]))):
                target[col] = catalogueData[col][0]
            elif np.any(np.array(target[col] == -999)):
                target[col] = catalogueData[col]

        return target

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

        template = yanny.yanny(readPath('+etc/mangaInput_Default.par'),
                               np=True)

        template['locationid'] = self.locationid
        template['instrument'] = 'MANGA' if self.targettype != 'sky' \
            else 'MANGA_SINGLE'
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

    def merge(self, plateInput):
        """Merges the targets from a different plateInput."""

        assert plateInput.designid == self.designid, 'inconsistent designids'
        assert plateInput.targettype == self.targettype, 'inconsistent targettype'

        self.mangaInput = table.vstack((self.mangaInput, plateInput.mangaInput))

        return self

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
