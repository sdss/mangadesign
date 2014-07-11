#!/usr/bin/env python
# encoding: utf-8
"""
InputCatalogue.py

Created by José Sánchez-Gallego on 6 Feb 2014.
Licensed under a 3-clause BSD license.

Major revision history:
    13 Feb 2014 J. Sánchez-Gallego
      Initial version.
    19 Feb 2014  J. Sánchez-Gallego
      Modified to accept both location ids and catalogues.
    28 Feb 2014 J. Sánchez-Gallego
      Added documentation.
    22 May 2014 J. Sánchez-Gallego
      Heavily modified to use a version of the master catalogue that is closer
      to the final one.

"""

from astropy import table
import os
from astropysics.coords import ICRSCoordinates, separation_matrix
from astropysics.coords import AngularCoordinate
from astropysics.coords.funcs import match_coords
from .utils import Staralt
import warnings
from .exceptions import GohanUserWarning, GohanError, GohanNotImplemented
from .exceptions import GohanCollisionWarning
from . import config
from . import readPath
from sdss.utilities import yanny
import fitsio
import numpy as np
from . import log
from collections import OrderedDict


__ALL__ = ['InputCatalogue']
defaultMaNGAInput = yanny.yanny(
    readPath('+etc/mangaInput_Default.par'), np=True)
defaultColumns = defaultMaNGAInput['MANGAINPUT'].dtype.names

targettypeDic = {'SCI': 'science', 'STD': 'standard',
                 'STA': 'standard', 'SKY': 'sky'}

FOV = AngularCoordinate(config['decollision']['FOV'])
centreAvoid = AngularCoordinate(config['decollision']['centreAvoid'])
targetAvoid = AngularCoordinate(config['decollision']['targetAvoid'])

defaultValues = {
    'ifudesign': -999,
    'ifudesignsize': -999,
    'manga_target1': 0,
    'psfmag': [0.0, 0.0, 0.0, 0.0, 0.0],
    'sourcetype': None
}


class InputCatalogue(object):
    """The base InputCatalogue class.

    This class returns an instance of ``astropy.table.Table`` with the fields
    needed to create plateInput files.

    Parameters
    ----------
    tileID : int or None, optional
        The tile id of the field. If an integer, only the rows matching that
        tileID will be used. If None, all the rows in the master catalogue
        will be used.
    format : string, optional
        The formt of the sample catalogue. It can be `'file'` (default) to
        indicate that the catalogue is a file, or `'database'` if the data is
        to be read from mangaSampleDB.
    type : string, optional
        The type of the targets in the catalogue. It can be 'SCI' (the
        default), 'STD' or 'SKY'.
    file : string or None, optinal
        If `format='file'`, `file` can be set with the path of the catalogue to
        be used. If None, the default catalogue for the specific type, as
        defined in the configuration file, will be used.
    conversions : dict, optional
        A dictionary with conversions for the column names in the catalogue
        file. For example, if the columns for right ascension and declination
        in the catalogue are called `'ra'` and `'dec'`, `conversions` can be
        defined as `conversion={'RA': 'ra', 'Dec': 'dec'}`.
    fill : bool, optional
        If True, each row tries to be matched with its corresponding record
        in the NSA catalogue and the additional columns there are added to the
        InputCatalogue table.
    meta : dict, optional
        A dictionary with meta information that will update values in the meta
        attribute of `table.Table`.
    removeSuperfluous : bool, optional
        If True, all non-necessary columns will be removed from the input
        catalogue.
    verbose : bool, optional
        If False, only prints warnings. If None, only errors are raised.
    decollision : bool, optional
        If None, decollision is completely disabled.
        If False, no decollision is done for the targets in the input
        catalogue except checking that all the targets are within the FOV and
        that there are no collisions with the central post.
        If set to a numpy recarray-like object, the targets of the input
        catalogue that collide will be excluded. The input array must have,
        at least, two columns called 'ra' and 'dec'. Collision distances can be
        adjusted via the configuration file.
    failOnCollision : bool, optional
        If False (the default), a collision will raise a warning and the target
        will be rejected but no error will be raised. If True, the process will
        stop after the warning.
    autocomplete : bool, optional
        For science catalogues, if the number of assigned bundles is lower
        than the number of available bundles, the remaining bundles will be
        assigned to unused targets from the general exNSA catalogue.

    """

    def __init__(self, tileid=None, format='file', type='SCI', file=None,
                 conversions=None, fill=False, meta=None,
                 removeSuperfluous=False, verbose=True,
                 decollision=False, failOnCollision=False,
                 autocomplete=True, **kwargs):

        log.setVerbose(verbose)

        self.tileid = tileid
        self.format = format
        self.type = type.upper()
        self.file = file.lower() if file is not None else None
        self.conversions = {} if conversions is None else conversions
        self.fill = fill
        self._meta = {} if meta is None else meta
        self.removeSuperfluous = removeSuperfluous
        self.kwargs = kwargs

        if self.tileid is None and self.format == 'file' and self.file is None:
            raise GohanError('tileid has to be specified if format=\'file\' '
                             'and file=None.')

        if self.format not in ['file', 'database']:
            raise GohanError('format={0} is not valid.'.format(self.format))

        if self.type not in ['STD', 'SCI', 'SKY']:
            raise GohanError('type={0} is not valid.'.format(self.type))

        if self.type == 'SCI':
            if self.file is None:
                self.file = readPath(config['files']['sciSample'])
            self.nBundles = np.sum(config['IFUs'].values())
        elif self.type == 'STD':
            if self.file is None:
                self.file = readPath(config['files']['stdSample'])
            self.nBundles = np.sum(config['miniBundles'].values())
        elif self.type == 'SKY':
            if self.file is None:
                self.file = readPath(config['files']['skySample'])
            self.nBundles = np.sum(config['skies'].values())

        defaultValues['sourcetype'] = self.type

        if not os.path.exists(self.file):
            raise GohanError('catalogue file {0} does not exist.'.format(
                             self.file))

        if self.fill is True:
            warnings.warn('fill=True not yet implemented. Setting fill=False',
                          GohanUserWarning)

        if self.format == 'database':
            raise GohanNotImplemented('The DB access is not yet implemented.')
        elif self.format == 'file':
            log.info('Creating input catalogue from file {0}'.format(
                os.path.basename(self.file)))
            self._createTableFromFile()

        self._createMeta(self._meta)

        if self.removeSuperfluous:
            self.removeSuperfluous()

        if decollision is not None:
            self.decollision(decollision, failOnCollision=failOnCollision)

        if autocomplete and self.type == 'SCI' and len(self) < self.nBundles:
            self.autocomplete()

        self.cropCatalogue()
        self.checkData()

    def removeSuperfluous(self):
        removedCols = []
        for colname in self.colnames:
            if colname not in defaultColumns:
                self.data.remove_column(colname)
                removedCols.append(colname)
        log.debug('Removed superfluous columns {0}'.format(removedCols))

    def autocomplete(self):

        if self.type != 'SCI':
            raise GohanError('autocomplete can only be used for SCI '
                             'inputCatalogues')

        if len(self) >= self.nBundles:
            log.info('the are at least as many targets as bundles. No need '
                     'to run autocomplete.')
        else:
            nBundlesToAssign = self.nBundles - len(self)
            log.info('autocompleting {0} bundles.'.format(nBundlesToAssign))

        exNSATargets = self.getExNSATargets()
        unassignedBundleSizes = self.getUnassignedBundleSizes()

        for bundleSize in unassignedBundleSizes:

            target = self._getOptimalTarget(exNSATargets, bundleSize)

            if target is not None:
                self.addTarget(target, bundleSize)
            else:
                warnings.warn('no valid replacement found for '
                              'ifudesignsize={0}'.format(int(bundleSize)),
                              GohanUserWarning)

            if len(self) >= self.nBundles:
                break

        if len(self) < self.nBundles:
            log.important('even after autocomplete, {0} '.format(
                self.nBundles-len(self)) + 'bundles are still unassigned.')

    def _getOptimalTarget(self, candidateTargets, bundleSize):

        reffField = config['reffField'].upper()
        if reffField in candidateTargets.colnames:
            return self._getOptimalTargetSize(candidateTargets, bundleSize)
        else:
            for target in candidateTargets:
                if self._checkExNSATarget(target):
                    return target
            return None

    def _getOptimalTargetSize(self, candidateTargets, bundleSize):

        reffField = config['reffField'].upper()
        candidateTargets.sort(reffField)

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
                if self._checkExNSATarget(target):
                    if 'IFUDESIGNSIZE' in target.colnames():
                        target.add_column(
                            table.Column(
                                -999, name='IFUDESIGNSIZE', dtype=int))

                    target['IFUDESIGNSIZE'] = bundleSize

                    return target

        # If it comes to this, it means that there are no target with valid
        # size. We fall back to use whatever target that does not collide.
        log.info('no target of right size ({0})'.format(bundleSize) +
                 ' can be used to autocomplete. '
                 'Trying to find a suitable target.')

        for target in candidateTargets:
                if self._checkExNSATarget(target):
                    return target

        return None

    def getUnassignedBundleSizes(self):

        if self.type == 'SCI':
            bundles = config['IFUs'].copy()
        elif self.type == 'STD':
            bundles = config['miniBundles'].copy()
        elif self.type == 'SKY':
            bundles = config['skies'].copy()

        for target in self:
            if target['ifudesignsize'] > 0:
                bundles[target['ifudesignsize']] -= 1

        unassignedBundleSizes = []
        for size in bundles:
            for nn in range(bundles[size]):
                unassignedBundleSizes.append(size)

        return unassignedBundleSizes

    def getExNSATargets(self):

        exNSACat = fitsio.read(config['files']['sciSample'])

        raCen = self.meta['racen']
        decCen = self.meta['deccen']

        matchA, matchB = match_coords(
            raCen, decCen, exNSACat['RA'], exNSACat['DEC'],
            eps=FOV.degrees, mode='mask')

        validExNSATargets = table.Table(exNSACat[matchB])
        validExNSATargets = validExNSATargets[validExNSATargets[
            'IFUDESIGNSIZE'] <= 0]

        return validExNSATargets

    def _checkExNSATarget(self, target):

        if target['MANGAID'] in self.data['mangaid']:
            return False

        centralPost = ICRSCoordinates(self.data.meta['racen'],
                                      self.data.meta['deccen'])

        targetCoords = ICRSCoordinates(target['RA'], target['DEC'])

        if (targetCoords - centralPost).degrees < centreAvoid.degrees:
            return False
        if (targetCoords - centralPost).degrees > FOV.degrees:
            return False

        for inputTarget in self.data:

            inputCoords = ICRSCoordinates(
                inputTarget['ra'], inputTarget['dec'])

            if (inputCoords - targetCoords).degrees < targetAvoid.degrees:
                return False

        return True

    def addTarget(self, target, bundleSize):

        newRow = []

        for column in self.data.colnames:
            if column.upper() == 'IFUDESIGNSIZE':
                newRow.append(bundleSize)
            elif column.upper() == 'PLATERA':
                newRow.append(self.data.meta['racen'])
            elif column.upper() == 'PLATEDEC':
                newRow.append(self.data.meta['deccen'])
            elif column.upper() in target.dtype.names:
                newRow.append(target[column.upper()])
            elif column.upper() == 'PRIORITY':
                newRow.append(np.max(self.data['priority'])+1)
            elif column.lower() in defaultValues:
                newRow.append(defaultValues[column.lower()])
            else:
                newRow.append(-999)

        self.data.add_row(newRow)

        log.important('autocomplete: added target with mangaid=' +
                      target['MANGAID']  + ' (ifudesignsize=' +
                      str(int(bundleSize)) + ')')

    def decollision(self, decollCatalogue, failOnCollision=False):

        self.data.sort('priority')

        centralPost = ICRSCoordinates(self.data.meta['racen'],
                                      self.data.meta['deccen'])

        targetsToRemove = []
        for ii in range(len(self.data)):

            inputTarget = self.data[ii]
            inputCoord = ICRSCoordinates(inputTarget['ra'],
                                         inputTarget['dec'])

            if (centralPost - inputCoord).degrees < centreAvoid.degrees:
                self.logCollision(
                    'mangaid={0} '.format(inputTarget['mangaid']) +
                    'rejected because collides with the central post')
                targetsToRemove.append(ii)

            if (centralPost - inputCoord).degrees > FOV.degrees:
                self.logCollision(
                    'mangaid={0} '.format(inputTarget['mangaid']) +
                    'rejected because it\'s outside the FOV')
                targetsToRemove.append(ii)

            for jj in range(ii+1, len(self.data)):
                otherTarget = self.data[jj]
                otherCoord = ICRSCoordinates(
                    otherTarget['ra'], otherTarget['dec'])
                if (otherCoord - inputCoord).degrees < targetAvoid.degrees:
                    self.logCollision(
                        'mangaid={0} '.format(inputTarget['mangaid']) +
                        'rejected because collides with target ' +
                        'mangaid={0}'.format(otherTarget['mangaid']))

        if len(targetsToRemove) > 0 and failOnCollision:
            raise GohanError('exiting because there was a target rejection '
                             'and failOnCollision=True')

        self.data.remove_rows(targetsToRemove)

        if decollCatalogue is False:
            return

        if isinstance(decollCatalogue, InputCatalogue):
            decollCatalogue = decollCatalogue.data.copy()

        if not isinstance(decollCatalogue, (table.Table, np.ndarray)):
            raise GohanError('decollision catalogue must be an astropy.Table'
                             ' or a numpy recarray.')

        if 'ra' not in decollCatalogue.dtype.names or \
                'dec' not in decollCatalogue.dtype.names:
            raise GohanError('decollision catalogue must have ra and dec '
                             'fields')

        inputCatCoords = np.array(
            [ICRSCoordinates(self['ra'][ii], self['dec'][ii])
             for ii in range(len(self))]
        )

        decollCatCoords = np.array(
            [ICRSCoordinates(decollCatalogue['ra'][ii],
                             decollCatalogue['dec'][ii])
             for ii in range(len(decollCatalogue))]
        )

        sepMatrix = separation_matrix(inputCatCoords, decollCatCoords)

        targetsToRemove = []
        for ii in range(len(self.data)):
            sepRow = sepMatrix[ii]
            for jj in range(len(sepRow)):
                if sepRow[jj].degrees < targetAvoid.degrees:
                    targetsToRemove.append(ii)
                    self.logCollision(
                        'mangaid={0} '.format(self.data['mangaid'][ii]) +
                        'rejected because collides with target with ' +
                        'RA={0:.5f}, Dec={0:.5f}'.format(
                            decollCatCoords[jj].ra.degrees,
                            decollCatCoords[jj].dec.degrees))
                    break

        self.data.remove_rows(targetsToRemove)

        if len(targetsToRemove) > 0 and failOnCollision:
            raise GohanError('exiting because there was a target rejection '
                             'and failOnCollision=True')

    def cropCatalogue(self):
        self.data.sort('priority')

        nRows = len(self.data)
        if nRows > self.nBundles:
            self.data = self.data[0:self.nBundles]
            log.debug('Cropped {0} targets with priority > {1}'.format(
                (nRows-len(self.data)), np.max(self['priority'])))

    def _createTableFromFile(self):

        ext = os.path.splitext(self.file)[1]

        if 'fits' in ext.lower():
            log.debug('File is a FITS file')
            data = fitsio.read(self.file)

        elif ext.lower() == 'par':
            log.info('File is a Yanny file')
            structureName = self.kwargs['structure'] if 'structure' in \
                self.kwargs else 'STRUCT1'

            data = yanny.yanny(self.file, np=True)[structureName]

        else:
            log.info('File format unknown. Trying astropy.Table ...')
            try:
                data = table.Table.read(self.file)
            except:
                raise GohanError('file format not understood. '
                                 'Use fits or par.')

        tileIDCol = self.conv('tileid')

        if self.tileid is not None:
            log.info('Only rows with tileid={0} have been selected'.format(
                self.tileid))
            if tileIDCol not in data.dtype.names:
                raise GohanError('No tileid column found in the data. Try '
                                 'using conversions to specify it.')
            data = data[np.where(data[tileIDCol] == self.tileid)]

        self.data = table.Table(data)
        self.checkColumns()

    def checkColumns(self):

        self.convertColumns()
        self.loweriseColumns()

        for col in defaultColumns:
            if col in self.colnames:
                continue
            elif col == 'priority':
                self.add_col(np.arange(1, len(self)+1), 'priority', int)
                log.debug('Automatically added column for priority.')
            elif col == 'sourcetype':
                self.add_col(defaultValues['sourcetype'], 'sourcetype', 'S10')
                log.debug('Automatically added column for sourcetype.')
            elif col == 'ifudesign':
                self.add_col(defaultValues['ifudesign'], 'ifudesign', int)
                log.debug('Automatically added column for ifudesign.')
            elif col == 'ifudesignsize':
                self.add_col(
                    defaultValues['ifudesignsize'], 'ifudesignsize', int)
                log.debug('Automatically added column for ifudesignsize.')
            elif col == 'manga_target1':
                self.add_col(
                    defaultValues['manga_target1'], 'manga_target1', int)
                log.debug('Automatically added column for manga_target1.')
            elif col == 'psfmag':
                self.add_col([defaultValues['psfmag']
                              for ii in range(len(self))], 'psfmag')
                log.debug('Automatically added psfmag for priority.')
            else:
                raise GohanError('mandatory column {0} not found.'.format(col))

    def checkData(self):

        self._fixDtypes()
        self.reorder()
        self.checkIFUCoords()

        if len(self) < self.nBundles:
            log.important('number of targets < number of IFUs '
                          '({0}<{1})'.format(len(self), self.nBundles))

    def checkIFUCoords(self):

        for coordType in ['ifu', 'target']:
            for col in ['ra', 'dec']:

                if not '{0}_{1}'.format(coordType, col) in self.data.colnames:
                    self.add_col([-999. for nn in range(len(self))],
                                 '{0}_{1}'.format(coordType, col), dtype=float)
                    log.debug('Adding {0}_{1} column'.format(coordType, col))

        for col in ['ra', 'dec']:
            for target in self:

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

    def convertColumns(self):
        for key, value in self.conversions.items():
            if value in self.colnames and key != value:
                self.data.rename_column(value, key)
                log.debug('Column {0} --> {1}'.format(value, key))

    def loweriseColumns(self):
        for colname in self.colnames:
            if colname != colname.lower():
                self.data.rename_column(colname, colname.lower())
                log.debug('Column {0} --> {1}'.format(
                    colname, colname.lower()))

    def reorder(self):
        newOrder = []
        for col in defaultColumns:
            if col in self.colnames:
                newOrder.append(col)
        for col in self.colnames:
            if col not in newOrder:
                newOrder.append(col)
        self.data = self[newOrder]
        log.debug('Column order rearranged')

    def _fixDtypes(self):
        formats = defaultMaNGAInput['MANGAINPUT'].dtype
        for name in formats.names:
            dd = formats[name]
            if len(dd.shape) > 0:
                dd = dd.subdtype[0]
            if self[name].dtype != dd:
                oldDtype = self[name].dtype
                tmpData = self.data[name]
                self.data.remove_column(name)
                self.data.add_column(table.Column(data=tmpData.astype(dd),
                                                  name=name))
                log.debug('Column {0}: dtype {1} --> {2}'.format(
                    name, oldDtype, dd))

    def add_col(self, data, name, dtype=None):
        if not isinstance(data, (list, tuple, np.ndarray)):
            data = len(self) * [data]
        self.data.add_column(table.Column(data=data, name=name, dtype=dtype))

    def _createMeta(self, meta):

        defaultPairs = defaultMaNGAInput.pairs()
        defaultMeta = OrderedDict([(pair, defaultMaNGAInput[pair])
                                  for pair in defaultPairs])

        self.data.meta = defaultMeta
        self.data.meta.update(meta)

        log.debug('Saved pairs to InputCatalogue metadata.')

        for pair in self.meta:

            if pair.lower() == 'racen':
                platera = self.conv('platera')
                if platera in self.colnames:
                    self.data.meta['racen'] = self.data[0][platera]
                    log.debug('Pair racen set from {0}'.format(platera))

            elif pair.lower() == 'deccen':
                platedec = self.conv('platedec')
                if platedec in self.colnames:
                    self.data.meta['deccen'] = self.data[0][platedec]
                    log.debug('Pair deccen set from {0}'.format(platedec))

            elif pair.lower() == 'targettype':
                self.data.meta['targettype'] = targettypeDic[self.type]
                log.debug('Pair targettype set by using type keyword.')

            elif pair.lower() == 'locationid':
                if self.tileid is not None:
                    self.data.meta['locationid'] = self.tileid

            elif pair.lower() == 'tileid':
                if self.tileid is not None:
                    self.data.meta['tileid'] = self.tileid

            if pair.lower() == 'designid':
                continue

            if self.meta[pair] != 'None':
                continue

            raise GohanError('mandatory pair {0} not set.'.format(
                pair.lower()))

        if self.data.meta['fieldname'] == 'NA':
            self.data.meta['fieldname'] = \
                'MJ{0:.5f}{1:+.5f}'.format(
                    self.data.meta['racen'], self.data.meta['deccen'])

    def conv(self, key):
        return self.conversions[key] if key in self.conversions else key

    def write(self, *args, **kwargs):
        if os.path.exists(args[0]):
            os.remove(args[0])
        self.data.write(*args, **kwargs)
        log.debug('InputCatalogue saved to {0}.'.format(args[0]))

    @classmethod
    def read(cls, path, format=None, decollision=None, **kwargs):

        if not os.path.exists(path):
            raise GohanError('The file to read does not exists.')

        data = table.Table.read(path, format=format)
        tileid = data.meta['tileid']
        type = [key for key in targettypeDic
                if targettypeDic[key] == data.meta['targettype']][0]

        return cls(tileid=tileid, format='file', type=type,
                   decollision=decollision, **kwargs)

    def logCollision(self, text):

        if self.type == 'SCI':
            warnings.warn(text, GohanCollisionWarning)
        else:
            log.info(text)

    def plotStaralt(self, date, filename=None):
        """Plots a staralt graph for the catalogue centre position.

        Parameters
        ----------
        date : ``astropy.time.Time`` instance, str
            The date for which the plot is created. It can be either a Time
            instance or a string containing an ISO time and date that can be
            understood by ``astropy.time.Time``.
        filename : str, optional
            The filename of the output plot. The default is Staralt-XXXX.pdf
            where XXXX is the locationID.

        """

        if Staralt is None:
            raise ImportError(
                'No Staralt functionality. Check your instalation.')

        if filename is None:
            filename = 'Staralt-{0:04d}.pdf'.format(self.locationID)

        ra = self.meta['racen']
        dec = self.meta['deccen']

        if self.tileid is not None:
            st = Staralt(ra, dec, date, name='tileID={0:04d}'.format(
                self.tileid))
        else:
            st = Staralt(ra, dec, date)

        st.save(filename)

    def copy(self):
        return self.data.copy()

    def __repr__(self):
        return self.data.__repr__()

    def __str__(self):
        return self.data.__str__()

    def __getitem__(self, slice):
        return self.data[slice]

    def __len__(self):
        return len(self.data)

    @property
    def colnames(self):
        return self.data.colnames

    @property
    def meta(self):
        return self.data.meta
