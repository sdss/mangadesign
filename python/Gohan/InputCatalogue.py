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

    """

    def __init__(self, tileid=None, format='file', type='SCI', file=None,
                 conversions=None, fill=False, meta=None,
                 removeSuperfluous=False, verbose=True,
                 decollision=False, failOnCollision=False, **kwargs):

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

        self.cropCatalogue()
        self.checkData()

    def removeSuperfluous(self):
        removedCols = []
        for colname in self.colnames:
            if colname not in defaultColumns:
                self.data.remove_column(colname)
                removedCols.append(colname)
        log.debug('Removed superfluous columns {0}'.format(removedCols))

    def decollision(self, decollCatalogue, failOnCollision=False):

        self.data.sort('priority')

        centralPost = ICRSCoordinates(self.data.meta['racen'],
                                      self.data.meta['deccen'])

        FOV = AngularCoordinate(config['decollision']['FOV'])
        centreAvoid = AngularCoordinate(config['decollision']['centreAvoid'])
        targetAvoid = AngularCoordinate(config['decollision']['targetAvoid'])

        targetsToRemove = []
        for ii in range(len(self.data)):

            inputTarget = self.data[ii]
            inputCoord = ICRSCoordinates(inputTarget['ra'],
                                         inputTarget['dec'])

            if (centralPost - inputCoord).degrees < centreAvoid.degrees:
                warnings.warn(
                    'mangaid={0} '.format(inputTarget['mangaid']) +
                    'rejected because collides with the central post',
                    GohanCollisionWarning)
                targetsToRemove.append(ii)

            if (centralPost - inputCoord).degrees > FOV.degrees:
                warnings.warn(
                    'mangaid={0} '.format(inputTarget['mangaid']) +
                    'rejected because it\'s outside the FOV',
                    GohanCollisionWarning)
                targetsToRemove.append(ii)

            for jj in range(ii+1, len(self.data)):
                otherTarget = self.data[jj]
                otherCoord = ICRSCoordinates(
                    otherTarget['ra'], otherTarget['dec'])
                if (otherCoord - inputCoord).degrees < targetAvoid.degrees:
                    warnings.warn(
                        'mangaid={0} '.format(inputTarget['mangaid']) +
                        'rejected because collides with target ' +
                        'mangaid={0}'.format(otherTarget['mangaid']),
                        GohanCollisionWarning)

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
                    warnings.warn(
                        'mangaid={0} '.format(self.data['mangaid'][ii]) +
                        'rejected because collisions with target with ' +
                        'RA={0:.5f}, Dec={0:.5f}'.format(
                            decollCatCoords[jj].ra.degrees,
                            decollCatCoords[jj].dec.degrees),
                        GohanCollisionWarning)
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
                self.add_col(self.type, 'sourcetype', 'S10')
                log.debug('Automatically added column for sourcetype.')
            elif col == 'ifudesign':
                self.add_col(-666, 'ifudesign', int)
                log.debug('Automatically added column for ifudesign.')
            elif col == 'ifudesignsize':
                self.add_col(-666, 'ifudesignsize', int)
                log.debug('Automatically added column for ifudesignsize.')
            elif col == 'manga_target1':
                self.add_col(0.0, 'manga_target1', int)
                log.debug('Automatically added column for manga_target1.')
            elif col == 'psfmag':
                self.add_col([[0.0, 0.0, 0.0, 0.0, 0.0]
                              for ii in range(len(self))], 'psfmag')
                log.debug('Automatically added psfmag for priority.')
            else:
                raise GohanError('mandatory column {0} not found.'.format(col))

    def checkData(self):

        self._fixDtypes()
        self.reorder()

        if len(self) < self.nBundles:
            warnings.warn('number of targets < number of IFUs '
                          '({0}<{1})'.format(len(self), self.nBundles))

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
            if self.tileid is not None:
                self.data.meta['fieldname'] = self.tileid
            elif self.data.meta['tileid'] != 'None':
                self.data.meta['fieldname'] = self.data.meta['tileid']

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
