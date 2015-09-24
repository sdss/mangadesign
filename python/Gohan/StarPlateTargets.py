#!/usr/bin/env python
# encoding: utf-8
"""
StarPlateTargets.py

Created by José Sánchez-Gallego on 26 Mar 2015.
Licensed under a 3-clause BSD license.

Revision history:
    26 Mar 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

from sdss.utilities import yanny
from astropy import table
import numpy as np

from Gohan.exceptions import GohanPlateTargetsError, GohanPlateTargetsWarning
from Gohan.PlateTargets import PlateTargets, _toLowerCase
from Gohan.utils import utils
from Gohan import readPath, log

import warnings
import os


class StarPlateTargets(PlateTargets):

    def __init__(self, arg=None, **kwargs):
        """A class to handle stellar library plateTargets files.

        Parameters
        ----------
        arg : string or None.
            The path to the plateTargets file to use. If None, the template
            will be used.

        """

        self.template = False
        self._nAppended = 0
        self.catalogid = 'MaSTAR'

        if arg is not None:
            self.path = arg
            if not os.path.exists(self.path):
                raise GohanPlateTargetsError('path {0} cannot be found'
                                             .format(self.path))

        else:
            starPlateTargetsPath = os.path.join(
                os.path.dirname(utils.getPlateTargetsPath(1)),
                'starPlateTargets.par')

            if os.path.exists(starPlateTargetsPath):
                self.path = starPlateTargetsPath
            else:
                self.path = readPath('+plateTargets/starPlateTargets.template')
                if self.path is None:
                    raise GohanPlateTargetsError(
                        'neither the plateTargets nor the template '
                        'for catalogid={0} can be found'.format(catalogid))
                warnings.warn('using template for starPlateTargets',
                              GohanPlateTargetsWarning)
                self.template = True

        data = yanny.yanny(self.path, np=True)

        self.comments = self._getComments(data)
        self.structure = table.Table(data['PLTTRGT'])

        if self.template:
            self.structure.remove_row(0)

    def write(self):
        """Writes the current instance to a Yanny file."""

        path = os.path.join(os.path.dirname(utils.getPlateTargetsPath(1)),
                            'starPlateTargets.par')

        path, nAppended = super(StarPlateTargets, self).write(path=path,
                                                              useCatID=False)

        log.debug('starPlateTargets saved to {0}'.format(path))
        log.debug('{0} targets appended to starPlateTargets'.format(nAppended))

        return path, nAppended

    def addTargets(self, **kwargs):
        """Adds targets to the stellar library plate targets."""

        addedIndices = []

        overwrite = kwargs.get('overwrite', False)

        commonData, mangaids, plateid = super(StarPlateTargets,
                                              self).addTargets(**kwargs)

        specificData = self.getSpecificData(mangaids, plateid)

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
                        warnings.warn('replacing target mangaid={0} '
                                      'in plateid={1}'
                                      .format(mangaid, plateid),
                                      GohanPlateTargetsWarning)
                    else:
                        # If overwrite is False, skips this target.
                        log.debug('skipping mangaid={0} because it is already '
                                  'in plateTargets.'.format(mangaid))
                        continue

            # Cleans up values
            targetData = self._cleanupTargetData(targetData)

            # Applies target fixes
            targetData = self._applyTargetFix(targetData)

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

            log.debug('mangaid={0} added to starPlateTargets'.format(mangaid))

        return self.structure[addedIndices]

    def getSpecificData(self, mangaids, plateid):
        """Gathers star parameters from mangaScience."""

        if plateid is None:
            raise GohanPlateTargetsError('plateid required to retrieve '
                                         'specific data')

        requiredFields = utils.getRequiredPlateTargetsColumns()
        specificFields = [field for field in self.structure.colnames
                          if field not in requiredFields]

        mangaSciencePath = utils.getMangaSciencePath(plateid, format='plateid')
        mangaScience = table.Table(
            yanny.yanny(mangaSciencePath, np=True)['MANGAINPUT'])
        mangaScience = _toLowerCase(mangaScience)
        mangaScience['mangaid'] = map(lambda xx: xx.strip(),
                                      mangaScience['mangaid'])

        specificData = {}
        for mangaid in mangaids:

            specificData[mangaid] = {}

            row = mangaScience[mangaScience['mangaid'] == mangaid.strip()]

            for field in specificFields:
                if field in mangaScience.colnames:
                    specificData[mangaid][field] = row[field]
                else:
                    if field != 'extinction':
                        specificData[mangaid][field] = -999.
                    else:
                        specificData[mangaid][field] = 0.

        return specificData
