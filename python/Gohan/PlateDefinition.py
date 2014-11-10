#!/usr/bin/env python
# encoding: utf-8
"""
PlateDefinition.py

Created by José Sánchez-Gallego on 10 Nov 2014.
Licensed under a 3-clause BSD license.

Revision history:
    10 Nov 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function

import numpy as np
from collections import OrderedDict
import os
import shutil as sh
from numbers import Real

from astropy import time

from Gohan import log, config, readPath
from Gohan.utils import yanny


platePlansTemplate = """PLATEPLANS {plateID} {designID} {locationID} -1
{platedesignversion} {{ 0.0 0.00000 0.00000 0.00000 0.00000 0.00000 }}
{temp} {epoch} {raCen} {decCen} {survey} {programname} {drillstyle} \" \"
{plateRun} {chunk} \"{name}\" \"{comment}\"
"""


class PlateDefinition(object):
    """A plateDefinition-like class."""

    def __init__(self, plateInputs, designid=None, surveyMode=None,
                 plateRun=None, locationid=None, **kwargs):

        assert isinstance(plateInputs, (list, tuple, np.ndarray))
        assert len(plateInputs) > 0

        self.plateInputs = plateInputs

        self.designid = designid if designid is not None \
            else plateInputs[0].designid
        self.surveyMode = surveyMode if surveyMode is not None \
            else plateInputs[0].surveyMode
        self.plateRun = plateRun if plateRun is not None \
            else plateInputs[0].plateRun
        self.locationid = locationid if locationid is not None \
            else plateInputs[0].locationid

    def getPlateDefinitionYanny(self):
        """Return the plateDefinition yanny object."""

        plateDefinitionTemplate = readPath('+etc/mangaDefinition_Default.par')

        raCen = self.plateInputs[0].raCen
        decCen = self.plateInputs[0].decCen
        nInputs = len(self.plateInputs)
        priority = ' '.join([str(ii+1) for ii in range(nInputs)])

        defDict = OrderedDict(
            [['raCen', raCen], ['decCen', decCen], ['nInputs', nInputs],
             ['priority', priority], ['designID', self.designid]])

        defDict['plateType'] = \
            config['plateTypes'][self.surveyMode]['plateType']
        defDict['plateLead'] = \
            config['plateTypes'][self.surveyMode]['plateLead']
        defDict['platedesignversion'] = \
            config['plateTypes'][self.surveyMode]['platedesignversion']
        defDict['defaultSurveyMode'] = \
            config['plateTypes'][self.surveyMode]['defaultSurveyMode']

        inputs = []
        for nn, plateInput in enumerate(self.plateInputs):
            inputFilename = 'manga/{0}/{1}'.format(
                self.plateRun, self.plateInputs[nn].getDefaultFilename())
            inputs += [['plateInput{0:d}'.format(nn+1), inputFilename]]

        inputs = OrderedDict(inputs)
        defDict.update(inputs)

        template = yanny.yanny(plateDefinitionTemplate)

        for key in defDict:
            template[key] = defDict[key]

        return template

    def write(self, toRepo=False):
        """Writes the plateDefinition file."""

        filename = 'plateDefinition-{0:06d}.par'.format(self.designid)

        pathRoot = '{0:04d}XX/'.format(int(self.designid / 100))
        plateDefinitionPath = os.path.join(
            readPath(config['platelist']),
            'definitions', pathRoot)
        plateDefinitionRepoFile = os.path.join(plateDefinitionPath, filename)

        if self.surveyMode == 'apogeeLead':
            plateDefinition = self._getPlateDefinitionAppend(
                plateDefinitionRepoFile)

            blob = open(filename, 'w')
            for line in plateDefinition:
                blob.write(line + '\n')
            blob.close()

        else:
            definitionYanny = self.getPlateDefinitionYanny()
            definitionYanny.set_filename(filename)

            if os.path.exists(filename):
                os.remove(filename)

            definitionYanny.write()

        log.info('plateDefinition file {0} saved.'.format(filename))

        if toRepo:

            if not os.path.exists(plateDefinitionPath):
                os.makedirs(plateDefinitionPath)

            sh.copy(filename, plateDefinitionRepoFile)
            log.info('plateDefinition file {0} '.format(filename) +
                     'saved to platelist repo.')

        return filename

    def getPlatePlans(self, epoch, temp=5.0, name=None, comment=' ',
                      plateid=1000, **kwargs):
        """Returns the platePlans text. `epoch` muct be the epoch of the
        observation as a float or a date in the format `'YYYY-MM-DD'`"""

        plansDic = {}

        plansDic['raCen'] = self.plateInputs[0].raCen
        plansDic['decCen'] = self.plateInputs[0].decCen

        plansDic['plateID'] = plateid
        plansDic['comment'] = comment
        plansDic['temp'] = str(temp)
        if name is not None:
            plansDic['name'] = name
        else:
            plansDic['name'] = self.plateInputs[0].fieldName

        plansDic['locationID'] = self.locationid
        plansDic['designID'] = self.designid
        plansDic['plateRun'] = self.plateRun
        plansDic['chunk'] = self.plateRun

        if isinstance(epoch, Real):
            pass
        else:
            epoch = time.Time(epoch.strip() + ' 00:00:00', format='iso').byear
            plansDic['epoch'] = str(epoch)
        plansDic['epoch'] = str(epoch)

        plansDic.update(config['plateTypes'][self.surveyMode])

        platePlans = platePlansTemplate.format(**plansDic)
        platePlans = platePlans.replace('\n', ' ')

        return platePlans
