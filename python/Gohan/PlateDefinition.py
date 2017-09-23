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
from astropy.coordinates import SkyCoord

from Gohan import log, config, readPath
from Gohan.exceptions import GohanError
# from Gohan.utils import getPlateTemperature


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

    def _addPlateInput(self, lines, plateInput, path):
        """Adds a plateInput to the lines of a new plateDefinition."""

        # fullPath = os.path.join(
        #     readPath(config['platelist']), 'inputs', path)

        # if not os.path.exists(fullPath):
        #     raise GohanError('input {0} does not exist'.format(path))

        # If plateInput does not include the id, uses the following to the
        # existing ones.
        if plateInput == 'plateInput':
            maxPlateInput = 0
            for line in lines:
                if line.strip() == '':
                    continue
                elif line.strip()[0] == '#':
                    continue
                elif 'plateInput' in line.strip():
                    pp = int(line.strip().split()[0][10:])
                    if pp > maxPlateInput:
                        maxPlateInput = pp
            plateInput = 'plateInput' + str(maxPlateInput + 1)

        plateInputID = int(plateInput[10:])

        lastLine = None
        for ii, line in enumerate(lines):
            if 'plateInput' in line.strip():
                lastLine = ii

        if lastLine is None:
            if lines[-1].strip() != '':
                lines.append('')
            lines.append(plateInput + ' ' + str(path))
        else:
            lines.insert(lastLine + 1, plateInput + ' ' + str(path))

        for ii, line in enumerate(lines):
            if line.strip() == '':
                continue
            elif line.strip()[0] == '#':
                continue
            elif 'nInputs' in line:
                nInputs = int(line.strip().split()[1])
                lines[ii] = 'nInputs ' + str(nInputs + 1)
            elif 'priority' in line:
                lines[ii] = lines[ii].strip() + ' ' + str(plateInputID)

    def createPlateDefinition(self):
        """Return the plateDefinition yanny object."""

        if self.surveyMode == 'mangaOnly':
            plateDefinitionTemplate = readPath(
                '+etc/mangaDefinition_Default.par')
        elif self.surveyMode == 'mangaLead':
            plateDefinitionTemplate = readPath(
                '+etc/mangaDefinition_coDesigned_Default.par')
        else:
            raise GohanError(
                'surveyMode {0} is invalid'.format(self.surveyMode))

        raCen = self.plateInputs[0].raCen
        decCen = self.plateInputs[0].decCen

        coords = SkyCoord(raCen, decCen, unit='deg')
        apogeeCoords = '{0:03d}{1:+03d}'.format(
            int(np.round(coords.galactic.galactic.l.deg)),
            int(np.round(coords.galactic.galactic.b.deg)))

        nInputs = len(self.plateInputs)
        priority = ' '.join([str(ii + 1) for ii in range(nInputs)])

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
            inputs += [['plateInput{0:d}'.format(nn + 1), inputFilename]]

        inputs = OrderedDict(inputs)
        defDict.update(inputs)

        template = open(plateDefinitionTemplate, 'r').read().splitlines()

        for ii, line in enumerate(template):
            if line.strip() == '':
                continue
            elif line.strip()[0] == '#':
                continue

            delete_keys = []

            for key in defDict:
                if key in line:
                    template[ii] = key + ' ' + str(defDict[key])
                    delete_keys.append(key)

            for key in delete_keys:
                if key in defDict:
                    del defDict[key]

        if len(defDict) > 0:
            if template[-1].strip() != '':
                template.append('')
            for key in defDict:
                if 'plateInput' in key:
                    self._addPlateInput(template, key, defDict[key])
                else:
                    template.append(key + ' ' + str(defDict[key]))

        # If this is a MaNGA-APOGEE design, adds the APOGEE inputs
        if self.surveyMode == 'mangaLead':
            for inputType in ['STA', 'SCI', 'SKY']:
                path = 'apogee/{0}/plateInput_{1}_MGA_{2}_{3}.par'.format(
                    self.plateRun, apogeeCoords, inputType, self.designid)
                self._addPlateInput(template, 'plateInput', path)

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
            definitionYanny = self.createPlateDefinition()
            # definitionYanny.set_filename(filename)

            unit = open(filename, 'w')
            for line in definitionYanny:
                unit.write(line + '\n')
            unit.close()

        log.info('plateDefinition file {0} saved.'.format(filename))

        if toRepo:

            if not os.path.exists(plateDefinitionPath):
                os.makedirs(plateDefinitionPath)

            sh.copy(filename, plateDefinitionRepoFile)
            log.info('plateDefinition file {0} '.format(filename) +
                     'saved to platelist repo.')

        return filename

    def getPlatePlans(self, epoch, temp=None, name=None, comment=' ',
                      plateid=None, plateRun=None, **kwargs):
        """Returns the platePlans text. `epoch` muct be the epoch of the
        observation as a float or a date in the format `'YYYY-MM-DD'`"""

        plansDic = {}

        if isinstance(epoch, Real):
            pass
        else:
            epoch = time.Time(epoch.strip() + ' 00:00:00', format='iso').byear
            plansDic['epoch'] = str(epoch)
        plansDic['epoch'] = str(epoch)

        plansDic['raCen'] = self.plateInputs[0].raCen
        plansDic['decCen'] = self.plateInputs[0].decCen

        plansDic['plateID'] = '@plateid' if plateid is None else plateid
        plansDic['comment'] = comment

        if temp is None:
            # Calculates year fraction from the epoch.
            # yearFraction = epoch - int(epoch)
            plansDic['temp'] = 5  # getPlateTemperature(yearFraction)
        else:
            plansDic['temp'] = temp

        if name is not None:
            plansDic['name'] = name
        else:
            plansDic['name'] = self.plateInputs[0].fieldName

        plansDic['locationID'] = self.locationid
        plansDic['designID'] = self.designid
        plansDic['plateRun'] = '@platerun' if plateRun is None else plateRun
        plansDic['chunk'] = '@platerun' if plateRun is None else plateRun

        plansDic.update(config['plateTypes'][self.surveyMode])

        platePlans = platePlansTemplate.format(**plansDic)
        platePlans = platePlans.replace('\n', ' ')

        return platePlans
