#!/usr/bin/env python
# encoding: utf-8
#
# generate_designs.py
#
# Created by José Sánchez-Gallego on 22 Apr 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import warnings

import astropy.table as table

from Gohan import log, config, readPath
from Gohan.core.colourPrint import _color_text
from Gohan.exceptions import GohanUserWarning
from Gohan.utils import utils
from Gohan.utils.selectSkies import selectSkies

from Gohan.PlateInput import PlateInput


def do_one_apogee2_manga(plateRun, platerun_dir, field, reject_science=[], reject_standard=[],
                         niter=1, exclude=[]):
    """Does one apogee2-manga design and returns the list of allocated mangaids."""

    reject_science = list(reject_science) + exclude
    reject_standard = list(reject_standard) + exclude

    designID = int(field['DesignID'])
    fieldName = field['FieldName']

    log.important(_color_text('designID: {0}, fieldName: {1} (iteration {2})'
                              .format(designID, fieldName, niter), 'lightred'))

    raCen = float(field['RA'])
    decCen = float(field['Dec'])

    sciPath = os.path.join(platerun_dir, 'inputs/{0}_target.fits'.format(fieldName))

    sciCat = table.Table.read(sciPath)
    stdCat = table.Table.read(os.path.join(platerun_dir,
                                           'inputs/{0}_standards.fits'.format(fieldName)))

    for column in ['FIELD', 'PS_APASS']:
        if column in sciCat.colnames:
            sciCat.remove_column(column)

    skyRaw = os.path.join(readPath(config['skiesPath']), 'sky_{0}.fits'.format(fieldName))
    assert(os.path.exists(skyRaw))

    mangaScience = PlateInput(designID, 'science', catalogues=sciCat,
                              surveyMode='apogeeLead', plateRun=plateRun,
                              silentOnCollision=True, sort=False,
                              raCen=raCen, decCen=decCen,
                              rejectTargets=reject_science)
    mangaScience.write(toRepo=False)

    # Checks how many special targets survive
    try:
        mangaScience_path = mangaScience.getDefaultFilename()
        special = utils.calculate_special_targets(mangaScience_path, sciPath)
        if special is None or special[1] == 0:
            log.important('no special targets found for this design.')
        else:
            survived, total = special
            log.important('{0} out of {1} special targets survived '
                          'decollision.'.format(survived, total))
    except Exception as ee:
        warnings.warn('failed determining how many special targets '
                      'survived for {0}: {1}.'.format(mangaScience_path, ee), GohanUserWarning)

    mangaStandard = PlateInput(designID, 'standard',
                               catalogues=stdCat,
                               surveyMode='apogeeLead', plateRun=plateRun,
                               silentOnCollision=False, sort=False,
                               decollidePlateInputs=[mangaScience],
                               raCen=raCen, decCen=decCen,
                               rejectTargets=reject_standard)
    mangaStandard.write(toRepo=False)

    result = selectSkies(skyRaw, designID, fieldName, raCen, decCen, raise_error=False)

    if isinstance(result, table.Table):
        mangaSky = PlateInput(designID, 'sky', catalogues=result,
                              surveyMode='apogeeLead', plateRun=plateRun,
                              silentOnCollision=True, sort=False,
                              raCen=raCen, decCen=decCen)
        mangaSky.write(toRepo=False)
    else:
        assert result[0] is False, 'cannot parse the result of selectSkies'
        warnings.warn(result[1], GohanUserWarning)
        mangaID, ifuDesign = result[2:]

        if ifuDesign < 1000:
            reject_standard.append(mangaID)
            log.debug('added {0} to the reject_standard list.'.format(mangaID))
        else:
            reject_science.append(mangaID)
            log.debug('added {0} to the reject_science list.'.format(mangaID))

        return do_one_apogee2_manga(plateRun, platerun_dir, field,
                                    reject_science=reject_science,
                                    reject_standard=reject_standard,
                                    niter=niter + 1)

    return mangaScience.getMangaIDs()


def update_platedefinition_apogee2manga(field, plateRun, platerun_dir):
    """Creates an updated copy of the plateDefinition file committed by APOGEE."""

    designID = int(field['DesignID'])
    fieldName = field['FieldName']
    plateDefinitionPath = utils.getPlateDefinition(designID)

    if not os.path.exists(plateDefinitionPath):
        raise RuntimeError('{0} not found'.format(plateDefinitionPath))

    plateDefinition = open(plateDefinitionPath, 'r').read().splitlines()

    updatedPlateDefinition = []
    for line in plateDefinition:
        if line.strip() == '':
            updatedPlateDefinition.append(line)
        elif line.strip()[0] == '#':
            updatedPlateDefinition.append(line)
        elif 'platedesignVersion' in line:
            updatedPlateDefinition.append('platedesignVersion v1')
        elif 'plateType' in line:
            updatedPlateDefinition.append('plateType APOGEE2-MANGA')
        elif 'priority' in line:
            updatedPlateDefinition.append('priority  1 2 3 4 5 6')
        elif 'plateInput3' in line:
            updatedPlateDefinition.append(
                line.replace('plateInput3', 'plateInput6'))
        elif 'nInputs' in line:
            updatedPlateDefinition.append('nInputs 6')
        elif 'plateInput2' in line:
            updatedPlateDefinition.append(line)
            updatedPlateDefinition.append(
                'plateInput3 manga/{0}/mangaScience_{1}_{2}.par'
                .format(plateRun, fieldName, designID))
            updatedPlateDefinition.append(
                'plateInput4 manga/{0}/mangaStandard_{1}_{2}.par'
                .format(plateRun, fieldName, designID))
            updatedPlateDefinition.append(
                'plateInput5 manga/{0}/mangaSky_{1}_{2}.par'
                .format(plateRun, fieldName, designID))
        else:
            updatedPlateDefinition.append(line)

    updatedPlateDefinition += ['',
                               'defaultSurveyMode apogeeLead',
                               'plateLead APOGEE',
                               '',
                               '# Skies for MaNGA',
                               'respectIFUID yes',
                               'skyType NONE',
                               'plateDesignSkies NONE',
                               '']

    updatedPlateDefinitionString = '\n'.join(updatedPlateDefinition)

    plateDefinition_path_new = os.path.join(platerun_dir, os.path.basename(plateDefinitionPath))
    unitPlateDefinition = open(plateDefinition_path_new, 'w')
    unitPlateDefinition.write(updatedPlateDefinitionString)
    unitPlateDefinition.close()

    log.important('saved updated plateDefinition as {0}'.format(plateDefinition_path_new))


def designs_apogee2manga(platerun, platerun_dir, plate_data_path,
                         special=False, exclude=[], **kwargs):
    """Designs an APOGEE2-MaNGA platerun."""

    fields = table.Table.read(plate_data_path, format='ascii.commented_header')

    if special is True:
        log.important('not designing field. Only printing summary of special targets.')
        utils.print_special_summary(plate_data_path)
        return False

    mastar_targets = utils.getStellarLibraryTargets()
    mastar_mangaid = set(mastar_targets['mangaid'].tolist())

    allocated_run = set([])

    reject_science = mastar_mangaid.union(allocated_run)

    for field in fields:

        allocated_one = do_one_apogee2_manga(platerun, platerun_dir, field,
                                             reject_science=list(reject_science),
                                             exclude=exclude)

        update_platedefinition_apogee2manga(field, platerun, platerun_dir)

        reject_science = reject_science.union(set(allocated_one))

    utils.print_special_summary(plate_data_path)

    return True


def generate_designs(platerun, plate_data, **kwargs):
    """Generates plate inputs, plate definitons, and platePlans lines for a list of plates.

    Parameters:
        platerun (str):
            A string with the format ``'{year}.{month}.{id}.{plateType}'`` where ``{id}`` is
            a unique identifier for a given ``{year}.{month}`` (usually ``'x'``), and
            ``{plateType}`` is one of ``'manga-apogee2'``, ``'apogee2-manga'``, or ``'manga'``.
        plate_data (str):
            The path to the file containing the data for the plates to design. For MaNGA-led
            runs this will usually be a list of fields selected using Totoro. For APOGEE-led
            runs it will be a table with format ascii.commented_header and, at least, the
            following columns: ``FieldName, DesignID, RA, Dec, HA, Epoch``.

    """

    platerun = platerun.strip().lower()
    plate_type = platerun.split('.')[-1]
    assert plate_type in ['manga-apogee2', 'apogee2-manga', 'manga'], 'invalid plateType'

    platerun_dir = os.path.join(os.environ['MANGAWORK_DIR'], 'manga/platedesign/plateruns',
                                platerun)
    assert os.path.exists(platerun_dir), 'cannot find platerun dir {0}'.format(platerun_dir)

    assert os.path.exists(plate_data), 'cannot find path {0}'.format(plate_data)

    log.info('Plate Run: {0}'.format(platerun))
    log.info('Plate Type: {0}'.format(plate_type))
    log.info('Plate Data file: {0}'.format(plate_data))
    log.info('Working on platerun directory: {0}'.format(platerun_dir))

    if plate_type == 'apogee2-manga':
        keep_log = designs_apogee2manga(platerun, platerun_dir, plate_data, **kwargs)

    if keep_log:
        log_path = os.path.join(platerun_dir, platerun + '.log')
        log.info('copying log to {0}'.format(log_path))
        log.saveLog(log_path)

    return
