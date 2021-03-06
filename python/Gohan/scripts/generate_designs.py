#!/usr/bin/env python
# encoding: utf-8
#
# generate_designs.py
#
# Created by José Sánchez-Gallego on 22 Apr 2017.


from __future__ import absolute_import, division, print_function

import os
import pathlib
import re

import astropy.table as table
import numpy
from pydl.pydlutils import yanny

from Gohan import config, log, readPath
from Gohan.core.color_print import _color_text
from Gohan.exceptions import GohanUserWarning
from Gohan.PlateDefinition import PlateDefinition
from Gohan.PlateInput import PlateInput
from Gohan.utils import utils
from Gohan.utils.selectSkies import selectSkies

from .update_addendas import update_addendas


def read_catalogue(fn):

    sciCat = table.Table.read(str(fn))

    for column in ['FIELD', 'PS_APASS']:
        if column in sciCat.colnames:
            sciCat.remove_column(column)

    return sciCat


def do_one_apogee2_manga(plateRun, platerun_dir, field, reject_science=[],
                         reject_standard=[], niter=1, exclude=[]):
    """Does one apogee2-manga design and returns the list of allocated mangaids."""

    reject_science = list(reject_science) + exclude
    reject_standard = list(reject_standard) + exclude

    designID = int(field['DesignID'])
    fieldName = field['FieldName']

    log.important(_color_text('designID: {0}, fieldName: {1} (iteration {2})'
                              .format(designID, fieldName, niter), 'lightred'))

    raCen = float(field['RA'])
    decCen = float(field['Dec'])
    tileRad = float(field['TileRad']) \
        if 'TileRad' in field.colnames else config['decollision']['FOV']

    science_target = pathlib.Path(platerun_dir) / f'inputs/{fieldName}_target.fits'
    science_target_high = pathlib.Path(platerun_dir) / f'inputs/{fieldName}_target_high.fits'
    science_target_low = pathlib.Path(platerun_dir) / f'inputs/{fieldName}_target_low.fits'

    mangaScience_high = None
    exclude_science_ifudesigns = []
    n_targets_high = 0
    if science_target_high.exists():

        # assert science_target_low.exists(), \
        #     'found {!r} file without a _low version.'.format(str(science_target_high.name))

        assert not science_target.exists(), \
            'found incompatible {!r} and {!r}'.format(str(science_target),
                                                      str(science_target_low))

        log.warning(f'found high priority target catalogue {str(science_target_high.name)!r}',
                    GohanUserWarning)

        science_cat_high = read_catalogue(science_target_high)

        # Runs high priority targets. We allocate all high priority targets so no
        # rejects list is passed. But if a high priority target gets rejected
        # because there are not enough skies around it, we carry that reject.

        mangaScience_high = PlateInput(designID, 'science', catalogues=science_cat_high,
                                       surveyMode=None, plateRun=plateRun,
                                       decollideExternal=False,
                                       silentOnCollision=True, sort=False,
                                       fieldName=fieldName,
                                       raCen=raCen, decCen=decCen,
                                       rejectTargets=reject_science,
                                       FOV=tileRad)
        n_targets_high = len(mangaScience_high.getMangaIDs())

        if hasattr(mangaScience_high, 'mangaInput') and n_targets_high > 0:
            science_target = science_target_low
            exclude_science_ifudesigns = mangaScience_high.mangaInput['ifudesign']
        else:
            log.warning('no high priority targets survive design!!!', GohanUserWarning)
            mangaScience_high = None

    if n_targets_high == 17:

        mangaScience = mangaScience_high
        mangaScience_low = None

    elif science_target.exists() or science_target_low.exists():

        if not science_target.exists():
            science_target = science_target_low

        science_cat = read_catalogue(science_target)

        mangaScience = PlateInput(designID, 'science', catalogues=science_cat,
                                  surveyMode='apogeeLead', plateRun=plateRun,
                                  silentOnCollision=True, sort=False,
                                  raCen=raCen, decCen=decCen,
                                  rejectTargets=reject_science,
                                  exclude_ifudesigns=exclude_science_ifudesigns,
                                  FOV=tileRad)

        mangaScience_low = mangaScience.copy()

        if mangaScience_high is not None:
            log.info('merging high priority targets')
            mangaScience.merge(mangaScience_high)
            assert len(mangaScience.mangaInput) == 17

    elif not science_target.exists() and not science_target_low.exists():

        # If no low priority targets exists we just check that the high
        # priority catalogue allocated all 17 IFUs.

        assert mangaScience_high is not None, 'no low priority target files found'
        assert len(mangaScience_high.mangaInput) == 17,\
            'no low priority target files found and not enough high priority targets allocated'

        log.info('no low priority targets found. Using only high priority targets.')
        mangaScience = mangaScience_high
        mangaScience_low = None

    mangaScience.write(toRepo=False)

    # Checks how many special targets survive
    # try:
    #     mangaScience_path = mangaScience.getDefaultFilename()
    #     special = utils.calculate_special_targets(mangaScience_path, sciPath)
    #     if special is None or special[1] == 0:
    #         log.important('no special targets found for this design.')
    #     else:
    #         survived, total = special
    #         log.important('{0} out of {1} special targets survived '
    #                       'decollision.'.format(survived, total))
    # except Exception as ee:
    #     log.warning('failed determining how many special targets '
    #                   'survived for {0}: {1}.'.format(mangaScience_path, ee), GohanUserWarning)

    stdCat = read_catalogue(os.path.join(platerun_dir,
                                         'inputs/{0}_standards.fits'.format(fieldName)))

    mangaStandard = PlateInput(designID, 'standard',
                               catalogues=stdCat,
                               surveyMode='apogeeLead', plateRun=plateRun,
                               silentOnCollision=False, sort=True,
                               decollidePlateInputs=[mangaScience],
                               raCen=raCen, decCen=decCen,
                               rejectTargets=reject_standard,
                               plotIFUs=False,
                               FOV=tileRad)
    mangaStandard.write(toRepo=False)

    skyRaw = os.path.join(readPath(config['skiesPath']), 'sky_{0}.fits'.format(fieldName))
    assert(os.path.exists(skyRaw))

    result = selectSkies(skyRaw, designID, fieldName, raCen, decCen, raise_error=False,
                         FOV=tileRad)

    if isinstance(result, table.Table):
        mangaSky = PlateInput(designID, 'sky', catalogues=result,
                              surveyMode='apogeeLead', plateRun=plateRun,
                              silentOnCollision=True, sort=False,
                              raCen=raCen, decCen=decCen)
        mangaSky.write(toRepo=False)
    else:
        assert result[0] is False, 'cannot parse the result of selectSkies'
        log.warning(result[1], GohanUserWarning)
        mangaID, ifuDesign = result[2:]

        if ifuDesign < 1000:
            reject_standard.append(mangaID)
            log.debug('added {0} to the reject_standard list.'.format(mangaID))
        else:
            reject_science.append(mangaID)
            # if mangaScience_high is not None and mangaID in mangaScience_high.getMangaIDs():
            #     reject_high.append(mangaID)
            log.debug('added {0} to the reject_science list.'.format(mangaID))

        return do_one_apogee2_manga(plateRun, platerun_dir, field,
                                    reject_science=reject_science,
                                    reject_standard=reject_standard,
                                    # reject_high=reject_high,
                                    niter=niter + 1)

    return mangaScience_low, mangaScience_high


def update_platedefinition_apogee2manga(field, plateRun, platerun_dir, high_priority=False):
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
        elif 'nInputs' in line:
            updatedPlateDefinition.append('nInputs 6')
        elif 'plateInput3' in line:
            updatedPlateDefinition.append(
                line.replace('plateInput3', 'plateInput6'))
        elif 'plateInput2' in line and not high_priority:
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
        elif 'plateInput2' in line and high_priority:
            updatedPlateDefinition.append(
                line.replace('plateInput2', 'plateInput3'))
            updatedPlateDefinition.append(
                'plateInput4 manga/{0}/mangaStandard_{1}_{2}.par'
                .format(plateRun, fieldName, designID))
            updatedPlateDefinition.append(
                'plateInput5 manga/{0}/mangaSky_{1}_{2}.par'
                .format(plateRun, fieldName, designID))
        elif 'plateInput1' in line and high_priority:
            updatedPlateDefinition.append(line.replace('plateInput1', 'plateInput2'))
            updatedPlateDefinition.insert(
                -1, 'plateInput1 manga/{0}/mangaScience_{1}_{2}.par'.format(plateRun,
                                                                            fieldName,
                                                                            designID))
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


def design_apogee2manga(platerun, platerun_dir, plate_data_path,
                        special=False, exclude=[], repeats=False, **kwargs):
    """Designs an APOGEE2-MaNGA platerun."""

    fields = table.Table.read(plate_data_path, format='ascii.commented_header')

    if special is True:
        log.important('not designing field. Only printing summary of special targets.')
        utils.print_special_summary(plate_data_path)
        return False

    mastar_targets = utils.getStellarLibraryTargets()
    psfmag = mastar_targets['psfmag']
    mastar_reject = mastar_targets[(psfmag[:, 1] < 14.6) | (psfmag[:, 3] < 14.8)]
    mastar_mangaid = set(mastar_reject['mangaid'].tolist())

    reject_science = mastar_mangaid

    for field in fields:

        mangaScience_low, mangaScience_high = do_one_apogee2_manga(
            platerun, platerun_dir, field, reject_science=list(reject_science),
            exclude=exclude)

        update_platedefinition_apogee2manga(field, platerun, platerun_dir,
                                            high_priority=mangaScience_high)

        if mangaScience_low is None:
            mangaScience = mangaScience_high
        elif mangaScience_high is None:
            mangaScience = mangaScience_low
        else:
            mangaScience = mangaScience_low.merge(mangaScience_high)

        if repeats is False:

            # Default psfmag limit
            psfmag_lim = 14.7

            # Determine if this is a short exposure plate
            plate_add_path = utils.get_path('plateDefinitionAddenda', designid=field['DesignID'])
            if os.path.exists(plate_add_path):
                plate_add = yanny.yanny(plate_add_path)
                if 'MANGA_exposure_time' in plate_add:
                    if int(plate_add['MANGA_exposure_time']) == 30:
                        psfmag_lim = 13.

            psfmag = mangaScience.mangaInput['psfmag']
            brightest_psfmag = numpy.min(psfmag[:, [1, 3]], axis=1)
            allocated_reject = mangaScience.mangaInput[brightest_psfmag < psfmag_lim]
            mangaids_reject = [mangaid.strip() for mangaid in allocated_reject['mangaid']]

            log.info('removing {} out of {} targets from future allocation '
                     'because their PSFMAG < {:.1f}'.format(len(mangaids_reject),
                                                            len(mangaScience.mangaInput),
                                                            psfmag_lim))

        else:

            mangaids_reject = [row['mangaid'].strip() for row in mangaScience.mangaInput
                               if row['repeats'] == 0]

            log.info('removing {} out of {} targets from future allocation '
                     'because REPEATS=0'.format(len(mangaids_reject),
                                                len(mangaScience.mangaInput)))

        reject_science = reject_science.union(set(mangaids_reject))

    # utils.print_special_summary(plate_data_path)
    update_addendas(platerun)

    return True


def design_manga_apogee2(platerun, platerun_dir, plate_data,
                         obs_date=None, target_version=None, std_path=None,
                         **kwargs):
    """Designs a platerun with MaNGA-led plates."""

    platerun_dir = pathlib.Path(platerun_dir)
    assert platerun_dir.exists(), 'platerun directory {!s} does not exist'.format(platerun_dir)

    assert obs_date is not None, 'obs_date must be specified for MaNGA-led plateruns.'
    assert re.match(r'^[0-9]{4}\-[0-9]{2}\-[0-9]{2}$', obs_date), 'invalid obs_date format.'

    target_version = target_version if target_version is not None else config['targets']['version']
    target_dir = pathlib.Path(os.path.expandvars(config['targets']['path'])) / target_version
    assert target_dir.exists(), 'invalid directory {!s}'.format(target_dir)

    if std_path is not None:
        std_path = pathlib.Path(std_path)
    else:
        std_path = platerun_dir / 'stds'
    assert std_path.exists(), 'invalid directory {!s}'.format(std_path)

    assert pathlib.Path(plate_data).exists()
    field_list = table.Table.read(plate_data, format='ascii.commented_header')

    platePlans_path = pathlib.Path(platerun_dir) / '{0}-platePlans.dat'.format(platerun)
    platePlansBlob = open(str(platePlans_path), 'w')

    reject_targets = []
    reject_standards = []

    for field in field_list:

        mangaTileID = int(field['manga_tileid'])
        locationID = int(field['location_id'])
        designID = int(field['design_id'])

        print(_color_text('\ndesignID: {0}, mangaTileID: {1}'.format(designID, mangaTileID),
                          'lightred'))

        raCen = float(field['ra'])
        decCen = float(field['dec'])

        # Checks whether the default MaNGA_targets_extNSA file exists, or the replacement file.
        sci_cat = target_dir / ('MaNGA_targets_extNSA_tiled_ancillary_{0:d}_ra.fits'
                                .format(mangaTileID))
        if not sci_cat.exists():
            sci_cat = target_dir / ('MaNGA_targets_extNSA_tiled_ancillary_{0:d}.fits'
                                    .format(mangaTileID))
            assert sci_cat.exists(), ('designID: {0}: neither the MaNGA_targets_extNSA file '
                                      'or the _ra version exist.'.format(designID))

        std_cat = std_path / 'manga_stds_{0}.fit'.format(mangaTileID)
        assert std_cat.exists(), 'cannot find standard file {!s}'.format(std_cat)

        mangaScience = PlateInput(designID, 'science', catalogues=str(sci_cat),
                                  surveyMode='mangaLead', plateRun=platerun,
                                  silentOnCollision=True, sort=False,
                                  raCen=raCen, decCen=decCen,
                                  manga_tileid=mangaTileID, locationid=locationID,
                                  rejectTargets=reject_targets, autocomplete=True)
        mangaScience.write()
        reject_targets += mangaScience.getMangaIDs()

        mangaStandard = PlateInput(designID, 'standard', catalogues=str(std_cat),
                                   surveyMode='mangaLead', plateRun=platerun,
                                   silentOnCollision=True, sort=False,
                                   decollidePlateInputs=[mangaScience],
                                   raCen=raCen, decCen=decCen,
                                   manga_tileid=mangaTileID, locationid=locationID,
                                   rejectTargets=reject_standards)
        mangaStandard.write()

        plateDefinition = PlateDefinition([mangaScience, mangaStandard])
        plateDefinition.write()

        platePlans = plateDefinition.getPlatePlans(obs_date)
        platePlansBlob.write(platePlans + '\n')

    return True


def design_mastar(platerun, platerun_dir, plate_data, obs_date=None, exclude=[], **kwargs):
    """Designs a platerun with MaStar plates."""

    reject_science = exclude
    reject_standard = exclude

    platerun_dir = pathlib.Path(platerun_dir)
    assert platerun_dir.exists(), 'platerun directory {!s} does not exist'.format(platerun_dir)

    assert obs_date is not None, 'obs_date must be specified for MaNGA-led plateruns.'
    assert re.match('^[0-9]{4}\-[0-9]{2}\-[0-9]{2}$', obs_date), 'invalid obs_date format.'

    field_list = table.Table.read(plate_data, format='ascii.commented_header')

    platePlans_path = pathlib.Path(platerun_dir) / '{0}-platePlans.dat'.format(platerun)
    platePlansBlob = open(str(platePlans_path), 'w')

    nn = 0
    niter = 1

    while nn < len(field_list):

        field = field_list[nn]

        designID = int(field['design_id'])
        fieldName = field['field_name']

        log.important(_color_text('designID: {0}, fieldName: {1} (iteration {2})'
                                  .format(designID, fieldName, niter), 'lightred'))

        raCen = float(field['RA'])
        decCen = float(field['Dec'])

        science_target = pathlib.Path(platerun_dir) / f'inputs/{fieldName}_target.fits'

        mangaScience = PlateInput(designID, 'science',
                                  catalogues=read_catalogue(science_target),
                                  surveyMode='MaStar', plateRun=platerun,
                                  silentOnCollision=True, sort=False,
                                  raCen=raCen, decCen=decCen,
                                  rejectTargets=reject_science,
                                  fieldName=fieldName,
                                  manga_tileid=-999, locationid=-999)

        mangaScience.write(toRepo=False)

        standards_target = pathlib.Path(platerun_dir) / f'inputs/{fieldName}_standards.fits'

        mangaStandard = PlateInput(designID, 'standard',
                                   catalogues=read_catalogue(standards_target),
                                   surveyMode='MaStar', plateRun=platerun,
                                   silentOnCollision=False, sort=True,
                                   decollidePlateInputs=[mangaScience],
                                   raCen=raCen, decCen=decCen,
                                   rejectTargets=reject_standard,
                                   plotIFUs=False, manga_tileid=-999,
                                   fieldName=fieldName,
                                   locationid=-999)
        mangaStandard.write(toRepo=False)

        skyRaw = os.path.join(readPath(config['skiesPath']), 'sky_{0}.fits'.format(fieldName))
        assert(os.path.exists(skyRaw))

        result = selectSkies(skyRaw, designID, fieldName, raCen, decCen,
                             raise_error=False, use_apogee=False, limitTo=40)

        if isinstance(result, table.Table):
            mangaSky = PlateInput(designID, 'sky', catalogues=result,
                                  surveyMode='MaStar', plateRun=platerun,
                                  silentOnCollision=True, sort=False,
                                  raCen=raCen, decCen=decCen, manga_tileid=-999,
                                  fieldName=fieldName, locationid=-999)
            mangaSky.write(toRepo=False)
        else:
            assert result[0] is False, 'cannot parse the result of selectSkies'
            log.warning(result[1], GohanUserWarning)
            mangaID, ifuDesign = result[2:]

            if ifuDesign < 1000:
                reject_standard.append(mangaID)
                log.debug('added {0} to the reject_standard list.'.format(mangaID))
            else:
                reject_science.append(mangaID)
                # if mangaScience_high is not None and mangaID in mangaScience_high.getMangaIDs():
                #     reject_high.append(mangaID)
                log.debug('added {0} to the reject_science list.'.format(mangaID))

            continue

        plateDefinition = PlateDefinition([mangaScience, mangaStandard, mangaSky])
        plateDefinition.write()

        platePlans = plateDefinition.getPlatePlans(obs_date)
        platePlansBlob.write(platePlans + '\n')

        nn += 1
        niter = 0

    return True


def generate_designs(platerun, plate_data, plate_type=None, **kwargs):
    """Generates plate inputs, plate definitions, and platePlans lines for a list of plates.

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

    if plate_type is None:
        plate_type = platerun.split('.')[-1]

    assert plate_type in ['manga-apogee2', 'apogee2-manga', 'manga', 'mastar'], 'invalid plateType'

    platerun_dir = os.path.join(os.environ['MANGA_SANDBOX'], 'platedesign/plateruns', platerun)
    assert os.path.exists(platerun_dir), 'cannot find platerun dir {0}'.format(platerun_dir)

    assert os.path.exists(plate_data), 'cannot find path {0}'.format(plate_data)

    log.info('Plate Run: {0}'.format(platerun))
    log.info('Plate Type: {0}'.format(plate_type))
    log.info('Plate Data file: {0}'.format(plate_data))
    log.info('Working on platerun directory: {0}'.format(platerun_dir))

    if plate_type == 'apogee2-manga':
        keep_log = design_apogee2manga(platerun, platerun_dir, plate_data, **kwargs)
    elif plate_type == 'manga-apogee2':
        keep_log = design_manga_apogee2(platerun, platerun_dir, plate_data, **kwargs)
    elif plate_type == 'mastar':
        keep_log = design_mastar(platerun, platerun_dir, plate_data, **kwargs)

    if keep_log:
        log_path = os.path.join(platerun_dir, platerun + '.log')
        log.info('copying log to {0}'.format(log_path))
        log.save_log(log_path)

    return
