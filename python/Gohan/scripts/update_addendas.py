#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2019-12-06
# @Filename: update_addendas.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os
import pathlib
import re

import astropy.table
from sdss_access import SDSSPath


def update_addendas(platerun):
    """Updates plateDefinitionAddenda files for a given platerun.

    Assumes that the input files contain files in the form
    ``{fieldName}_target_low.fits`` which includes an ``EXPTIME`` column.

    Parameters
    ----------
    platerun : str
        Either the platerun (it expects to find a directory called
        ``$PLATERUNS_DIR/{platerun}/inputs) or the full path to the
        platerun directory.

    """

    if os.path.exists(platerun):
        path = pathlib.Path(platerun)
    else:
        path = pathlib.Path(os.environ['PLATERUNS_DIR']) / platerun

    plate_data_file = path / f'plate_data_{platerun}.dat'
    assert plate_data_file.exists(), f'cannot find {plate_data_file!s}'

    plate_data = astropy.table.Table.read(plate_data_file, format='ascii')

    sdss_path = SDSSPath()

    for row in plate_data:
        design_id = row['DesignID']
        field_name = row['FieldName']

        input_file = path / 'inputs' / f'{field_name}_target_low.fits'
        assert input_file.exists(), f'cannot find input file for field {field_name}'

        input_data = astropy.table.Table.read(input_file)
        exp_time = input_data['EXPTIME'][0]

        plate_addenda = sdss_path.full('plateDefinitionAddenda', designid=design_id)
        assert os.path.exists(plate_addenda), f'cannot find file {plate_addenda}'

        plate_addenda_data = open(plate_addenda).read().rstrip()
        if 'MANGA_exposure_time' not in plate_addenda_data:
            plate_addenda_data += '\n\nMANGA_exposure_time ' + str(int(exp_time)) + '\n\n'
        else:
            plate_addenda_data = re.sub(r'MANGA_exposure_time [0-9]+\.*[0-9]*',
                                        r'MANGA_exposure_time ' + str(int(exp_time)),
                                        plate_addenda_data)
            plate_addenda_data += '\n'

        with open(plate_addenda, 'w') as ff:
            ff.write(plate_addenda_data)
