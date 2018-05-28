#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2018-05-26
# @Filename: __init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# @Copyright: José Sánchez-Gallego

import os
import warnings

from astropy.wcs import FITSFixedWarning

from Gohan.core.configuration import get_config
from Gohan.core.logger import log
from Gohan.exceptions import GohanWarning


def readPath(path):
    """Expands a path.

    Paths are expanded depending on the first characher of the input string.
    If the first character is '+', the path is considered to be
    Gohan-internal relative to the root of the package. Otherwise, the path
    is expanded using os.path.expandvars and os.path.expanduser. So,
    environment variables and user tildes '~' are valid in any path.

    """

    if path[0] == '+':
        return os.path.realpath(
            os.path.join(
                os.path.dirname(__file__), path[1:]))

    else:
        return os.path.realpath(os.path.expanduser(os.path.expandvars(path)))


# warnings.simplefilter('ignore', fits.file.AstropyUserWarning)
# warnings.simplefilter('ignore', AstropyDeprecationWarning)
warnings.simplefilter('ignore', FITSFixedWarning)
warnings.filterwarnings('ignore',
                        'This figure includes Axes that are not compatible with tight_layout')
warnings.filterwarnings('ignore', 'Module argparse was already imported')

warnings.filterwarnings('always', category=GohanWarning)


# Reads the configuration file
config = get_config('~/.gohan/gohan.yaml')

# Creates the custom logging system
log.start_file_logger('.gohan/gohan')
log.sh.setLevel(10)
