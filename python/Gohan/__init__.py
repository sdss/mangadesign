
def readPath(path):
    """Get a Totoro-formatted path and returns the real path.

    Paths are expanded depending on the first characher of the input string.
    If the first character is '+', the path is considered to be
    Totoro-internal relative to the root of the package. Otherwise, the path
    is expanded using os.path.expandvars and os.path.expanduser. So,
    environment variables and user tildes '~' are valid in any path.

    """

    if path[0] == '+':
        return os.path.realpath(
            os.path.join(
                os.path.dirname(__file__), path[1:]))

    else:
        return os.path.realpath(os.path.expanduser(os.path.expandvars(path)))


from astropy.wcs import FITSFixedWarning
# from astropy.io import fits

import warnings
# warnings.simplefilter('ignore', fits.file.AstropyUserWarning)
# warnings.simplefilter('ignore', AstropyDeprecationWarning)
warnings.simplefilter('ignore', FITSFixedWarning)
warnings.filterwarnings(
    'ignore',
    'This figure includes Axes that are not compatible with tight_layout')
warnings.filterwarnings('ignore', 'Module argparse was already imported')

from Gohan.exceptions import GohanWarning
warnings.filterwarnings('always', category=GohanWarning)

# from .utils.readPath import readPath
import os

__DEFAULT_CONFIG_FILE__ = readPath('+defaults.yaml')
__GOHAN_CONFIG_PATH__ = readPath('~/.gohan/gohan.yaml')

# Reads the configuration file
from Gohan.core.configuration import GohanConfig
config = GohanConfig(__DEFAULT_CONFIG_FILE__)

if os.path.exists(__GOHAN_CONFIG_PATH__):
    # If a custom configuration file exists, updates default values.
    config.update(GohanConfig(__GOHAN_CONFIG_PATH__))

# Creates the custom logging system
from Gohan.core.logger import initLog
log = initLog()

from Gohan.PlateInput import PlateInput
from Gohan.PlateMags import PlateMags
from Gohan.PlateDefinition import PlateDefinition
from Gohan.PlateTargets import PlateTargets

__all__ = ['PlateInput', 'PlateMags', 'PlateTargets',
           'log', 'config', 'PlateDefinition']
