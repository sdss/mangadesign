import pywcsgrid2 as pywcsgrid2
import colorama
from getSDSSRun import getSDSSRun
from .runAll import runAll
from .Skies import Skies

try:
    from Staralt import Staralt
except:
    # print '{0}[WARNING]: {1}no Staralt functionality available.'.format(
    #     colorama.Fore.YELLOW, colorama.Style.RESET_ALL)
    Staralt = None

__all__ = ['loadDictionary', 'loadTemplateTable',
           'Staralt', 'pywcsgrid2', 'colorama', 'getSDSSRun',
           'runAll', 'readPath', 'Skies']
