#!/usr/bin/env python
# encoding: utf-8
"""
createPlateMags.py

Created by José Sánchez-Gallego on 17 Oct 2014.
Licensed under a 3-clause BSD license.

Revision history:
    17 Oct 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from Gohan import PlateMags
import glob


mangaScienceFiles = sorted(glob.glob('./mangaScience*.par'))

for mangaScience in mangaScienceFiles:
    plateMags = PlateMags(mangaScience)
    plateMags.write()
    plateMags.plot(useRepo=False)
