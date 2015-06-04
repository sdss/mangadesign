#!/usr/bin/env python
# encoding: utf-8
"""
selectFieldsJun2015.py

Created by José Sánchez-Gallego on 16 May 2015.
Licensed under a 3-clause BSD license.

Revision history:
    16 May 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.scheduler import Planner
import cPickle


startDate = None
# endDate = 2457315.002083  # October 19th
endDate = 2457285.830  # 2015-09-20

planner = Planner(startDate=None, endDate=endDate, usePlatesNotAtAPO=True)
planner.schedule(goodWeatherFraction=1.0, efficiency=0.85)

fileOut = open('plannerJun2015.pckl', 'w')
cPickle.dump(planner, fileOut)
fileOut.close()
