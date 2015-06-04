#!/usr/bin/env python
# encoding: utf-8
"""
createFieldList.py

Created by José Sánchez-Gallego on 28 Apr 2015.
Licensed under a 3-clause BSD license.

Revision history:
    28 Apr 2015 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
from sdss.internal.manga.Totoro.dbclasses import Plate
from sdss.internal.manga.Totoro import TotoroDBConnection, site, utils, config
from Gohan import utils
import cPickle
from matplotlib import pyplot as plt
from astropy import table
import os
import numpy as np


drillRun = 'Jun2015'


plt.ioff()

db = TotoroDBConnection()
session = db.Session()

tilingCat = table.Table.read(os.environ['MANGA_TILING_CATALOGUE'])
eff = config['planner']['efficiency']

fieldList = table.Table(None, names=['manga_tileid', 'location_id',
                                     'design_id', 'RA', 'Dec'],
                        dtype=[int, int, int, float, float])


def plotFields(plates, fields):

    raPlate = []
    decPlate = []
    for plate in plates:
        if not plate.isMock:
            plateFromDB = Plate(plate.pk)
            if not plateFromDB.isComplete and plateFromDB.priority > 2:
                raPlate.append(plateFromDB.ra)
                decPlate.append(plateFromDB.dec)
        else:
            raPlate.append(plate.ra)
            decPlate.append(plate.dec)

    raField = [field.ra for field in fields
               if len(field.getTotoroExposures()) > 0]
    decField = [field.dec for field in fields
                if len(field.getTotoroExposures()) > 0]

    fig, ax = plt.subplots()

    ax.scatter(tilingCat['RA'], tilingCat['DEC'], marker='x', color='gray')
    ax.scatter(raPlate, decPlate, marker='x', color='r', label='drilled')
    ax.scatter(raField, decField, marker='x', color='b',
               label='selected fields')

    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('Dec [deg]')

    plt.legend()
    plt.savefig('fieldDistribution.pdf')

    return


def plotHistograms(plates, fields):

    binSize = 0.125

    ranges = np.arange(0.0, 24.1, binSize)
    histo = np.zeros(len(ranges)-1)
    histoPlates = np.zeros(len(ranges)-1)
    histoFields = np.zeros(len(ranges)-1)

    for night in blocks:

        lst0 = site.localSiderealTime(night['JD0'])
        lst1 = site.localSiderealTime(night['JD1'])

        for ii in range(len(ranges)-1):

            lstRange = (ranges[ii], ranges[ii+1])
            ll = utils.getIntervalIntersectionLength(
                (lst0, lst1), lstRange, wrapAt=24.)

            if ll > 0.0:
                histo[ii] += ll

    for plate in plates:
        exposures = plate.getTotoroExposures(onlySets=True)
        for exp in exposures:
            lst0 = site.localSiderealTime(exp.getJD()[0])
            lst1 = lst0 + (site.localSiderealTime(exp.getJD()[1]) - lst0)
            for ii in range(len(ranges)-1):
                lstRange = (ranges[ii], ranges[ii+1])
                ll = utils.getIntervalIntersectionLength(
                    (lst0, lst1), lstRange, wrapAt=24.)

                histoPlates[ii] += ll

    for field in fields:
        exposures = field.getTotoroExposures(onlySets=True)
        for exp in exposures:
            lst0 = site.localSiderealTime(exp.getJD()[0])
            lst1 = lst0 + (site.localSiderealTime(exp.getJD()[1]) - lst0)
            for ii in range(len(ranges)-1):
                lstRange = (ranges[ii], ranges[ii+1])
                ll = utils.getIntervalIntersectionLength(
                    (lst0, lst1), lstRange, wrapAt=24.)

                histoFields[ii] += ll

    print('Number of hours available: {0:.1f}\nNumber of hours used: {1:.1f}'
          .format(np.sum(histo), np.sum(histoPlates)+np.sum(histoFields)))

    plt.clf()
    plt.cla()

    fig, ax = plt.subplots()

    ax.bar(ranges[0:-1], histo, binSize, edgecolor='None', color='y',
           linewidth=0.0, alpha=0.5, label='Available time')
    ax.bar(ranges[0:-1], histoFields+histoPlates, binSize, color='b',
           edgecolor='None', linewidth=0.0, alpha=0.4,
           label='Observed time (selected)')
    # ax.bar(ranges[0:-1], histoPlates, binSize, color='r', edgecolor='None',
    #        linewidth=0.0, alpha=0.4, label='Observed time (drilled)')

    ax.legend(frameon=True, loc='upper left')
    ax.set_xlim(0.0, 24)
    # ax.set_ylim(0.0, 10)
    ax.set_xlabel('LST [hours]')
    ax.set_ylabel('Time [hours]')
    plt.savefig('allocatedTime.pdf')

    return


def createFieldList(fields, drillRun):

    lastLocationID = utils.getLastLocationID()
    locationID = lastLocationID + 1

    for field in fields:
        if len(field.getTotoroExposures()) == 0:
            continue
        fieldList.add_row((field.manga_tileid, locationID, -1, field.ra,
                           field.dec))
        locationID += 1

    fieldList.write('fieldList{0}.dat'.format(drillRun),
                    format='ascii.commented_header')


planner = cPickle.load(open('planner{0}.pckl'.format(drillRun), 'r'))
blocks = planner.blocks

plates = planner.plates
fields = planner.fields

plotFields(plates, fields)
plotHistograms(plates, fields)
createFieldList(fields, drillRun)
