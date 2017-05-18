#!/usr/bin/env python
# encoding: utf-8
#
# update_cmm.py
#
# Created by José Sánchez-Gallego on 17 May 2017.


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from Gohan import log

from sdss.common import addPlateHoles2db
from sdss.internal.database.connections import APODatabaseAdminLocalConnection as db
from sdss.internal.database.apo.platedb import ModelClasses as plateDB

from fitPlugPlateMeas.cmmPlateMeas import PlateMeas

import fnmatch
import itertools
import numpy
import os
import warnings


def getCmmMeas(plate, pm, session):
    """Returns a CMM Measurement object that can be loaded to the DB.

    Modified from Petunia. ``pm`` is a fitPlugPlateMeas.PlateMeas object

    """

    def focalXY(xPos, yPos):
        # determined by least squares fitting radialFlat to radialFocal in plDrillPos files
        pp = numpy.asarray([-1.02073456e-16, 8.54520957e-14, -6.58010409e-11,
                            2.67344224e-08, -4.72650487e-06, 9.99726411e-01,
                            2.40696326e-06])
        rr = numpy.sqrt(xPos**2 + yPos**2)
        theta = numpy.arctan2(yPos, xPos)
        rNew = numpy.polyval(pp, rr)
        xNew = rNew * numpy.cos(theta)
        yNew = rNew * numpy.sin(theta)
        return xNew, yNew

    pmMeasDict = pm.export()

    cmmMeas = plateDB.CmmMeas()
    cmmMeas.plate = plate
    cmmMeas.date = pmMeasDict['measDate']
    cmmMeas.cmmfilename = pm.fileName
    cmmMeas.fitoffsetx = pmMeasDict['xyOff'][0]
    cmmMeas.fitoffsety = pmMeasDict['xyOff'][1]
    cmmMeas.fitrot = pmMeasDict['rotAngle']
    cmmMeas.fitscale = pmMeasDict['scale']
    cmmMeas.fitqpmag = pmMeasDict['quadrupoleMag']
    cmmMeas.fitqpang = pmMeasDict['quadrupoleAngle']

    plateHoleInds = []
    radialHoleErrs = []
    for xPos, yPos in itertools.izip(pmMeasDict['xPos'], pmMeasDict['yPos']):
        # discover indicies matching measured holes to dbholes
        xPos, yPos = focalXY(xPos, yPos)
        dbHolesXY = numpy.asarray([[float(hole.xfocal) - xPos, float(hole.yfocal) - yPos]
                                   for hole in plate.plateHolesFile[0].plateHole])
        holePosErr = numpy.linalg.norm(dbHolesXY, axis=1)
        matchInd = numpy.argmin(holePosErr)
        plateHoleInds.append(matchInd)
        he = holePosErr[matchInd]
        radialHoleErrs.append(he)
    if not len(plateHoleInds) == len(set(plateHoleInds)):
        # not a unique matching, something's wrong
        raise RuntimeError('Could not find unique hole matching for '
                           'plate {0:d}, measurement file {0}'.format(plate.plate_id,
                                                                      str(pm.fileName)))

    holeMeasList = []
    for measInd, plateHoleInd in enumerate(plateHoleInds):
        holeMeas = plateDB.HoleMeas()
        plateHole = plate.plateHolesFile[0].plateHole[plateHoleInd]
        # if this plate is boss, explicitly Null out apogee fields
        if not ('apogee' in ''.join([survey.label.lower() for survey in plate.surveys])):
            plateHole.apogee_target1 = None
            plateHole.apogee_target2 = None
        holeMeas.plateHole = plateHole
        holeMeas.cmmMeas = cmmMeas
        holeMeas.nomdia = pmMeasDict['nomDia'][measInd]
        holeMeas.diaerr = pmMeasDict['diaErrAll'][measInd]
        holeMeas.dia_measured = pmMeasDict['measDia'][measInd]
        holeMeas.nomx = pmMeasDict['xPos'][measInd]
        holeMeas.nomy = pmMeasDict['yPos'][measInd]
        holeMeas.measx = pmMeasDict['measx'][measInd]
        holeMeas.measy = pmMeasDict['measy'][measInd]
        holeMeas.residx = pmMeasDict['xErr'][measInd]
        holeMeas.residy = pmMeasDict['yErr'][measInd]
        holeMeas.residr = pmMeasDict['radErr'][measInd]
        holeMeas.qpresidx = pmMeasDict['qpXErr'][measInd]
        holeMeas.qpresidy = pmMeasDict['qpYErr'][measInd]
        holeMeas.qpresidr = pmMeasDict['qpRadErr'][measInd]
        holeMeasList.append(holeMeas)

    return cmmMeas, holeMeasList


def update_cmm(path):
    """Ingests new CMM files into the DB.

    Recursively checks whether there are new CMM files in ``path`` that have
    not been ingested to the DB and loads them.

    """

    session = db.Session()

    cmm_filenames = zip(*session.query(plateDB.CmmMeas.cmmfilename).all())[0]

    path = os.path.realpath(path)

    if os.path.isdir(path):
        log.info('recursively looking for CMM files in {0}'.format(path))

        matches = []
        for root, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, 'D[0-9]*_[0-9]*'):
                matches.append(os.path.join(root, filename))

        log.info('found {0} CMM files'.format(len(matches)))

        if len(matches) == 0:
            return

    else:
        assert os.path.exists(path), 'file {0} not found'.format(path)
        log.info('using file {0}'.format(path))

    new_files = sorted([match for match in matches
                        if os.path.basename(match) not in cmm_filenames])

    if len(new_files) > 0:
        log.info('found {0} new CMM files'.format(len(new_files)))
        log.info('new files: {0}'.format(', '.join(map(os.path.basename, new_files))))
    else:
        log.info('all CMM files have already been loaded to the DB.')
        return

    log.info('loading files to DB ...')

    for mf in new_files:
        pm = PlateMeas(mf)
        plate = session.query(plateDB.Plate).filter(plateDB.Plate.plate_id == pm.plateID).one()

        if len(plate.plateHolesFile) == 0:
            log.info('adding holes for plate {0}'.format(plate.plate_id))
            addPlateHoles2db.run(plateId=plate.plate_id, session=session)
            plate = session.query(plateDB.Plate).filter(plateDB.Plate.plate_id == pm.plateID).one()

        with session.begin(subtransactions=True):

            cmm_meas, holeMeasList = getCmmMeas(plate, pm, session)

            for holeMeas in holeMeasList:
                session.add(holeMeas)

            session.add(cmm_meas)

            log.debug('loaded file {0}'.format(os.path.basename(mf)))

    log.info('all files have been loaded')

    return
