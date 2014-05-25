#!/usr/bin/env python
# encoding: utf-8
"""
Staralt.py

Created by José Sánchez-Gallego on 26 Nov 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

from astropy import table
from astropy import coordinates as coo
from astropy import units
from astropy import time
import ephem
from random import randint
import pylab as plt
import matplotlib.dates as dates
import datetime as dt
import numpy as np
import brewer2mpl


TWILIGHT_ANGLE = coo.Angle(-19., unit=units.degree)
RISE_SET_ANGLE = coo.Angle(-0.56, unit=units.degree)
COLOURS = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

INTERVAL = 900.
MOON_INTERVAL = 300. * 12


def jd2lmst(jd, longitude=254.179722):

    if isinstance(jd, time.Time):
        jd = jd.jd

    jd0 = int(jd) + 0.5
    dd0 = jd0 - 2451545.
    dd = jd - 2451545.
    tt = dd / 36525
    hh = (jd - jd0) * 24.

    gmst = 6.697374558 + 0.06570982441908 * dd0 + 1.00273790935 * hh + \
        0.000026 * tt**2
    gmstL = coo.Longitude(gmst * units.hour)

    if not isinstance(longitude, coo.Longitude):
        longitude = coo.Longitude(longitude * units.degree)

    return coo.Longitude(gmstL + longitude)


class Staralt(object):

    def __init__(self, ra, dec, date, name=None, ha=0.0,
                 long=-105.820278, lat=32.780278, elev=2800,
                 haLimits='MaNGA', plotHALimits=True,
                 plotMoonDistance=True, plotMoon=True,
                 figID=None, lst=True):

        self.coords = coo.ICRS(ra=ra, dec=dec,
                               unit=(units.degree, units.degree))
        try:
            self.nObjects = len(ra)
        except:
            self.nObjects = 1

        self.haLimits = haLimits

        if isinstance(date, time.Time):
            self.date = date
        elif isinstance(date, basestring):
            self.date = time.Time(date, format='iso', scale='utc')
        else:
            raise ValueError('date format not understood.')

        self.long = long
        self.lat = lat
        self.elev = elev

        if figID is not None:
            self.figID = figID
        else:
            self.figID = randint(0, 99999)

        self.setObservatory()
        self.defineObjects(name)
        self.getStaraltData()

        self.createCanvas()

        self.plotObjectAltitudes()

        if plotHALimits:
            self.plotObjectHA()

        if plotMoonDistance:
            self.plotMoonDistance()

        if plotMoon:
            self.plotMoon()

        self.createLegend()

        if lst:
            self.plotLST()

    def save(self, filename, **kwargs):
        plt.savefig(filename, **kwargs)
        plt.close('all')

    def plotLST(self):

        tAx = [time.Time(dates.num2date(aa),
                         scale='utc', lon=self.long)
               for aa in self.ax.xaxis.get_majorticklocs()]
        lstLabels = []
        for tt in tAx:
            tt.delta_ut1_utc = 0.0
            sidereal = jd2lmst(tt.jd)
            lstLabels.append(sidereal.to_string(fields=2))

        lstLabels = ['LST'] + lstLabels
        ticks = [dates.date2num(self.sunSet.datetime)] + \
            self.ax.get_xticks().tolist()

        self.ax3.xaxis.set_ticks(ticks)
        self.ax3.xaxis.set_ticklabels(lstLabels)

        for label in self.ax3.get_xticklabels():
            label.set_fontsize(8)
        self.ax3.tick_params(axis='x', pad=20)

    def defineObjects(self, name):

        if self.nObjects > 8:
            raise ValueError('maximum of 8 objects per plot.')

        if self.nObjects == 1:
            self.objects = StaraltObject(self.coords, name=name,
                                         observatory=self.obsEphem)
            self.objects.haLimits = self.haLimits
        else:

            if name is None:
                name = self.nObjects * [None]
            if self.haLimits is None or isinstance(self.haLimits, basestring):
                self.haLimits = self.nObjects * [self.haLimits]

            self.objects = [StaraltObject(self.coords[nn], name=name[nn],
                                          observatory=self.obsEphem)
                            for nn in range(self.nObjects)]

            for nn, object in enumerate(self.objects):
                object.haLimits = self.haLimits[nn]

    def getStaraltData(self):

        self.interval = time.TimeDelta(INTERVAL, format='sec')
        self._nIntervals = int((self.sunRise - self.sunSet).sec /
                               self.interval.sec) + 1

        self.times = self.sunSet + np.arange(self._nIntervals) * self.interval

    def plotObjectAltitudes(self):

        for nn, object in enumerate(self.objects):
            object.getAltitudes(self.times)
            if object.colour is None:
                object.colour = COLOURS[nn]
            object.plotLine, = self.ax.plot(
                object.altitudes['time'],
                object.altitudes['altitudes'],
                color=object.colour, ls='solid', lw=1.0, zorder=10)

    def plotMoon(self):

        moon = ephem.Moon()

        middleNight = self.sunSet + (self.sunRise - self.sunSet) / 2.
        self.obsEphem.date = middleNight.datetime

        moon.compute(self.obsEphem)
        moonCoords = coo.ICRS(ra=moon.ra, dec=moon.dec,
                              unit=(units.radian, units.radian))

        self.moon = StaraltObject(moonCoords, name='Moon',
                                  observatory=self.obsEphem)
        self.moon.getAltitudes(self.times)
        self.moon.colour = 'k'

        self.moon.plotLine, = self.ax.plot(
            self.moon.altitudes['time'],
            self.moon.altitudes['altitudes'],
            color=self.moon.colour, ls='dashed', lw=1.0, zorder=10)

    def plotObjectHA(self):

        for nn, object in enumerate(self.objects):
            object.getHA(self.times)
            if object.haLimits is None:
                continue
            if object.colour is None:
                object.colour = COLOURS[nn]

            indices = (object.HA['HA'] > object.haLimits[0]) * \
                      (object.HA['HA'] < object.haLimits[1])
            self.ax.plot(
                object.altitudes['time'][indices],
                object.altitudes['altitudes'][indices],
                color=object.colour, ls='solid', lw=4.0, zorder=20)

    def plotMoonDistance(self):

        self.moonInterval = time.TimeDelta(MOON_INTERVAL, format='sec')
        self._nMoonIntervals = int((self.sunRise - self.sunSet).sec /
                                   self.moonInterval.sec) + 1

        self.moonTimes = self.sunSet + \
            np.arange(self._nMoonIntervals) * self.moonInterval

        for nn, object in enumerate(self.objects):
            object.getMoonDistance(self.moonTimes)
            if object.colour is None:
                object.colour = COLOURS[nn]
            for row in object.moonDistance:
                tt = row['time']
                aa = row['altitude']
                if aa < 5. or tt < (self.sunSet +
                                    time.TimeDelta(600.,
                                                   format='sec')).datetime:
                    continue
                distance = row['moonDistance']
                self.ax.text(tt, aa, int(distance), fontsize=9,
                             horizontalalignment='center',
                             verticalalignment='center',
                             bbox=dict(facecolor='white',
                                       edgecolor='None',
                                       boxstyle='round, pad=0'),
                             zorder=30, color=object.colour)

    def createLegend(self):
        if self.nObjects == 1:
            raHMS = self.coords.ra.to_string(unit=units.hour, precision=1,
                                             fields=2)
            decDMS = self.coords.dec.to_string(unit=units.degree, precision=1,
                                               fields=2)
            coords = r'$(\rm\alpha={0},\ \delta={1})$'.format(
                raHMS, decDMS)
            self.ax.set_title(self.objects.name + ' ' + coords, fontsize=11)
        else:

            lines = [object.plotLine for object in self.objects]
            labels = [object.name for object in self.objects]

            self.legend = self.ax.legend(lines, labels,
                                         labelspacing=0.7, loc='upper right',
                                         bbox_to_anchor=(1., 1.02),
                                         handletextpad=0.2, mode='expand')
            self.legend.draw_frame(False)
            plt.setp(self.legend.get_texts(), fontsize=9)

    def createCanvas(self):

        self.fig = plt.figure(self.figID)

        rightMargin = 0.9
        if self.nObjects > 1:
            rightMargin = 0.8

        self.fig.subplots_adjust(hspace=0.0, bottom=0.1,
                                 top=0.92, left=0.1, right=rightMargin)

        self.ax = self.fig.add_subplot('111')
        self.ax2 = self.fig.add_axes(self.ax.get_position(),
                                     sharey=self.ax, frameon=False)
        self.ax3 = self.fig.add_axes(self.ax.get_position(),
                                     sharey=self.ax, frameon=False)

        self.ax2.xaxis.tick_top()
        self.ax3.xaxis.tick_top()

        self.ax.cla()
        self.ax2.cla()
        self.ax3.cla()

        self.ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        self.ax.set_xlim(self.sunSet.datetime, self.sunRise.datetime)
        self.ax2.set_xlim(self.sunSet.datetime, self.sunRise.datetime)
        self.ax3.set_xlim(self.sunSet.datetime, self.sunRise.datetime)
        self.ax.set_ylim(0.0, 90.0)
        self.ax.set_xlabel('UT Starting on ' +
                           self.sunSet.datetime.strftime('%d-%m-%Y'))
        self.ax.set_ylabel('Elevation (degrees)')
        self.ax.axvline(self.sunSetTwilight.datetime, ls='--', c='k')
        self.ax.axvline(self.sunRiseTwilight.datetime, ls='--', c='k')
        self.ax.grid()

        middlePoint = (self.sunRise.datetime -
                       self.sunRise.datetime).seconds / 2.0
        middlePoint = self.sunSet.datetime + dt.timedelta(seconds=middlePoint)

        self.ax2.xaxis.set_ticks([self.sunSet.datetime,
                                  self.sunSetTwilight.datetime,
                                  self.sunRiseTwilight.datetime,
                                  self.sunRise.datetime])
        self.ax2.xaxis.set_ticklabels(['Sunset',
                                       'Ev. Tw.',
                                       'Mor. Tw.',
                                       'Sunrise'])
        self.ax2.yaxis.set_visible(False)

        for label in self.ax2.get_xticklabels():
            label.set_fontsize(8)
        self.ax.get_yticklabels()[-1].set_visible(False)

    def setObservatory(self):

        self.obsEphem = ephem.Observer()
        self.obsEphem.long = str(self.long)
        self.obsEphem.lat = str(self.lat)

        if self.elev is None:
            self.obsEphem.elevation = 0.0
        else:
            self.obsEphem.elevation = self.elev

        self.obsEphem.date = self.date.datetime
        self.obsEphem.pressure = 0
        self.computeSun()

    def computeSun(self):

        sun = ephem.Sun()
        self.obsEphem.horizon = RISE_SET_ANGLE.radian
        sun.compute(self.obsEphem)

        sunAltitude = coo.Angle(sun.alt, unit=units.radian)

        if sunAltitude.degree < 0:
            sunSetting = self.obsEphem.previous_setting(sun)
            sunRising = self.obsEphem.next_rising(sun)
        else:
            sunSetting = self.obsEphem.next_setting(sun)
            sunRising = self.obsEphem.next_rising(sun)

        self.sunSet = time.Time(self._formatISO(sunSetting),
                                format='iso', scale='utc')
        self.sunRise = time.Time(self._formatISO(sunRising),
                                 format='iso', scale='utc')

        self.computeTwilights()

    def computeTwilights(self):
        sun = ephem.Sun()
        self.obsEphem.horizon = TWILIGHT_ANGLE.radian

        self.obsEphem.date = self.sunSet.datetime
        settingTw = self.obsEphem.next_setting(sun)
        self.sunSetTwilight = time.Time(self._formatISO(settingTw),
                                        format='iso', scale='utc')

        self.obsEphem.date = self.sunRise.datetime
        risingTw = self.obsEphem.previous_rising(sun)
        self.sunRiseTwilight = time.Time(self._formatISO(risingTw),
                                         format='iso', scale='utc')

        self.obsEphem.horizon = RISE_SET_ANGLE.radian

    def _formatISO(self, value):
        return str(value).replace('/', '-')


class StaraltObject(object):

    def __init__(self, coords, name=None, observatory=None, haLimits=None):
        self.coords = coords
        self.name = name
        self.ephem = ephem.readdb(
            self.formatObjectEphem(self.coords, self.name))

        self.colour = None
        self.plotLine = None

        self.haLimits = haLimits

        if observatory is not None:
            self.observatory = observatory

    def formatObjectEphem(self, coords, name):
        return str(name) + ',f,' + \
            coords.to_string(precision=0, sep=':').replace(' ', ',') + \
            ',0.0,2000.0'

    def getAltitudes(self, obsTimes, onlyReturn=False,
                     observatory=None):
        if observatory is None:
            observatory = self.observatory
        altitudes = np.zeros(len(obsTimes))
        for nn, obsTime in enumerate(obsTimes):
            observatory.date = obsTime.datetime
            self.ephem.compute(observatory)
            altitudes[nn] = coo.Angle(self.ephem.alt,
                                      unit=units.radian).degree

        if not onlyReturn:
            self.altitudes = table.Table([obsTimes.datetime, altitudes],
                                         names=['time', 'altitudes'])

        return altitudes

    def getHA(self, obsTimes, observatory=None, onlyReturn=False):

        if observatory is None:
            observatory = self.observatory

        HA = np.zeros(len(obsTimes))
        observatory.date = obsTimes[0].datetime
        self.transit = time.Time(self._formatISO(
            observatory.next_transit(self.ephem)),
            format='iso', scale='utc')

        for nn, obsTime in enumerate(obsTimes):
            HA[nn] = (obsTime - self.transit).sec / 3600.

        if not onlyReturn:
            self.HA = table.Table([obsTimes.datetime, HA],
                                  names=['time', 'HA'])

        return HA

    def getMoonDistance(self, obsTimes, observatory=None):

        if observatory is None:
            observatory = self.observatory

        moon = ephem.Moon()
        moonDistances = np.zeros(len(obsTimes))
        altitudes = self.getAltitudes(obsTimes, onlyReturn=True)
        for nn, obsTime in enumerate(obsTimes):
            observatory.date = obsTime.datetime
            moon.compute(observatory)
            moonCoords = coo.ICRS(ra=moon.ra, dec=moon.dec,
                                  unit=(units.radian, units.radian))
            moonDistances[nn] = moonCoords.separation(self.coords).degree

        self.moonDistance = table.Table(
            [obsTimes.datetime, altitudes, moonDistances],
            names=['time', 'altitude', 'moonDistance'])

        return moonDistances

    def _formatISO(self, value):
        return str(value).replace('/', '-')

    @property
    def haLimits(self):
        return self._haLimits

    @haLimits.setter
    def haLimits(self, value):
        if value is None:
            self._haLimits = value
        elif value == 'MaNGA':
            self._haLimits = self.mangaHaLimits()
        else:
            try:
                if len(value) == 2:
                    self._haLimits = np.sorted(value)
                else:
                    raise ValueError('haLimits has length != 2')
            except:
                self._haLimits = None

    def mangaHaLimits(self):
        """
        Calculates the maximum HAs acceptable for a list of declinations.
        Uses the polynomial fit by David Law and a omega limit of 0.5.
        """

        funcFit = np.array([1.59349, 0.109658, -0.00607871,
                            0.000185393, -2.54646e-06, 1.16686e-08])[::-1]

        dec = self.coords.dec.degree

        if dec < -10 or dec > 80:
            haLim = 0.0
        else:
            haLim = np.polyval(funcFit, dec)

        return np.array([-haLim, haLim])


    def __iter__(self):
        return iter([self])

    def __len__(self):
        return 1
