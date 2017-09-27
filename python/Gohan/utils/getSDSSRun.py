#!/usr/bin/env python
# encoding: utf-8
"""
getSDSSRun.py

Created by José Sánchez-Gallego on 24 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    24 Feb 2014 J. Sánchez-Gallego
      Initial version

"""

import sys
import numpy as np

from six.moves.html_parser import HTMLParser

from six.moves.urllib.request import urlopen


URL = 'http://data.sdss3.org/fields/raDec?ra={0:.4f}&dec={1:.4f}'


def getSDSSRun(ra, dec):
    """
    Returns run, rerun, camcol and field for a certain RA and Dec.

    Parameters
    ----------
    RA : float or list of floats
        Right ascension.
    Dec : float or list of floats
        Declination of the field.

    Returns
    -------
    result : tuple or list of tuples.
        A tuple (run, rerun, camcol, field) for the SDSS imaging field
        matching the input coordinates or None if the coordinates are not
        in the footprint.

    """

    if isinstance(ra, (list, tuple, np.ndarray)):
        if not isinstance(dec, (list, tuple, np.ndarray)) or \
                len(ra) != len(dec):
            raise ValueError('length of ra != lenght of dec.')

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    result = []

    for ii in range(len(ra)):
        aa = ra[ii]
        dd = dec[ii]

        result.append(getPlate(aa, dd))

    if len(ra) == 1:
        return result[0]
    else:
        return result


def getPlate(aa, dd):

    url = URL.format(aa, dd)
    ff = urlopen(url)

    parser = TableParser()
    parser.feed(ff.read().decode('utf-8'))

    result = [getValue(parser._dd, value)
              for value in ['run', 'rerun', 'camcol', 'field']]

    if None in result:
        return None
    else:
        return tuple(map(int, result))


def getValue(table, value):
    for ii, td in enumerate(table):
        if td.strip() == value:
            return table[ii + 1]
    return None


class TableParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.in_td = False

    def feed(self, text):
        self._dd = []
        HTMLParser.feed(self, text)

    def handle_starttag(self, tag, attrs):
        if tag == 'dt' or tag == 'dd':
            self.in_td = True

    def handle_data(self, data):
        if self.in_td:
            self._dd.append(data)

    def handle_endtag(self, tag):
        self.in_td = False


if __name__ == '__main__':
    getSDSSRun(float(sys.argv[1]), float(sys.argv[2]))
