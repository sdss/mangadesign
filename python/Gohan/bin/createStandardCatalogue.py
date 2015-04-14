#!/usr/bin/env python
# encoding: utf-8
"""
createStandardCatalogue.py

Created by José Sánchez-Gallego on 21 May 2014.
Licensed under a 3-clause BSD license.

Revision history:
    21 May 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import sys
import glob
import os
from astropy import table
from astropy.io import ascii
import numpy as np


def processField(tt, tileID, catID, idx=1):

    tNew = tt.copy()

    psfColnames = [col for col in tt.colnames if '_psf' in col]
    tNew.remove_columns(psfColnames)

    psfs = [row._data for row in tt[psfColnames]]
    tNew.add_column(table.Column(data=psfs, name='psfmag', dtype='5f8'),
                    index=3)

    sdss = tt['SDSS']
    tNew.remove_column('SDSS')
    tNew.add_column(sdss)

    try:
        tNew.add_column(table.Column(data=len(tt) * [tileID],
                                     name='tileid', dtype=int), index=0)
    except ValueError:
        tNew.add_column(table.Column(data=len(tt) * [tileID],
                                     name='tileid', dtype='S10'), index=0)

    tNew.add_column(table.Column(data=len(tNew) * [2**23],
                                 name='manga_target2',
                                 dtype=int))

    tNew.add_column(table.Column(data=np.arange(len(tt))+1, name='priority',
                                 dtype=int))

    mangaIDs = ['{0}-{1}'.format(catID, ii) for ii in
                range(idx, idx+len(tt))]
    tNew.add_column(table.Column(mangaIDs, name='mangaid', dtype='S10'),
                    index=0)
    tNew.add_column(table.Column(
        data=np.zeros((len(tt), 7)),
        name='extinction'))

    return tNew, idx+len(tt)


def createStandardCatalogue(catID, outName, files=[]):

    if len(files) == 0:
        files = glob.glob('./manga_standards_*.txt')

    idx = 1
    for nn, file in enumerate(files):

        tileID = os.path.splitext(file)[0].split('_')[-1]
        tt = ascii.read(file, format='commented_header')
        tNew, idx = processField(tt, tileID, catID, idx=idx)

        if nn == 0:
            std = tNew.copy()
        else:
            for row in tNew:
                std.add_row(row)

    for col in std.colnames:
        if col != col.upper():
            std.rename_column(col, col.upper())

    if os.path.exists(outName):
        os.remove(outName)

    std.write(outName)

    print('\nRemember to add\n\n%s %s\n\nto catalog_ids.dat!\n' %
          (catID, outName))


if __name__ == '__main__':
    createStandardCatalogue(sys.argv[1], sys.argv[2], sys.argv[3:])
