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
import os
from astropy.io import ascii
from astropy import table
import numpy as np


def createAPASSCatalogue(catID, outName, cat):

    catTable = ascii.read(cat, names=['RA', 'DEC', 'G_UNCOR', 'R_UNCOR',
                                      'I_UNCOR', 'EBV'],
                          format='no_header')

    psf = []
    for row in catTable:
        psf.append([0.0, row['G_UNCOR'], row['R_UNCOR'], row['I_UNCOR'], 0.0])

    catTable.add_column(table.Column(psf, 'PSFMAG'))
    # catTable.add_column(table.Column(np.arange(1, len(catTable)+1),
    #                                  'PRIORITY'))

    mangaIDs = ['{0}-{1}'.format(catID, ii) for ii in
                range(1, 1+len(catTable))]
    catTable.add_column(table.Column(mangaIDs, name='MANGAID', dtype='S10'),
                        index=0)
    catTable.add_column(table.Column(len(catTable) * [2**25],
                                     name='MANGA_TARGET2', dtype=int))

    extData = []
    extFactors = np.array([0.0, 3.793, 2.751, 2.086, 0.0])
    for ebv in catTable['EBV']:
        extData.append(extFactors * ebv)

    catTable.add_column(table.Column(extData, name='EXTINCTION'))

    if os.path.exists(outName):
        os.remove(outName)

    catTable.write(outName)

    print('\nRemember to add\n\n%s %s\n\nto catalog_ids.dat!\n' %
          (catID, outName))


if __name__ == '__main__':
    createAPASSCatalogue(sys.argv[1], sys.argv[2], sys.argv[3])
