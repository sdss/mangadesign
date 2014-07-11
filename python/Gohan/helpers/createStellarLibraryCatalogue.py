#!/usr/bin/env python
# encoding: utf-8
"""
createStellarLibraryCatalogue.py

Created by José Sánchez-Gallego on 6 Jul 2014.
Licensed under a 3-clause BSD license.

Revision history:
    6 Jul 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import sys
import glob
from astropy import table
import os


def createStellarLibraryCatalogue(catID, outName, files=[]):

    if len(files) == 0:
        files = glob.glob('./*.cat')

    names = ['RA', 'Dec', 'PriorityType', 'FieldName']
    stellTable = table.Table(None, names=names,
                             dtype=[float, float, 'S10', 'S10'])

    for file in files:
        tmpTable = table.Table.read(file, format='ascii', names=names[0:3])
        fieldName = os.path.basename(file).split('_')[0]
        for row in tmpTable:
            stellTable.add_row(list(row) + [fieldName.strip()])

    mangaID = ['{0:d}-{1:d}'.format(int(catID), nn+1)
               for nn in range(len(stellTable))]
    stellTable.add_column(table.Column(mangaID, name='mangaid'), 0)

    if os.path.exists(outName):
        os.remove(outName)

    stellTable.write(outName)


if __name__ == '__main__':
    createStellarLibraryCatalogue(sys.argv[1], sys.argv[2], sys.argv[3:])
