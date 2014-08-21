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

    names = ['RA', 'DEC', 'PRIORITYTYPE', 'STDTYPE', 'SUBTYPE', 'PRIORITYSUB',
             'MANGA_TARGET2']
    stellTable = table.Table(None, names=names,
                             dtype=[float, float, 'S10', 'S3', 'S1', int,
                                    int])

    for file in files:
        tmpTable = table.Table.read(file, format='ascii',
                                    names=['RA', 'DEC', 'PRIORITYTYPE'])

        priorityType = tmpTable['PRIORITYTYPE']
        stdType = [std.split('_')[0] for std in priorityType]
        subType = [std.split('_')[1][0] for std in priorityType]
        priorityInSubType = [int(std.split('_')[1][1:])
                             for std in priorityType]

        tmpTable.add_column(table.Column(stdType, 'STDTYPE'))
        tmpTable.add_column(table.Column(subType, 'SUBTYPE'))
        tmpTable.add_column(table.Column(priorityInSubType, 'PRIORITYSUB'))

        mangaTarget2 = []
        for ss in stdType:
            if 'OPT' in ss:
                maskbit = 2
            elif 'NIR' in ss:
                maskbit = 3
            elif 'PAR' in ss:
                maskbit = 4
            else:
                raise ValueError('Priority type not recognised.')
            mangaTarget2.append(2**maskbit)

        tmpTable.add_column(table.Column(mangaTarget2, 'MANGA_TARGET2'))

        tmpTable.sort(['MANGA_TARGET2', 'SUBTYPE', 'PRIORITYSUB'])

        for nn, row in enumerate(tmpTable):
            stellTable.add_row(list(row))

    mangaID = ['{0:d}-{1:d}'.format(int(catID), nn+1)
               for nn in range(len(stellTable))]
    stellTable.add_column(table.Column(mangaID, name='MANGAID'), 0)

    if os.path.exists(outName):
        os.remove(outName)

    stellTable.write(outName)


if __name__ == '__main__':
    createStellarLibraryCatalogue(sys.argv[1], sys.argv[2], sys.argv[3:])
