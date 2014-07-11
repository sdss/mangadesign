#!/usr/bin/env python
# encoding: utf-8
"""
PlateInput.py

Created by José Sánchez-Gallego on 6 Feb 2014.
Licensed under a 3-clause BSD license.

Major revision history:
    3 Feb 2014 J. Sánchez-Gallego
      Initial version.
    13 Feb 2014 J. Sánchez-Gallego
      Modified to use InputCatalogue class as input method.
    28 Feb 2014 J. Sánchez-Gallego
      Added documentation.
    22 May 2014 J. Sánchez-Gallego
      Partially rewritten to work with the changes to InputCatalogue.

"""

from astropy import table
import numpy as np
from sdss.utilities import yanny
import os
import shutil as sh
from .exceptions import GohanUserWarning, GohanError
import warnings
# from PlateMags import PlateMags
from astropy import coordinates as coo
from astropy import units as uu
from . import config, log
from collections import OrderedDict
from . import readPath

__ALL__ = ['PLATEPLANS_TEMPLATE', 'PlateInput', 'PlateInputBase']

platePlansTemplate = """PLATEPLANS {plateID} {designID} {locationID} -1
{platedesignversion} {{ 0.0 0.00000 0.00000 0.00000 0.00000 0.00000 }}
{temp} {epoch} {raCen} {decCen} {survey} {programname} {drillstyle} \" \"
{plateRun} {chunk} \"{name}\" \"{comment}\"
"""

plateDefinitionTemplate = readPath('+etc/mangaDefinition_Default.par')


cartmap = yanny.yanny(os.path.expandvars(config['cartmap']), np=True)['CMAP']


def ifuDesign2Size(ifudesign):
    row = cartmap[cartmap['ifuDesign'] == ifudesign][0]
    if len(row) == 0:
        raise GohanError('Invalid ifudesign {0}'.format(ifudesign))
    return row['ifusize']


class PlateInput(list):
    """A class to construct plateInput files.

    This files uses an ``InputCatalogue`` instance to create an
    object with the information necessary to write a plateInput file. The
    plateInput file has two components: a structure, that correspond to the
    records in the InputCatalogue object; and a pairs element that is populated
    from the metadata in the input catalogue. Additional key/value pairs can
    be defined manually.

    PlateInput objects are lists. They can be initialised by providing one or
    more InputCatalogue instances. Each element of the resulting PlateInput
    object correspond to one plateInput file that can be saved to disk.
    Additional methos allow to create plateDefinition and platePlans files
    for the corresponding set of plateInputs.

    Note that the values defined for `pairs` and `reassignFerrules` affect to
    all the elements in then PlateInput object.

    Parameters
    ----------
    designid : int
        The designID for the plate.
    plateRun : str
        A string with the plateRun (e.g., '2014.02.x.manga').
    catalogues : a single ``InputCatalogue`` or a path or a list of them
        The input catalogue data, either as a single InputCatalogue
        or as a list of InputCatalogue instances. If strings (or list of
        strings) are provided, InputCatalogues will be read from them.
    plateType : str, optional
        The type of plate to be designed. By default it assumes MaNGA leading
        ('mangaLeading').
    pairs : dict, optional
        Additional keyword/value pairs to be added to the plateInput files.
        Note that if a key is defined in pairs, it supersedes the metadata
        for each input catalogue.
    reassignFerrules : bool
        If True, ifudesigns are automatically assigned to optimise their
        distribution per anchor block. This will still happen for each target
        for which ifudesign is undefined.
    verbose : bool, optional
        If False, only prints warnings. If None, only errors are raised.

    """

    def __init__(self, designid, plateRun, catalogues,
                 plateType='mangaLeading', pairs={}, reassignFerrules=False,
                 verbose=True, **kwargs):

        log.setVerbose(verbose)

        self.designid = designid
        self.plateType = plateType
        self.plateRun = plateRun
        self.kwargs = kwargs

        catalogues = np.atleast_1d(catalogues)

        inputs = [
            PlateInputBase(self.designid, catalogues[ii],
                           pairs=pairs,
                           reassignFerrules=reassignFerrules,
                           plateRun=self.plateRun)
            for ii in range(len(catalogues))
        ]

        list.__init__(self, inputs)

    @property
    def struct1(self):
        """Returns a stacked table with every struct1 table in self."""

        if len(self) == 1:
            return self[0].struct1
        return table.vstack([pp.struct1 for pp in self],
                            metadata_conflicts='silent')

    @struct1.setter
    def struct1(self, value):
        raise ValueError('It is not possible to set struct1.')

    # def getPlateMags(self, **kwargs):
    #     """Returns a PlateMags instance from this plateInput instance."""
    #     return PlateMags(self, **kwargs)

    def plotIFUs(self, filename=None, **kwargs):
        """ Plots the location of the IFUs in the plate.

        Parameters
        ----------
        filename : str, optional
            The filename of the output plot. If not defined, the default
            filename is plateIFUs-XXXX.pdf with XXXX the designID.
        kwars
            Other parameters to be passed to ``PlateInputBase.plotIFUs```.

        """

        tileid = self[0].meta['tileid'] if 'tileid' in self[0].meta else 0

        if filename is None:
            filename = 'plateIFUs_{0:04d}_{1:04d}.pdf'.format(
                int(tileid), self.designid)

        self[0].plotIFUs(filename, struct1=self.struct1, **kwargs)

    def check(self):
        """
        To be written. This method will provide a number of sanity checks to
        the plateInput data, such as checking that the RA and Dec centres are
        the same for all plateInput files, checking ifuDesign and ifuSize, etc.
        """

    def write(self, toRepo=False, appendDefinition=False, platePlans=True):
        """Writes the plateInput, platePlan and plateDefinition files.

        This methods creates the plateInput, plateDefinition and platePlans
        files with the appropriate filenames.

        Parameters
        ----------
        toRepo : bool, optional
            Copies the plateInput files and plateDefinition to the repository.
        appendDefinition : bool, optional
            If the plateDefinition file exists in the repo, only appends the
            plateInput files, but does not overwrite other information. Note
            that this parameter is only relevant if copyToRepo or moveToRepo
            are True. [TO BE IMPLEMENTED]
        platePlans : bool, optional
            If True, write also the platePlans line to a text file.

        """

        plateInputFilenames = []
        for plateInput in self:
            filename = plateInput.write(toRepo=toRepo)
            plateInputFilenames.append(filename)

        plateDefinitionFilename = self.writePlateDefinition(
            toRepo=toRepo, append=appendDefinition)

        platePlansFilename = None
        if platePlans:
            platePlansFilename = self.writePlatePlans()

        return plateInputFilenames, plateDefinitionFilename, platePlansFilename

    def getPlateDefinitionYanny(self):
        """Return the platePlans text."""

        raCen = self[0].pairs['racen']
        decCen = self[0].pairs['deccen']
        nInputs = len(self)
        priority = ' '.join([str(ii+1) for ii in range(nInputs)])

        defDict = OrderedDict(
            [['raCen', raCen], ['decCen', decCen], ['nInputs', nInputs],
             ['priority', priority], ['designID', self.designid]])

        defDict['plateType'] = \
            config['plateTypes'][self.plateType]['plateType']
        defDict['plateLead'] = \
            config['plateTypes'][self.plateType]['plateLead']
        defDict['platedesignversion'] = \
            config['plateTypes'][self.plateType]['platedesignversion']
        defDict['defaultSurveyMode'] = \
            config['plateTypes'][self.plateType]['defaultSurveyMode']

        inputs = []
        for nn, plateInput in enumerate(self):

            if plateInput.filename is None:
                raise GohanError('filename not set for plateInput with ' +
                                 'designID={0}. '.format(plateInput.designID) +
                                 'Write it to disk or set the filename ' +
                                 'before writing the plateDefinition file.')

            inputFilename = 'manga/{0}/{1}'.format(self.plateRun,
                                                   self[nn].filename)
            inputs += [['plateInput{0:d}'.format(nn+1), inputFilename]]

        inputs = OrderedDict(inputs)
        defDict.update(inputs)

        template = yanny.yanny(plateDefinitionTemplate)

        for key in defDict:
            template[key] = defDict[key]

        return template

    def writePlateDefinition(self, toRepo=False, append=False):

        definitionYanny = self.getPlateDefinitionYanny()

        filename = 'plateDefinition-{0:06d}.par'.format(self.designid)
        definitionYanny.set_filename(filename)

        if os.path.exists(filename):
            os.remove(filename)

        definitionYanny.write()
        log.info('plateDefinition file {0} saved.'.format(filename))

        if toRepo:
            pathRoot = '{0:04d}XX/'.format(int(self.designid / 100))
            plateDefinitionPath = os.path.join(
                os.path.expandvars(config['platelist']),
                'definitions', pathRoot)

            if not os.path.exists(plateDefinitionPath):
                os.makedirs(plateDefinitionPath)

            sh.copy(filename, os.path.join(plateDefinitionPath, filename))
            log.info('plateDefinition file {0} '.format(filename) +
                     'saved to platelist repo.')

        return filename

    def getPlatePlans(self):
        """Returns the platePlans text."""

        plansDic = {}

        plansDic['raCen'] = self[0].pairs['racen']
        plansDic['decCen'] = self[0].pairs['deccen']

        plansDic['plateID'] = 1000
        if 'plateID' in self.kwargs:
            plansDic['plateID'] = self.kwargs['plateID']

        if 'name' in self.kwargs:
            plansDic['name'] = self.kwargs['name']
        else:
            plansDic['name'] = 'MJ{raCen:.5f}{decCen:+.5f}'.format(**plansDic)

        if 'comment' in self.kwargs:
            plansDic['comment'] = self.kwargs['comment']
        else:
            plansDic['comment'] = ' '

        plansDic['locationID'] = self[0].pairs['locationid']
        plansDic['designID'] = self.designid
        plansDic['plateRun'] = self.plateRun
        plansDic['chunk'] = self.plateRun

        if 'epoch' in self.kwargs:
            plansDic['epoch'] = self.kwargs['epoch']
        else:
            plansDic['epoch'] = '2014.07'

        if 'temp' in self.kwargs:
            plansDic['temp'] = self.kwargs['temp']
        else:
            plansDic['temp'] = '5.0'

        plansDic.update(config['plateTypes'][self.plateType])

        platePlans = platePlansTemplate.format(**plansDic)
        platePlans = platePlans.replace('\n', ' ')

        return platePlans

    def writePlatePlans(self):

        platePlans = self.getPlatePlans()

        filename = 'platePlans-{0}.par'.format(self.designid)

        unitPlatePlans = open(filename, 'w')
        print >>unitPlatePlans, platePlans
        unitPlatePlans.close()

        log.info('platePlans file {0} saved.'.format(filename))

        return filename


class PlateInputBase(object):
    """The base plate input class.

    This class is the base to construct ``PlateInput`` instances. Each
    PlateInputBase object contains the information for one plateInput
    file.

    This class is not intended to be called directly. Instead, PlateInput
    should be used.

    """

    def __init__(self, designid, catalogue, pairs={},
                 reassignFerrules=False, plateRun=None, **kwargs):

        self.designid = designid
        self.inputCatalogue = catalogue
        self._reassignFerrules = reassignFerrules
        self.plateRun = plateRun

        self.filename = None

        self.inputCatalogue.meta.update(pairs)
        self.inputCatalogue.meta.update({'designid': self.designid})

        log.info('Creating PlateInputBase instance for InputCatalogue for ' +
                 'designid={0} and targettype={1}'.format(
                     self.designid, self.meta['targettype']))

        self.checkIFUDesigns()

    def checkIFUDesigns(self, failOnIFUDesignSize=False):

        if self.inputCatalogue.meta['targettype'] == 'science':
            bundleSizes = config['IFUs'].copy()
        elif self.inputCatalogue.meta['targettype'] == 'standard':
            bundleSizes = config['miniBundles'].copy()
        elif self.inputCatalogue.meta['targettype'] == 'sky':
            bundleSizes = config['skies'].copy()

        for size in self.inputCatalogue['ifudesignsize']:
            if size > 0:
                bundleSizes[size] -= 1

        ifuDesignSizeAssigned = False
        for target in self.inputCatalogue:
            if target['ifudesignsize'] < 0:
                for size in bundleSizes:
                    if bundleSizes[size] > 0:
                        target['ifudesignsize'] = size
                        ifuDesignSizeAssigned = True
                        bundleSizes[size] -= 1
                        log.debug(
                            'mangaid={0} has been automatically '.format(
                                target['mangaid']) +
                            'assigned an ifudesignsize={0}'.format(
                                target['ifudesignsize']))

        if ifuDesignSizeAssigned:
            log.info('some ifudesignsizes have been assigned automatically.')

        # self.inputCatalogue.data.pprint(max_width=1000, max_lines=1000)
        # Checks if one or more ifuDesigns are missing.
        missingDesign = False
        for ifuDesign in self.inputCatalogue['ifudesign']:
            if ifuDesign < 0.:
                missingDesign = True
                break

        # If at least one ifuDesign is undefined, reassign ferrules. If
        # reassignFerrules is True, reassign all the ferrules even if they
        # are defined in the input catalogue.
        if missingDesign and not self._reassignFerrules:
            log.info('one or all ifudesigns are missing. Reassigning IFUs.')
            self.inputCatalogue.data = self.assignFerrulesAnchorBlock()
        elif self._reassignFerrules:
            log.info('Reassigning all IFUs.')
            self.inputCatalogue.data = self.assignFerrulesAnchorBlock(
                reassignAll=True)

    def assignFerrulesAnchorBlock(self, reassignAll=False):

        struct1 = self.inputCatalogue.copy()

        # Gets the targets for which IFUs need to be assigned.
        if not reassignAll:
            struct1ToKeep = struct1[struct1['ifudesign'] > 0]
            struct1ToAssign = struct1[struct1['ifudesign'] <= 0]
        else:
            struct1ToKeep = None
            struct1ToAssign = struct1.copy()
            struct1ToAssign['ifudesign'] = -999

        assignOrder, blockPositionsCoo = self._getAssignOrder(struct1ToAssign)

        # Gets already assigned ifudesigns
        usedIFUDesigns = []
        if struct1ToKeep is not None:
            for row in struct1ToKeep:
                usedIFUDesigns.append(int(row['ifudesign']))

        for idx in assignOrder:

            row = struct1ToAssign[idx]

            ra = row['ra']
            dec = row['dec']
            targetCoo = coo.ICRS(ra=ra, dec=dec, unit=(uu.degree, uu.degree))

            ifuSize = row['ifudesignsize']
            ifuDesign = self._assignFerrule(
                targetCoo, ifuSize, blockPositionsCoo, usedIFUDesigns)

            usedIFUDesigns.append(ifuDesign)
            struct1ToAssign[idx]['ifudesign'] = ifuDesign
            log.debug('mangaid={0} assigned ifudesign={1}'.format(
                struct1ToAssign[idx]['mangaid'],
                struct1ToAssign[idx]['ifudesign']))

        if struct1ToKeep is None:
            return struct1ToAssign
        else:
            return table.vstack((struct1ToKeep, struct1ToAssign))

    def _getAssignOrder(self, struct):
        """Calculates the distances of the targets to the anchor blocks.

        Returns the indices of the input structure sorted by distance to
        the anchor blocks (ordered from farther to closer). The metric used
        is the square root of the sum of the squares of the distance of the
        target to each block. The (RA, Dec) position of the blocks is also
        returned in `time` format.

        """

        if len(struct) == 0:
            return struct

        raCen = self.pairs['racen']
        decCen = self.pairs['deccen']

        distances = np.zeros(len(struct), dtype=float)

        blockPositions = self.getRADecAnchorBlocks(raCen, decCen)
        blockPositionsCoo = coo.ICRS(
            ra=list(zip(*blockPositions)[0]),
            dec=list(zip(*blockPositions)[1]),
            unit=(uu.degree, uu.degree))

        for nn, row in enumerate(struct):
            ra = row['ra']
            dec = row['dec']

            targetCoo = coo.ICRS(ra=ra, dec=dec, unit=(uu.degree, uu.degree))
            separations = targetCoo.separation(blockPositionsCoo)

            distances[nn] = np.sum(separations.degree**2)

        return np.argsort(distances)[::-1], blockPositionsCoo

    def _assignFerrule(self, targetCoo, ifuSize,
                       blockPositionsCoo, usedIFUDesigns):
        """
        Returns the optimum ifuDesign for a target based on the target
        position, the desired IFU size and the number of available ifuDesigns.

        """

        ifus = config['IFUs'].copy()
        ifus.update(config['miniBundles'])

        for usedIFUDesign in usedIFUDesigns:
            usedIFUSize = ifuDesign2Size(usedIFUDesign)
            ifus[usedIFUSize] -= 1

        if ifus[ifuSize] <= 0:
            raise ValueError('trying to assign too many '
                             'IFUs of size {0}.'.format(ifuSize))

        anchorBlockOrder = np.argsort(
            targetCoo.separation(blockPositionsCoo).degree)

        sortedBlocks = [config['anchorBlocks']['anchors'][anchorIdx]
                        for anchorIdx in anchorBlockOrder]
        ifuDesign = 0
        for anchorName in sortedBlocks:
            ferrules = config['anchorBlocks'][anchorName]
            for ferrule in ferrules:
                cartMapRow = cartmap[cartmap['frlPlug'] == ferrule][0]
                ferruleSize = cartMapRow['ifusize']
                ferruleDesign = cartMapRow['ifuDesign']
                if ferruleDesign not in usedIFUDesigns and \
                        ferruleSize == ifuSize:
                    ifuDesign = ferruleDesign
                    break
                else:
                    continue
            if ifuDesign != 0:
                break

        return ifuDesign

    @staticmethod
    def getRADecAnchorBlocks(ra, dec):
        """Calculates RA and Dec of the anchor blocks.

        Returns a list of tuples, each tuple containing RA and Dec.
        The order of blocks in the returned list is [NE, SE, NW, SE]

        """

        raAngToAnchor = 1.5 / 2. / np.cos(dec * np.pi / 180.)
        decAngToAnchor = 1.5 / 2

        return [
            (ra + raAngToAnchor, dec + decAngToAnchor),
            (ra + raAngToAnchor, dec - decAngToAnchor),
            (ra - raAngToAnchor, dec + decAngToAnchor),
            (ra - raAngToAnchor, dec - decAngToAnchor)]

    def plotIFUs(self, filename, struct1=None,
                 copyToRepo=False, moveToRepo=False):
        """Plots the position of the IFUs on the plate.

        This method plots the position of the IFUs (and, thus,
        of the targets) on the plate with information about
        MaNGAID and IFUDESIGN.

        Parameters
        ----------
        filename : str
            The plot filename.
        struct1 : `Table` or None, optional
            If specified, uses the data from the structure to
            determine the MaNGAIDs, RA, Dec and IFUDESIGN. Otherwise,
            uses the one in the current instance.
        copyToRepo, moveToRepo : bool, optional
            If True, copies or moves the plot to the platelist plate run
            directory. If both are True, copyToRepo supersedes moveToRepo.

        """

        if struct1 is None:
            struct1 = self.inputCatalogue.data

        from matplotlib import pyplot as plt
        from matplotlib.patches import Ellipse

        plt.cla()
        plt.clf()

        fig, ax = plt.subplots()
        fig.set_size_inches(8, 8)

        raCen = self.pairs['racen']
        decCen = self.pairs['deccen']

        plate = Ellipse((raCen, decCen),
                        height=3.,
                        width=3/np.cos(decCen*np.pi/180.),
                        linewidth=2,
                        edgecolor='k', facecolor='None')
        ax.add_patch(plate)

        for row in struct1:

            ra = row['ra']
            dec = row['dec']
            ifuSize = row['ifudesignsize']
            ifuDesign = row['ifudesign']
            mangaID = row['mangaid']

            cartMapInfo = cartmap[cartmap['ifuDesign'] == ifuDesign][0]
            ferrule = cartMapInfo['frlPlug']

            anchorName = ''
            for key in config['anchorBlocks']['anchors']:
                if ferrule in config['anchorBlocks'][key]:
                    anchorName = key
                    break

            closestAnchor = self._isClosestAnchor(
                ra, dec, raCen, decCen, anchorName)

            if closestAnchor:
                color = 'b'
            else:
                color = 'r'

            ax.scatter(ra, dec, s=ifuSize/2, marker='o',
                       edgecolor=color, facecolor='None')

            ax.text(ra, dec-0.06,
                    '{0} ({1})'.format(ifuDesign, anchorName),
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=8, color=color)

            ax.text(ra, dec+0.05, r'{0}'.format(mangaID),
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=6, color=color)

        ax.set_xlim(raCen + 1.6/np.cos(decCen*np.pi/180.),
                    raCen - 1.6/np.cos(decCen*np.pi/180.))
        ax.set_ylim(decCen - 1.6, decCen + 1.6)

        ax.set_xlabel(r'$\alpha_{2000}$')
        ax.set_ylabel(r'$\delta_{2000}$')

        plt.savefig(filename)

        plt.close('all')

        log.info('IFU plot for designid={0} saved to {1}'.format(
            self.designid, filename))

        if copyToRepo:
            self._toPlateListInputs(filename)
        elif moveToRepo:
            self._toPlateListInputs(filename, command=sh.move)

    def _isClosestAnchor(self, ra, dec, raCen, decCen, anchorName):
        """Returns True if an ifuDesign is from the closest anchor."""

        blockPositions = self.getRADecAnchorBlocks(raCen, decCen)
        blockPositionsCoo = coo.ICRS(
            ra=list(zip(*blockPositions)[0]),
            dec=list(zip(*blockPositions)[1]),
            unit=(uu.degree, uu.degree))

        targetCoo = coo.ICRS(
            ra=ra, dec=dec,
            unit=(uu.degree, uu.degree))

        closestAnchorIdx = np.argmin(
            targetCoo.separation(blockPositionsCoo).degree)
        closestAnchor = config['anchorBlocks']['anchors'][closestAnchorIdx]

        if closestAnchor == anchorName:
            return True
        else:
            return False

    def write(self, filename=None, toRepo=False):

        if filename is None:
            if self.filename is not None:
                filename = self.filename
            else:
                filename = self.getDefaultFilename()
                self.filename = filename

        if os.path.exists(filename):
            os.remove(filename)

        yanny.write_ndarray_to_yanny(filename, self.inputCatalogue.data,
                                     structname='MANGAINPUT',
                                     hdr=self.inputCatalogue.meta)

        log.info('{0} saved.'.format(filename))

        if toRepo:
            self._toRepo(filename)

        return filename

    def getDefaultFilename(self):

        try:
            field = '{0:04d}'.format(self.inputCatalogue.meta['fieldname'])
        except:
            field = self.inputCatalogue.meta['fieldname']

        template = 'manga{type}_{field:s}_{designid:04d}.par'.format(
            type=self.inputCatalogue.meta['targettype'].title(),
            field=field, designid=self.designid)

        return template

    def _toRepo(self, filename):
        """Copies a file to platelist."""

        inputPath = os.path.join(
            os.path.expandvars(config['platelist']), 'inputs/manga',
            self.plateRun)

        if not os.path.exists(inputPath):
            os.makedirs(inputPath)

        sh.copy(filename, inputPath)

        log.info(filename + ' copied $PLATELIST/inputs.')

    @property
    def struct1(self):
        return self.inputCatalogue.data

    @property
    def meta(self):
        return self.inputCatalogue.meta

    @property
    def pairs(self):
        return self.inputCatalogue.meta
