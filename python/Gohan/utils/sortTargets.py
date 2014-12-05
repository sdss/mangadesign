#!/usr/bin/env python
# encoding: utf-8
"""
sortTargets.py

Created by José Sánchez-Gallego on 6 Nov 2014.
Licensed under a 3-clause BSD license.

Revision history:
    6 Nov 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse


def simpleMesh(centre, radius, width=0.2, **kwargs):
    """Creates a simple mesh for the field. Returns the centres of the cells"""

    ra0 = centre[0] - radius / np.cos(np.deg2rad(centre[1]))
    ra1 = centre[0] + radius / np.cos(np.deg2rad(centre[1]))
    dec0 = centre[1] - radius
    dec1 = centre[1] + radius

    mRA, mDec = np.meshgrid(np.arange(ra0, ra1+width, width),
                            np.arange(dec0, dec1+width, width))

    coords = np.array([mRA.flatten(), mDec.flatten()]).T
    distanceToCentre = calculateSeparation(coords, centre)

    return coords[distanceToCentre <= radius]


def plotTargets(targets, centre, plotGrid=None, plotAllTargets=None,
                filename=None, **kwargs):
    """Creates a simple plot with the targets."""

    raCen, decCen = centre

    if filename is None:
        filename = 'sortedTargets.pdf'

    plt.clf()
    plt.cla()

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)

    plate = Ellipse((raCen, decCen),
                    height=3.,
                    width=3/np.cos(decCen*np.pi/180.),
                    linewidth=2,
                    edgecolor='k', facecolor='None')
    ax.add_patch(plate)

    ax.scatter(centre[0], centre[1], marker='o', s=25, color='k',
               edgecolor='k')

    if plotGrid is not None:
        ax.scatter(plotGrid[:, 0], plotGrid[:, 1], marker='x',
                   color='0.5', s=10.)

    if plotAllTargets is not None:
        ax.scatter(plotAllTargets[:, 0], plotAllTargets[:, 1],
                   marker='x', s=20, color='0.5')

    ax.scatter(targets[:, 0], targets[:, 1], marker='x', s=20, color='r')

    ax.set_xlim(raCen + 1.6/np.cos(decCen*np.pi/180.),
                raCen - 1.6/np.cos(decCen*np.pi/180.))
    ax.set_ylim(decCen - 1.6, decCen + 1.6)

    ax.set_xlabel(r'$\alpha_{2000}$')
    ax.set_ylabel(r'$\delta_{2000}$')

    plt.savefig(filename)

    plt.close('all')

    return


def calculateSeparation(coord1, coord2):
    """Calculates the separation in the sky between a list of targets. This
    function replaces `astropy.coordinates.SkyCoord.separation`, which is
    still very slow for a large number of points."""

    coord1 = np.atleast_2d(coord1)
    coord2 = np.atleast_2d(coord2)

    assert coord1.shape[0] == 1 or coord2.shape[0] == 1

    if coord1.shape[0] > 1:
        coord2, coord1 = coord1, coord2

    return np.sqrt(((coord1[:, 0] - coord2[:, 0]) *
                    np.cos(np.deg2rad(coord1[:, 1])))**2
                   + (coord1[:, 1]-coord2[:, 1])**2)


def getTargetIdx(coords, assigned, grid, centre):
    """Returns the index of the target to assign, based on what has already
    been assigned."""

    assignedCoords = coords[assigned]
    eField = np.zeros(len(grid))

    for ii in range(len(assignedCoords)):
        distance = calculateSeparation(assignedCoords[ii], grid)
        eField += 1. / distance

    eField += .1 / calculateSeparation(centre, grid)

    for alpha in range(0, 360, 10):
        ra = centre[0] + (1.49 * np.cos(np.deg2rad(alpha)) /
                          np.cos(np.deg2rad(centre[1])))
        dec = centre[1] + 1.49 * np.sin(np.deg2rad(alpha))
        eField += .5 / calculateSeparation(np.array([(ra, dec)]), grid)

    return np.argmin(eField)


def sortTargets(targets, centre=None, radius=1.49,
                limitTo=None, plot=False, **kwargs):
    """Sorts a list of targets, evenly distributing them in a plate.

    This routine takes a list of targets and distributes them as uniformly as
    in a circular field of radius `radius`. Targets are decollided against
    themselves. The routine uses an electrostatic approach, replacing targets
    with positive charges and calculating the electric field in a grid.

    Parameters
    ----------
    targets : numpy.ndarray
        A Numpy array of shape NxM where N>=1 and M>2. Each row must be a
        target. The first two columns in the array will be considered the RA
        and Dec for the targets. Any other column will be ignored.

    centre : tuple-like object, opional
        A tuple or array of shape 1x2 with the coordinates of the field centre.
        If None, the centre is calculated as the centre of mass of all the
        targets.

    radius : float, optional
        The radius of the field, in degrees

    limitTo : int or None, optional
        If set to an integer, returns only that number of sorted targets.

    plot : bool, optional
        If True, a plot is generated. See `plotTargets` for more information
        on which kwargs parameters can be passed to the function.

    Returns
    -------
    result : tubple
        A tuple containing, first, the same input `targetList` but reordered
        according to the sorting algorithm. If `limitTo` is not None,
        the length of the list is the same as `limitTo`. The second element is
        index order referred to the original array.

    """

    targets = np.atleast_2d(targets)

    if limitTo is None:
        limitTo = targets.shape[0]

    if centre is None:
        centre = np.mean(targets, axis=0)

    assert limitTo <= targets.shape[0]

    grid = simpleMesh(centre, radius, **kwargs)

    assigned = []

    while len(assigned) < limitTo:

        newGridIdx = getTargetIdx(targets, assigned, grid, centre)

        distancesToGrid = calculateSeparation(targets, grid[newGridIdx])
        for ii in np.argsort(distancesToGrid):
            if ii not in assigned:
                assigned.append(ii)
                break

    sortedTargets = np.array([targets[ii] for ii in assigned])

    if plot:
        plotTargets(sortedTargets, centre, plotAllTargets=targets,
                    plotGrid=None, **kwargs)

    return sortedTargets, assigned
