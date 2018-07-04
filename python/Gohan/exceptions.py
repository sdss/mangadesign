#!/usr/bin/env python
# encoding: utf-8
"""
exceptions.py

Created by José Sánchez-Gallego on 28 Feb 2014.
Licensed under a 3-clause BSD license.

Revision history:
    28 Feb 2014 J. Sánchez-Gallego
      Initial version

"""

from __future__ import absolute_import, division, print_function, unicode_literals


class GohanError(Exception):
    """Base exception for Gohan. Other exceptions should inherit this."""
    pass


class GohanPostDesignError(Exception):
    """An exception for MaNGA postdesign."""
    pass


class GohanNotImplemented(GohanError):
    """A class for exceptions about functionalities not yet implemented."""
    pass


class GohanPlateInputError(GohanError):
    """A class for error related with PlateInput."""
    pass


class GohanPlateTargetsError(Exception):
    """An exception for PlateTargets."""
    pass


class GohanWarning(Warning):
    """Base warning for Gohan."""


class GohanUserWarning(UserWarning, GohanWarning):
    """The primary warning class."""


class GohanPostDesignWarning(UserWarning, GohanWarning):
    """Warning for the postdesign process."""


class GohanPlateTargetsWarning(UserWarning, GohanWarning):
    """Warning for PlateTargets."""


class GohanCollisionWarning(UserWarning, GohanWarning):
    """A warning raised when a target collides."""
