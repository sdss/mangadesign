#!/usr/bin/env python
# encoding: utf-8
"""
configuration.py

Created by José Sánchez-Gallego on 10 Dec 2013.
Copyright (c) 2013. All rights reserved.
Licensed under a 3-clause BSD license.

"""

import yaml

try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib


def merge(user, default):
    """Merges a user configuration with the default one."""

    if not user:
        return default

    if isinstance(user, dict) and isinstance(default, dict):
        for kk, vv in default.items():
            if kk not in user:
                user[kk] = vv
            else:
                user[kk] = merge(user[kk], vv)

    return user


def get_config(user_path):
    """Returns a dictionary object with configuration options."""

    user_path = pathlib.Path(user_path).expanduser()
    user = user_path.exists() and yaml.load(open(str(user_path), 'r'), Loader=yaml.FullLoader)

    default_path = pathlib.Path(__file__).parents[0] / '../defaults.yaml'
    default = yaml.load(open(str(default_path), 'r'), Loader=yaml.FullLoader)

    return merge(user, default)
