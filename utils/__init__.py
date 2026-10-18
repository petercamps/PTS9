#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.utils General utilities
#
# This package includes utilities deemed sufficiently generic that they don’t belong in one of the other packages,
# and that are likely shared between multiple packages.
#

from .config import setInteractive, interactive
from .error import UserError
from .path import dataPath, ptsPath, ptsResourcesPath, projectParentPath, skirtPath, skirtResourcesPath, \
                  absPath, savePath, BundlePath, isBundlePath
from .envir import timestamp, username, hostname
