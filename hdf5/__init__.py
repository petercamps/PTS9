#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.hdf5 Support for SKIRT input and output files stored as bundles in an HDF5 file
#
# SKIRT can optionally read its input files from, and write its output files to, an HDF5 file. Each input or
# output file then corresponds to a \em bundle: an HDF5 group that holds one or more datasets with names and
# attributes defined by SKIRT. Bundles can be organized in a \em suite, i.e. a (possibly nested) group, so that
# a single HDF5 file can hold the input and output of many simulations.
#
# This package contains functions to write each type of bundle from data in memory and to read it into memory
# (see the pts.hdf5.bundle module), and to convert the regular SKIRT input and output file types to the
# corresponding bundle and back (see the pts.hdf5.convert module). The functions in this package require the
# \c h5py package to be installed.
#

from .bundle import splitBundlePath, joinBundlePath, resolveBundlePath, datasetNames, \
                    writeColumnBundle, writeAdaptiveMeshBundle, writeStoredTableBundle, \
                    writeFitsInputBundle, writeFitsOutputBundle, writeGridPlotBundle, writeTextBundle, \
                    bundleTypes, inferBundleType, findBundles, \
                    readColumnBundle, readAdaptiveMeshBundle, readStoredTableBundle, \
                    readFitsBundle, readGridPlotBundle, readTextBundle
from .convert import inferFileType, fileTypes, convertFile, convertFiles, \
                     convertTextColumnFile, convertAdaptiveMeshFile, convertStoredColumnsFile, \
                     convertStoredTableFile, convertFitsInputFile, convertFitsOutputFile, \
                     convertGridPlotFile, convertTextFile, \
                     convertBundle, convertBundles, filesToBeReplaced, \
                     convertTextColumnBundle, convertAdaptiveMeshBundle, convertStoredTableBundle, \
                     convertFitsBundle, convertGridPlotBundle, convertTextBundle
