#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.hdf5.bundle Reading and writing SKIRT bundles in an HDF5 file
#
# The functions in this module write data held in memory to an HDF5 file, in the form expected by SKIRT when it
# reads its input from an HDF5 file, and produced by SKIRT when it writes its output to an HDF5 file; and they read
# such data from an HDF5 file into memory. They can be used, for example, by a script that prepares the input for a
# large number of SKIRT simulations, or by a script that processes the output of these simulations.
#
# Bundles and suites
# ------------------
#
# The contents of a SKIRT input or output file are represented in an HDF5 file by a \em bundle: an HDF5 group
# holding one or more datasets, with dataset names and attributes defined by SKIRT for each type of file. There is a
# function in this module to write and to read each type of bundle:
#
# | Bundle type | Writing | Reading | Corresponding regular file
# |-------------|---------|---------|---------------------------
# | Text column file bundle | writeColumnBundle() | readColumnBundle() | Column text file (input and output)
# | AMR text file bundle | writeAdaptiveMeshBundle() | readAdaptiveMeshBundle() | Adaptive mesh text file (input)
# | Stored table bundle | writeStoredTableBundle() | readStoredTableBundle() | SKIRT stored table file (input)
# | FITS file bundle for input | writeFitsInputBundle() | readFitsBundle() | FITS file with just pixel values (input)
# | FITS file bundle for output | writeFitsOutputBundle() | readFitsBundle() | FITS file with data and axes (output)
# | Spatial grid plot file bundle | writeGridPlotBundle() | readGridPlotBundle() | Spatial grid plot file (output)
# | Unstructured text file bundle | writeTextBundle() | readTextBundle() | Log, parameter, or free-form text file
#
# A bundle resides in the root group of the HDF5 file or in a \em suite, i.e. an HDF5 group that holds bundles
# and that itself resides in the root group or in another group. For example, the input for each of the simulations
# in a campaign can be stored in a separate suite of the same HDF5 file.
#
# Bundle paths
# ------------
#
# The functions in this module identify a bundle by its <em>bundle path</em>, a string of the form
#
#     <filepath>:<suite>/<bundle>
#
# where \c \<filepath\> is the absolute or relative path of the HDF5 file, including the mandatory \c .hdf5 filename
# extension; \c \<suite\> is the optional path of the suite within the file, with nested groups separated by forward
# slashes; and \c \<bundle\> is the name of the bundle. The file path is interpreted as described for the
# pts.utils.absPath() function. For example:
#
#     in/data.hdf5:stars.txt
#     ~/campaign/data.hdf5:campaign7/galaxy042/stars.txt
#
# The name of a bundle usually is the name of the regular file it replaces, including the filename extension,
# because that is the name under which SKIRT looks for the bundle.
#
# Creating and replacing bundles
# ------------------------------
#
# The HDF5 file as a whole is never replaced. If the file does not exist, it is created. If it exists, the bundle
# is added to it. If the file already holds a bundle with the same path, that bundle is replaced; all other bundles
# in the file are left untouched. Suites are created as needed. Note that HDF5 does not release the space taken up
# by a replaced bundle; use the \c h5repack tool provided with the HDF5 library to compact a file if needed.
#
# Every bundle written by this module carries three attributes that describe its origin: \c producer (a string
# identifying PTS), \c format (the version of the bundle format, currently 1), and \c created (the date and time
# at which the bundle was written, in ISO 8601 format with milliseconds).
#
# Dataset names
# -------------
#
# In a text column file bundle, an AMR text file bundle, and a stored table bundle, the datasets are named for
# the columns, axes, or quantities they represent. Because a column in a column text file is identified by a
# free-form description, and several columns can have the same description, these names are adjusted where needed
# so that they are valid and unique within the bundle. See the datasetNames() function for the rules. The reading
# functions undo these adjustments, so that they return the original descriptions.
#
# Reading bundles
# ---------------
#
# Each of the reading functions returns a dictionary with the contents of a bundle. The keys of this dictionary are
# the names of the arguments of the corresponding writing function, other than the bundle path and the compression
# level. For example,
#
#     table = readColumnBundle("in/data.hdf5:stars.txt")
#
# returns a dictionary with the keys \c columns, \c names, \c units, and \c description, and
#
#     writeColumnBundle("in/other.hdf5:stars.txt", **table)
#
# writes a bundle with the same contents to another HDF5 file. Just like SKIRT, the reading functions use the
# datasets and attributes defined for the type of bundle and ignore any other ones, and they convert numbers to
# 64-bit floating point numbers, regardless of the type in which they are stored. Compressed datasets are read
# transparently.
#
# A bundle does not record its type. The inferBundleType() function derives the type of a bundle from the datasets
# and attributes it holds, and returns one of the following names:
#
# | Type | Bundle type | Reading function
# |------|-------------|-----------------
# | columns | Text column file bundle | readColumnBundle()
# | amr | AMR text file bundle | readAdaptiveMeshBundle()
# | storedtable | Stored table bundle | readStoredTableBundle()
# | fits | FITS file bundle for input or for output | readFitsBundle()
# | gridplot | Spatial grid plot file bundle | readGridPlotBundle()
# | text | Unstructured text file bundle | readTextBundle()
#
# The findBundles() function lists the bundles in a suite, or those with a name matching a wildcard pattern.
#
# Data types and compression
# --------------------------
#
# All floating point data is stored as 64-bit floating point numbers, integer data as 32-bit integers, and Boolean
# data in the form used by \c h5py for Booleans (an enumeration stored as an 8-bit integer). The data passed
# to the functions is converted to these types as needed.
#
# The functions that write potentially large datasets have an optional \em compression argument. If this argument
# is omitted, None, or zero, the datasets are stored without compression, as a single contiguous block. If it is
# an integer from 1 (fastest) to 9 (smallest), the datasets are stored in chunks that are compressed with the
# \c gzip filter at that level, after reordering the bytes with the shuffle filter. SKIRT reads compressed datasets
# transparently.
#
# \note These low-level functions do \em not support astropy quantities; they expect and return plain numpy arrays
#       (or sequences) with values in the units declared by the corresponding unit strings.
#

# -----------------------------------------------------------------

import contextlib
import datetime
import fnmatch
import logging
import re
import h5py
import numpy as np
import pts.utils as ut

# -----------------------------------------------------------------

# the filename extension that identifies an HDF5 file in a bundle path or in an output location
_extension = ".hdf5"

# the value of the "producer" attribute for the bundles written by this module
_producer = "PTS (Python Toolkit for SKIRT)"

# the value of the "format" attribute for the bundles written by this module
_format = 1

# the values allowed for the "type" attribute of an unstructured text file bundle
_textTypes = ("free form", "XML", "log")

# the names of the attributes that describe the observer in a FITS file bundle for output:
# those with a numeric value, and the one with a string value
_observerNumbers = ("inclination", "azimuth", "roll", "redshift", "luminosity_distance", "angular_diameter_distance")
_observerStrings = ("distance_unit",)

# the names of the bundle types returned by the inferBundleType() function, in the order of the table in the header
_bundleTypes = ("columns", "amr", "storedtable", "fits", "gridplot", "text")

# the syntax of a dataset name with the sequence number added by the datasetNames() function to make it unique
_sequenceNumberSyntax = re.compile(r"(.+) \((\d+)\)")

# the syntax of a name with wildcards, as supported by the standard fnmatch module
_wildcardSyntax = re.compile(r"[*?[]")

# -----------------------------------------------------------------

## This function splits the specified bundle path, of the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt>, and
# returns a tuple with two strings: the path of the HDF5 file (as specified, without converting it to an absolute
# path), and the path within that file, without leading or trailing slashes. The path within the file is empty if the
# argument consists of just the path of an HDF5 file. The path is split at the first occurrence of the \c .hdf5
# filename extension followed by a colon, in the same way as SKIRT does. The function raises a user error if the
# argument does not refer to an HDF5 file.
def splitBundlePath(bundlePath):
    bundlePath = str(bundlePath)
    index = bundlePath.find(_extension + ":")
    if index >= 0:
        index += len(_extension)
        return bundlePath[:index], _cleanGroupPath(bundlePath[index+1:])
    if bundlePath.endswith(_extension):
        return bundlePath, ""
    raise ut.UserError("Path does not refer to an HDF5 file with the {} filename extension: '{}'" \
                       .format(_extension, bundlePath))

## This function returns the bundle path, of the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt>, for the HDF5 file
# with the specified path and the group within that file obtained by joining the remaining arguments with forward
# slashes. Empty arguments are ignored, so that the function can be called with an empty suite. For example:
#
#     joinBundlePath("in/data.hdf5", "campaign7/galaxy042", "stars.txt")
#
# returns the string "in/data.hdf5:campaign7/galaxy042/stars.txt".
def joinBundlePath(filePath, *groupNames):
    return "{}:{}".format(filePath, _cleanGroupPath("/".join([ str(name) for name in groupNames ])))

## This function returns the bundle path for a bundle with the specified name at the specified output location.
# It gives the functions and commands that convert files to bundles the same semantics for their output location
# as the \c -o command line option of SKIRT.
#
# The \em outPath argument specifies the output location in one of the following forms:
#
# | Form | Meaning
# |------|--------
# | None or empty string | The HDF5 file specified by the \em defaultFilePath argument
# | <tt>\<dir\></tt> | The HDF5 file with the name of \em defaultFilePath, in the specified existing directory
# | <tt>\<dir\>/\<hdf\></tt> | The specified HDF5 file, with a name ending in \c .hdf5, in the specified directory
# | <tt>\<dir\>/\<hdf\>:\<suite\></tt> | The specified suite in the specified HDF5 file
#
# In all cases, the directory that holds the HDF5 file must exist; the HDF5 file itself is created when the first
# bundle is written to it. The paths are interpreted as described for the pts.utils.absPath() function. The
# \em bundleName argument can include a suite, separated from the name by a forward slash, which is then located
# within the suite of the output location.
def resolveBundlePath(outPath, bundleName, defaultFilePath):
    outPath = "" if outPath is None else str(outPath)
    suite = ""
    if not outPath:
        filePath = ut.absPath(defaultFilePath)
    elif ut.absPath(outPath).is_dir():
        filePath = ut.absPath(outPath) / ut.absPath(defaultFilePath).name
    else:
        try:
            filePath, suite = splitBundlePath(outPath)
        except ut.UserError:
            raise ut.UserError("Output location is not an existing directory, and does not specify an HDF5 file " \
                               "with the {} filename extension: '{}'".format(_extension, outPath)) from None
        filePath = ut.absPath(filePath)

    if not filePath.parent.is_dir():
        raise ut.UserError("Output directory does not exist: '{}'".format(filePath.parent))
    return joinBundlePath(filePath, suite, bundleName)

# -----------------------------------------------------------------

## This function converts a sequence of column descriptions, or of other names, to a list with the same number of
# dataset names that are valid in an HDF5 file and that differ from each other. It defines the names of the
# datasets in a text column file bundle, an AMR text file bundle, and a stored table bundle. The following rules
# are applied to each description, in this order:
#
#  - White space at the start and at the end is removed, and consecutive white space characters are replaced
#    by a single space. This is the form in which SKIRT compares column descriptions.
#  - A forward slash is replaced by a vertical bar, because HDF5 uses the forward slash to separate the names
#    of nested groups. For example, the description "n_He/n_H" becomes the dataset name "n_He|n_H".
#  - An empty description is replaced by the word "column" followed by the zero-based index of the description
#    in the sequence, as in "column3".
#  - If the resulting name already occurs earlier in the sequence, or is one of the names listed in the optional
#    \em reserved argument, a sequence number between parentheses is added to it. For example, three columns
#    described as "weight" result in the dataset names "weight", "weight (2)", and "weight (3)".
#
# A column description in a column text file cannot contain parentheses, and SKIRT does not use vertical bars in
# the descriptions it writes, so that the original description can be recognized in the dataset name.
def datasetNames(descriptions, reserved=()):
    names = []
    used = set(reserved)
    for index, description in enumerate(descriptions):
        name = " ".join(str(description).split()).replace("/", "|")
        if name in ("", ".", ".."):
            name = "column{}".format(index)
        unique = name
        count = 1
        while unique in used:
            count += 1
            unique = "{} ({})".format(name, count)
        used.add(unique)
        names.append(unique)
    return names

# -----------------------------------------------------------------

## This function writes a text column file bundle, the equivalent of a SKIRT column text file: a table with a row
# for each item (for example a particle, a cell, or a wavelength) and a column for each property.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em columns: a sequence of one-dimensional arrays of the same length, one for each column; or a single
#    two-dimensional array with shape (numColumns, numRows). There must be at least one column; the number of rows
#    can be zero.
#  - \em names: a sequence of strings with the description of each column. The datasets are named for these
#    descriptions, adjusted as described for the datasetNames() function.
#  - \em units: a sequence of SKIRT unit strings, one for each column, or omitted or None if the units are not
#    specified. An individual item in the sequence can be None as well. SKIRT uses its default unit for a column
#    without a unit; a dimensionless column has the unit string "1".
#  - \em description: an optional free-form string that describes the contents of the bundle.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds a one-dimensional dataset of 64-bit floating point numbers for each column, with two attributes:
# \c column, the zero-based index of the column, and \c unit, the unit string (if specified). Note that, in contrast,
# the column numbers in the header of a column text file are one-based. SKIRT matches the datasets to the columns
# it expects in the same way as for a column text file: by default in the order of the columns, which is given
# by the \c column attribute, and by name if the ski file specifies the \em useColumns property.
def writeColumnBundle(bundlePath, columns, names, units=None, *, description=None, compression=None):
    columns, names, units = _verifiedColumns(columns, names, units)
    compression = _verifiedCompression(compression)

    with _newBundle(bundlePath, "text column file") as bundle:
        _setAttributes(bundle, description=description)
        for index, (column, name, unit) in enumerate(zip(columns, datasetNames(names), units)):
            _addDataset(bundle, name, column, np.float64, compression, column=index, unit=unit)

# -----------------------------------------------------------------

## This function writes an AMR text file bundle, the equivalent of a SKIRT column text file that describes an
# adaptive mesh: a hierarchical tree of nodes in which each nonleaf node is subdivided into a regular grid
# of child nodes. The nodes are listed in the same order as the lines in the text file: the root node first,
# and the children of a nonleaf node immediately after that node, recursively, with the first index of the
# grid of children varying most rapidly.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em isLeaf: a one-dimensional array of Booleans with an item for each node, in the order described above:
#    True for a leaf node (a cell) and False for a nonleaf node.
#  - \em topology: a two-dimensional array of integers with shape (numNonLeafNodes, 3), listing the number of
#    subdivisions in the x, y, and z direction for each nonleaf node, in the order in which these nodes occur.
#  - \em columns, \em names, \em units: the properties of the leaf nodes, as described for the writeColumnBundle()
#    function, with a row for each leaf node, in the order in which these nodes occur.
#  - \em description: an optional free-form string that describes the contents of the bundle.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds the datasets \c is_leaf (Booleans) and \c topology (32-bit integers), and a dataset for each
# column, as in a text column file bundle. A column description equal to one of these two names receives a
# sequence number, as described for the datasetNames() function.
def writeAdaptiveMeshBundle(bundlePath, isLeaf, topology, columns, names, units=None, *,
                            description=None, compression=None):
    columns, names, units = _verifiedColumns(columns, names, units)
    compression = _verifiedCompression(compression)

    # verify the structure of the tree
    isLeaf = np.asarray(isLeaf)
    if isLeaf.dtype != bool or isLeaf.ndim != 1:
        raise ut.UserError("The leaf flags of an adaptive mesh must be a one-dimensional array of Booleans")
    topology = np.asarray(topology)
    if topology.size == 0: topology = topology.reshape(0, 3)
    if not np.issubdtype(topology.dtype, np.integer) or topology.ndim != 2 or topology.shape[1] != 3:
        raise ut.UserError("The topology of an adaptive mesh must be an array of integers with shape (N, 3)")
    if topology.size > 0 and (topology.min() < 1 or topology.max() > np.iinfo(np.int32).max):
        raise ut.UserError("The subdivision counts in the topology of an adaptive mesh must be positive " \
                           "32-bit integers")
    numLeaves = int(np.count_nonzero(isLeaf))
    if len(topology) != len(isLeaf) - numLeaves:
        raise ut.UserError("The topology of an adaptive mesh has {} rows for {} nonleaf nodes" \
                           .format(len(topology), len(isLeaf) - numLeaves))
    if len(columns[0]) != numLeaves:
        raise ut.UserError("The columns of an adaptive mesh have {} rows for {} leaf nodes" \
                           .format(len(columns[0]), numLeaves))

    # verify that the nodes form a single complete tree: each node fills one open position, and a nonleaf node
    # opens a position for each of its children; the tree is complete when no open positions remain
    opened = np.full(len(isLeaf), -1, dtype=np.int64)
    opened[~isLeaf] += np.prod(topology.astype(np.int64), axis=1)
    remaining = 1 + np.cumsum(opened)
    if len(remaining) == 0 or remaining[-1] > 0:
        raise ut.UserError("The adaptive mesh ends before all nodes of the tree are specified")
    if np.any(remaining[:-1] <= 0):
        raise ut.UserError("The adaptive mesh has superfluous nodes after all nodes of the tree are specified")

    with _newBundle(bundlePath, "AMR text file") as bundle:
        _setAttributes(bundle, description=description)
        _addDataset(bundle, "is_leaf", isLeaf, bool, compression)
        _addDataset(bundle, "topology", topology, np.int32, compression)
        datasets = datasetNames(names, reserved=("is_leaf", "topology"))
        for index, (column, name, unit) in enumerate(zip(columns, datasets, units)):
            _addDataset(bundle, name, column, np.float64, compression, column=index, unit=unit)

# -----------------------------------------------------------------

## This function writes a stored table bundle, the equivalent of a file in SKIRT stored table format: one or more
# quantities tabulated on the grid formed by one or more axes. Its arguments, other than the bundle path and the
# compression level, are the same as those of the pts.storedtable.io.writeStoredTable() function:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em axisNames: a sequence of axis name strings.
#  - \em axisUnits: a sequence of axis unit strings.
#  - \em axisScales: a sequence of axis scale strings ('lin' or 'log') controlling interpolation behavior.
#  - \em axisGrids: a sequence of one-dimensional arrays specifying the grid points for each axis.
#  - \em quantityNames: a sequence of quantity name strings.
#  - \em quantityUnits: a sequence of quantity unit strings.
#  - \em quantityScales: a sequence of quantity scale strings ('lin' or 'log') controlling interpolation behavior.
#  - \em quantityValues: a sequence of arrays specifying the values for each quantity; each array has indices in
#    the same order and with the same range as the specified axes.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds a one-dimensional dataset for each axis, with the attributes \c axis (the zero-based index of
# the axis), \c unit, and \c log (a Boolean that is True for a logarithmic scale); and a dataset for each quantity,
# with a dimension for each axis in the order of the axes, and with the attributes \c unit and \c log. The datasets
# are named for the axes and the quantities, adjusted as described for the datasetNames() function if needed.
#
# \note SKIRT locates each axis and each quantity by the name of its dataset. The names of the axes and quantities
# of a table should therefore all differ from each other, which a file in SKIRT stored table format does not
# require. If a quantity has the name of an axis, or if two axes or two quantities have the same name, the dataset
# for the second one receives a sequence number, and SKIRT cannot read that axis or quantity from the bundle. The
# readStoredTableBundle() function does read such a bundle. None of the stored tables in the SKIRT resources has
# axes or quantities with the same name.
def writeStoredTableBundle(bundlePath, axisNames, axisUnits, axisScales, axisGrids,
                           quantityNames, quantityUnits, quantityScales, quantityValues, *, compression=None):
    compression = _verifiedCompression(compression)

    # verify the specified data
    numAxes = len(axisNames)
    if numAxes < 1 or numAxes != len(axisUnits) or numAxes != len(axisScales) or numAxes != len(axisGrids):
        raise ut.UserError("Mismatch in number of axes for stored table")
    numQuantities = len(quantityNames)
    if numQuantities < 1 or numQuantities != len(quantityUnits) or numQuantities != len(quantityScales) \
                         or numQuantities != len(quantityValues):
        raise ut.UserError("Mismatch in number of quantities for stored table")
    axisGrids = [ np.asarray(grid, dtype=np.float64) for grid in axisGrids ]
    if any([ grid.ndim != 1 for grid in axisGrids ]):
        raise ut.UserError("The grid for each axis of a stored table must be a one-dimensional array")
    quantityValues = [ np.asarray(values, dtype=np.float64) for values in quantityValues ]
    shape = tuple([ len(grid) for grid in axisGrids ])
    if any([ values.shape != shape for values in quantityValues ]):
        raise ut.UserError("Mismatch in number of values for stored table")
    axisLogs = [ _isLogScale(scale) for scale in axisScales ]
    quantityLogs = [ _isLogScale(scale) for scale in quantityScales ]

    with _newBundle(bundlePath, "stored table") as bundle:
        names = datasetNames(list(axisNames) + list(quantityNames))
        for index in range(numAxes):
            _addDataset(bundle, names[index], axisGrids[index], np.float64, compression,
                        axis=index, unit=axisUnits[index], log=axisLogs[index])
        for index in range(numQuantities):
            _addDataset(bundle, names[numAxes+index], quantityValues[index], np.float64, compression,
                        unit=quantityUnits[index], log=quantityLogs[index])

# -----------------------------------------------------------------

## This function writes a FITS file bundle for input, the equivalent of a FITS file with a two-dimensional image
# or a three-dimensional data cube from which SKIRT reads just the pixel values, for example to define a geometry.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em data: a two-dimensional array with shape (ny, nx) or a three-dimensional array with shape (nz, ny, nx).
#    This is the index order of the array obtained when reading a FITS file with the \c astropy.io.fits package.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds a single dataset named \c data, without attributes. This is also the name of the dataset that
# holds the frame or data cube in a FITS file bundle for output, see the writeFitsOutputBundle() function. SKIRT
# reads just this dataset from either type of bundle.
def writeFitsInputBundle(bundlePath, data, *, compression=None):
    compression = _verifiedCompression(compression)
    data = np.asarray(data, dtype=np.float64)
    if data.ndim != 2 and data.ndim != 3:
        raise ut.UserError("An image must have two or three dimensions rather than {}".format(data.ndim))

    with _newBundle(bundlePath, "FITS file (input)") as bundle:
        _addDataset(bundle, "data", data, np.float64, compression)

# -----------------------------------------------------------------

## This function writes a FITS file bundle for output, the equivalent of a FITS file written by SKIRT: a
# two-dimensional frame or a three-dimensional data cube, with the grid points of its axes and the units
# of the data and of the axes.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em data: a two-dimensional array with shape (ny, nx) or a three-dimensional array with shape (nz, ny, nx).
#    This is the index order of the array obtained when reading a FITS file with the \c astropy.io.fits package.
#  - \em dataUnit: the unit string for the data values.
#  - \em x, \em y, \em z: one-dimensional arrays with the grid points of each axis, or None for an axis that has no
#    grid points. For a two-dimensional frame, \em z must be omitted or have a single grid point.
#  - \em xUnit, \em yUnit, \em zUnit: the unit strings for the grid points of each axis.
#  - \em quantities: an optional dictionary that provides the name of the physical quantity represented by the data
#    and by the axes, as used by SKIRT (for example "wavelength" or "surfacebrightness"), for one or more of the
#    keys "data", "x", "y" and "z".
#  - \em description: an optional free-form string that describes the contents of the bundle.
#  - \em observer: an optional dictionary describing the observer of a frame or data cube recorded by a distant
#    instrument, with one or more of the keys "inclination", "azimuth", "roll" (in degrees), "redshift",
#    "luminosity_distance", "angular_diameter_distance" (numbers), and "distance_unit" (a unit string).
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds the datasets \c data, \c x, \c y and \c z (the latter three only if grid points are provided),
# each with the attributes \c unit and \c quantity (if provided). The description and the observer information
# become attributes of the bundle. A data cube with a single frame is stored as a two-dimensional frame, with
# a single grid point for the z-axis.
#
# SKIRT can also read the bundle as input, in the same way as a FITS file bundle for input: it then uses just the
# \c data dataset.
def writeFitsOutputBundle(bundlePath, data, dataUnit, x=None, xUnit=None, y=None, yUnit=None, z=None, zUnit=None, *,
                          quantities=None, description=None, observer=None, compression=None):
    compression = _verifiedCompression(compression)

    # verify the data and the axes
    data = np.asarray(data, dtype=np.float64)
    if data.ndim != 2 and data.ndim != 3:
        raise ut.UserError("A frame or data cube must have two or three dimensions rather than {}".format(data.ndim))
    if data.ndim == 3 and data.shape[0] == 1 and (z is None or len(z) == 1):
        data = data[0]
    axes = {}
    for name, grid, unit, size in (("x", x, xUnit, data.shape[-1]), ("y", y, yUnit, data.shape[-2]),
                                   ("z", z, zUnit, data.shape[0] if data.ndim == 3 else 1)):
        if grid is not None:
            grid = np.asarray(grid, dtype=np.float64)
            if grid.ndim != 1 or len(grid) != size:
                raise ut.UserError("The {}-axis must have {} grid point(s) to match the data".format(name, size))
            axes[name] = (grid, unit)

    # verify the optional dictionaries
    quantities = dict(quantities) if quantities else {}
    if not set(quantities) <= set(("data", "x", "y", "z")):
        raise ut.UserError("Unsupported key in quantities for frame or data cube: {}".format(", ".join(quantities)))
    observer = dict(observer) if observer else {}
    if not set(observer) <= set(_observerNumbers + _observerStrings):
        raise ut.UserError("Unsupported key in observer information: {}".format(", ".join(observer)))
    for key in _observerNumbers:
        if key in observer: observer[key] = float(observer[key])
    for key in _observerStrings:
        if key in observer: observer[key] = str(observer[key])

    with _newBundle(bundlePath, "FITS file (output)") as bundle:
        _setAttributes(bundle, description=description, **observer)
        for name, (grid, unit) in axes.items():
            _addDataset(bundle, name, grid, np.float64, None, quantity=quantities.get(name), unit=unit)
        _addDataset(bundle, "data", data, np.float64, compression, quantity=quantities.get("data"), unit=dataUnit)

# -----------------------------------------------------------------

## This function writes a spatial grid plot file bundle, the equivalent of a text file with the data to plot a
# spatial grid, as written by the SpatialGridPlotProbe in SKIRT. The plot consists of polylines, in two dimensions
# (for a cut through the grid along a coordinate plane) or in three dimensions.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em points: a two-dimensional array with shape (numPoints, 2) or (numPoints, 3), listing the coordinates of
#    the consecutive points of all polylines.
#  - \em moveto: a one-dimensional array of Booleans with an item for each point: True for the first point of a
#    polyline, and False for a point that is connected to the previous point. The item for the first point must
#    be True.
#  - \em unit: the unit string for the coordinates, or omitted or None if the unit is not known.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds the datasets \c points and \c moveto; the unit becomes an attribute of the bundle.
def writeGridPlotBundle(bundlePath, points, moveto, *, unit=None, compression=None):
    compression = _verifiedCompression(compression)
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ut.UserError("The points of a spatial grid plot must be an array with shape (N, 2) or (N, 3)")
    moveto = np.asarray(moveto)
    if moveto.size == 0: moveto = moveto.astype(bool)    # an empty sequence does not have a Boolean type
    if moveto.dtype != bool or moveto.shape != (len(points),):
        raise ut.UserError("A spatial grid plot must have a Boolean 'moveto' flag for each point")
    if len(moveto) > 0 and not moveto[0]:
        raise ut.UserError("The first point of a spatial grid plot must start a new polyline")

    with _newBundle(bundlePath, "spatial grid plot file") as bundle:
        _setAttributes(bundle, unit=unit)
        _addDataset(bundle, "points", points, np.float64, compression)
        _addDataset(bundle, "moveto", moveto, bool, compression)

# -----------------------------------------------------------------

## This function writes an unstructured text file bundle, the equivalent of a text file that is stored as a whole,
# without interpreting its contents.
#
# The function expects the following arguments:
#  - \em bundlePath: the bundle path, as described in the header of this module.
#  - \em text: a string with the complete contents of the text file.
#  - \em textType: the type of text; one of "free form" (the default; for example a convergence report),
#    "XML" (for example a SKIRT parameter file), or "log" (a SKIRT log file).
#
# The bundle holds a single scalar dataset named \c text, a string of variable length in UTF-8 encoding. The type
# of text becomes the \c type attribute of the bundle.
def writeTextBundle(bundlePath, text, textType="free form"):
    if textType not in _textTypes:
        raise ut.UserError("The type of an unstructured text must be one of: {}".format(", ".join(_textTypes)))

    with _newBundle(bundlePath, "unstructured text file") as bundle:
        _setAttributes(bundle, type=textType)
        bundle.create_dataset("text", data=str(text), dtype=h5py.string_dtype(encoding="utf-8"))

# -----------------------------------------------------------------

## This function returns a tuple with the names of the bundle types, as listed in the header of this module.
def bundleTypes():
    return _bundleTypes

## This function returns the type of the bundle with the specified bundle path, as one of the type names listed in
# the header of this module, or None if the bundle is not recognized as one of the bundles defined by SKIRT. A
# bundle does not record its type, so the function derives the type from the datasets and attributes in the bundle,
# by applying the following rules in this order:
#
#  - \c amr: the bundle holds the datasets \c is_leaf and \c topology.
#  - \c storedtable: at least one of the datasets has an \c axis attribute.
#  - \c columns: at least one of the datasets has a \c column attribute.
#  - \c text: the bundle holds a dataset \c text that contains a string.
#  - \c gridplot: the bundle holds the datasets \c points and \c moveto.
#  - \c fits: the bundle holds a dataset \c data.
#
# The function raises a user error if the bundle path does not refer to a bundle.
def inferBundleType(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        datasets = { name: bundle[name] for name in bundle if bundle.get(name, getclass=True) is h5py.Dataset }

        if "is_leaf" in datasets and "topology" in datasets: return "amr"
        if any([ "axis" in dataset.attrs for dataset in datasets.values() ]): return "storedtable"
        if any([ "column" in dataset.attrs for dataset in datasets.values() ]): return "columns"
        if "text" in datasets and _holdsString(datasets["text"]): return "text"
        if "points" in datasets and "moveto" in datasets: return "gridplot"
        if "data" in datasets: return "fits"
        return None

## This function returns a list with the bundle paths of the bundles specified by its argument, in the order of
# their names. The argument has the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt> of a bundle path, described in
# the header of this module, and is interpreted in one of the following ways:
#
#  - If it refers to an HDF5 file, or to a suite in an HDF5 file, the function returns all bundles in the root group
#    of that file, or in that suite. This is the form of the value of the \c -i and \c -o command line options of
#    SKIRT. The bundles in nested suites are not included.
#  - If it refers to a bundle, the function returns just that bundle.
#  - Otherwise, the last segment of the path within the file can include wildcards, as supported by the standard
#    \c fnmatch module. The function then returns the bundles with a name that matches this pattern, in the root
#    group of the file or in the suite specified by the preceding segments. The list is empty if there are no
#    such bundles.
#
# A group in an HDF5 file is considered to be a suite if it holds at least one other group, and to be a bundle
# otherwise; the root group of the file is always a suite. The function raises a user error if the HDF5 file does
# not exist, or if the argument does not refer to a suite, to a bundle, or to a pattern in an existing suite.
def findBundles(location):
    filePath, groupPath = splitBundlePath(location)
    filePath = ut.absPath(filePath)

    with _openFile(filePath) as file:
        # a suite or a bundle
        target = _groupAt(file, groupPath)
        if target is not None:
            if groupPath and not _isSuite(target): return [ joinBundlePath(filePath, groupPath) ]
            suitePath = groupPath
            names = _bundleNames(target)

        # a pattern in a suite
        else:
            suitePath, _, pattern = groupPath.rpartition("/")
            suite = _groupAt(file, suitePath) if _wildcardSyntax.search(pattern) else None
            if suite is None or (suitePath and not _isSuite(suite)):
                raise ut.UserError("There is no bundle or suite '{}' in HDF5 file {}".format(groupPath, filePath))
            names = [ name for name in _bundleNames(suite) if fnmatch.fnmatchcase(name, pattern) ]

    return [ joinBundlePath(filePath, suitePath, name) for name in names ]

# -----------------------------------------------------------------

## This function reads a text column file bundle, the equivalent of a SKIRT column text file, and returns a
# dictionary with the following keys, which are the names of the arguments of the writeColumnBundle() function:
#  - \em columns: a list with a one-dimensional array of 64-bit floating point numbers for each column, in the
#    order of the columns; all arrays have the same length, the number of rows.
#  - \em names: a list of strings with the description of each column. This is the name of the corresponding
#    dataset, after undoing the adjustments described for the datasetNames() function where possible: a vertical
#    bar is replaced by a forward slash, and the sequence number added to distinguish columns with the same
#    description is removed.
#  - \em units: a list with the SKIRT unit string for each column, or None for a column without a unit, in which
#    case SKIRT uses its default unit for the column.
#  - \em description: the free-form string that describes the contents of the bundle, or None.
#
# The columns are the datasets in the bundle with a \c column attribute, which provides the zero-based index of the
# column. These indices must be consecutive and start at zero, and the datasets must be one-dimensional and have
# the same length. Other datasets are ignored, as SKIRT does. The function raises a user error if the bundle path
# does not refer to a bundle, or if the bundle does not conform to these rules.
def readColumnBundle(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        columns, names, units = _readColumns(bundle)
        return dict(columns=columns, names=names, units=units, description=_stringAttribute(bundle, "description"))

# -----------------------------------------------------------------

## This function reads an AMR text file bundle, the equivalent of a SKIRT column text file that describes an
# adaptive mesh, and returns a dictionary with the following keys, which are the names of the arguments of
# the writeAdaptiveMeshBundle() function:
#  - \em isLeaf: a one-dimensional array of Booleans with an item for each node of the mesh, in the order described
#    for the writeAdaptiveMeshBundle() function: True for a leaf node (a cell) and False for a nonleaf node.
#  - \em topology: a two-dimensional array of integers with shape (numNonLeafNodes, 3), listing the number of
#    subdivisions in the x, y, and z direction for each nonleaf node, in the order in which these nodes occur.
#  - \em columns, \em names, \em units: the properties of the leaf nodes, as described for the readColumnBundle()
#    function, with a row for each leaf node, in the order in which these nodes occur.
#  - \em description: the free-form string that describes the contents of the bundle, or None.
#
# The function raises a user error if the bundle path does not refer to a bundle, if the datasets do not have the
# expected form, or if the number of rows in the \c topology dataset or in the columns does not match the number of
# nonleaf or leaf nodes. It does not verify that the nodes form a complete tree.
def readAdaptiveMeshBundle(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        # read the structure of the tree
        isLeaf = _booleanData(_dataset(bundle, "is_leaf"))
        if isLeaf.ndim != 1:
            raise ut.UserError("Dataset 'is_leaf' in bundle {} is not one-dimensional".format(_label(bundle)))
        topology = _integerData(_dataset(bundle, "topology"))
        if topology.size == 0: topology = topology.reshape(0, 3)
        if topology.ndim != 2 or topology.shape[1] != 3:
            raise ut.UserError("Dataset 'topology' in bundle {} does not have shape (N, 3)".format(_label(bundle)))
        numLeaves = int(np.count_nonzero(isLeaf))
        if len(topology) != len(isLeaf) - numLeaves:
            raise ut.UserError("Dataset 'topology' in bundle {} has {} rows for {} nonleaf nodes" \
                               .format(_label(bundle), len(topology), len(isLeaf) - numLeaves))

        # read the properties of the leaf nodes
        columns, names, units = _readColumns(bundle, reserved=("is_leaf", "topology"))
        if len(columns[0]) != numLeaves:
            raise ut.UserError("The columns in bundle {} have {} rows for {} leaf nodes" \
                               .format(_label(bundle), len(columns[0]), numLeaves))
        return dict(isLeaf=isLeaf, topology=topology, columns=columns, names=names, units=units,
                    description=_stringAttribute(bundle, "description"))

# -----------------------------------------------------------------

## This function reads a stored table bundle, the equivalent of a file in SKIRT stored table format, and returns
# a dictionary with the following keys, which are the names of the arguments of the writeStoredTableBundle()
# function and of the pts.storedtable.io.writeStoredTable() function:
#  - \em axisNames, \em axisUnits, \em axisScales: lists of strings with the name, the unit, and the scale ('lin'
#    or 'log') of each axis, in the order of the axes.
#  - \em axisGrids: a list with a one-dimensional array of grid points for each axis.
#  - \em quantityNames, \em quantityUnits, \em quantityScales: lists of strings with the name, the unit, and the
#    scale ('lin' or 'log') of each quantity.
#  - \em quantityValues: a list with an array of values for each quantity; each array has indices in the same
#    order and with the same range as the axes.
#
# The axes are the datasets with an \c axis attribute, which provides the zero-based index of the axis; all other
# datasets in the bundle are quantities. A bundle does not record the order of its quantities, so they are listed
# in the alphabetical order of the names of their datasets; SKIRT locates a quantity by its name. The names are
# the names of the datasets, after undoing the adjustments described for the datasetNames() function where possible.
#
# The function raises a user error if the bundle path does not refer to a bundle, if the indices of the axes are
# not consecutive integers starting at zero, if there is no quantity, if a dataset lacks its \c unit or \c log
# attribute, or if the shape of a quantity does not match the axes.
def readStoredTableBundle(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        # distinguish the axes from the quantities
        axes = {}
        quantities = []
        for name in sorted(bundle):
            if bundle.get(name, getclass=True) is not h5py.Dataset: continue
            index = _integerAttribute(bundle[name], "axis")
            if index is None: quantities.append(bundle[name])
            elif index in axes:
                raise ut.UserError("Two datasets in bundle {} have axis index {}".format(_label(bundle), index))
            else: axes[index] = bundle[name]
        if len(axes) == 0 or sorted(axes) != list(range(len(axes))):
            raise ut.UserError("The axis indices of the datasets in bundle {} are not consecutive integers " \
                               "starting at zero".format(_label(bundle)))
        if len(quantities) == 0:
            raise ut.UserError("Bundle {} does not hold a dataset for a quantity".format(_label(bundle)))
        axes = [ axes[index] for index in range(len(axes)) ]

        # read the grids and the values
        axisGrids = [ _floatData(dataset) for dataset in axes ]
        for dataset, grid in zip(axes, axisGrids):
            if grid.ndim != 1: raise ut.UserError("The {} is not one-dimensional".format(_where(dataset)))
        shape = tuple([ len(grid) for grid in axisGrids ])
        quantityValues = [ _floatData(dataset) for dataset in quantities ]
        for dataset, values in zip(quantities, quantityValues):
            if values.shape != shape:
                raise ut.UserError("The {} has shape {}, which does not match the shape {} of the axes" \
                                   .format(_where(dataset), values.shape, shape))

        # read the metadata
        names = _descriptions([ dataset.name.split("/")[-1] for dataset in axes + quantities ])
        units = [ _requiredAttribute(dataset, "unit", _stringAttribute) for dataset in axes + quantities ]
        scales = [ "log" if _requiredAttribute(dataset, "log", _booleanAttribute) else "lin" \
                   for dataset in axes + quantities ]
        numAxes = len(axes)
        return dict(axisNames=names[:numAxes], axisUnits=units[:numAxes], axisScales=scales[:numAxes],
                    axisGrids=axisGrids, quantityNames=names[numAxes:], quantityUnits=units[numAxes:],
                    quantityScales=scales[numAxes:], quantityValues=quantityValues)

# -----------------------------------------------------------------

## This function reads a FITS file bundle, for input or for output, the equivalent of a FITS file with an image, a
# frame, or a data cube. It returns a dictionary with the following keys, which are the names of the arguments of
# the writeFitsOutputBundle() function:
#  - \em data: a two-dimensional array with shape (ny, nx) or a three-dimensional array with shape (nz, ny, nx).
#    This is the index order of the array obtained when reading a FITS file with the \c astropy.io.fits package.
#  - \em dataUnit: the unit string for the data values, or None.
#  - \em x, \em y, \em z: one-dimensional arrays with the grid points of each axis, or None for an axis without
#    grid points. For a two-dimensional array, \em z is None or has a single grid point; the latter indicates
#    a data cube with a single frame.
#  - \em xUnit, \em yUnit, \em zUnit: the unit strings for the grid points of each axis, or None.
#  - \em quantities: a dictionary that provides the name of the physical quantity represented by the data and by
#    the axes, as far as recorded in the bundle, for the keys "data", "x", "y" and "z".
#  - \em description: the free-form string that describes the contents of the bundle, or None.
#  - \em observer: a dictionary describing the observer, with the keys listed for the writeFitsOutputBundle()
#    function, as far as recorded in the bundle.
#
# A FITS file bundle for input holds just the \c data dataset. For such a bundle, all items other than \em data
# are None or an empty dictionary.
#
# If the optional \em readData argument is False, the function does not read the values of the \c data dataset,
# and the \em data item is None. This allows obtaining the axes and attributes of a large data cube without reading
# the data cube itself.
#
# The function raises a user error if the bundle path does not refer to a bundle, if the bundle has no dataset
# \c data with two or three dimensions, or if the number of grid points of an axis does not match the data.
def readFitsBundle(bundlePath, *, readData=True):
    with _existingBundle(bundlePath) as bundle:
        # read the data
        dataset = _dataset(bundle, "data")
        if dataset.ndim != 2 and dataset.ndim != 3:
            raise ut.UserError("The {} has {} dimension(s) rather than two or three" \
                               .format(_where(dataset), dataset.ndim))
        data = _floatData(dataset) if readData else None
        result = dict(data=data, dataUnit=_stringAttribute(dataset, "unit"))
        quantities = dict(data=_stringAttribute(dataset, "quantity"))

        # read the axes
        shape = dataset.shape
        sizes = dict(x=shape[-1], y=shape[-2], z=shape[0] if len(shape) == 3 else 1)
        for name, size in sizes.items():
            grid = unit = None
            if bundle.get(name, getclass=True) is h5py.Dataset:
                grid = _floatData(bundle[name])
                if grid.ndim != 1 or len(grid) != size:
                    raise ut.UserError("The {} does not have {} grid point(s) to match the data" \
                                       .format(_where(bundle[name]), size))
                unit = _stringAttribute(bundle[name], "unit")
                quantities[name] = _stringAttribute(bundle[name], "quantity")
            result[name] = grid
            result[name + "Unit"] = unit

        # read the attributes of the bundle
        observer = { key: _numberAttribute(bundle, key) for key in _observerNumbers }
        observer.update({ key: _stringAttribute(bundle, key) for key in _observerStrings })
        result["quantities"] = { key: value for key, value in quantities.items() if value is not None }
        result["description"] = _stringAttribute(bundle, "description")
        result["observer"] = { key: value for key, value in observer.items() if value is not None }
        return result

# -----------------------------------------------------------------

## This function reads a spatial grid plot file bundle, the equivalent of a text file with the data to plot a
# spatial grid, and returns a dictionary with the following keys, which are the names of the arguments of the
# writeGridPlotBundle() function:
#  - \em points: a two-dimensional array with shape (numPoints, 2) or (numPoints, 3), listing the coordinates of
#    the consecutive points of all polylines.
#  - \em moveto: a one-dimensional array of Booleans with an item for each point: True for the first point of a
#    polyline, and False for a point that is connected to the previous point.
#  - \em unit: the unit string for the coordinates, or None if the bundle does not record the unit.
#
# The function raises a user error if the bundle path does not refer to a bundle, or if the bundle does not hold
# these two datasets in the expected form.
def readGridPlotBundle(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        points = _floatData(_dataset(bundle, "points"))
        if points.ndim != 2 or points.shape[1] not in (2, 3):
            raise ut.UserError("Dataset 'points' in bundle {} does not have shape (N, 2) or (N, 3)" \
                               .format(_label(bundle)))
        moveto = _booleanData(_dataset(bundle, "moveto"))
        if moveto.shape != (len(points),):
            raise ut.UserError("Dataset 'moveto' in bundle {} does not have an item for each of the {} points" \
                               .format(_label(bundle), len(points)))
        return dict(points=points, moveto=moveto, unit=_stringAttribute(bundle, "unit"))

# -----------------------------------------------------------------

## This function reads an unstructured text file bundle, the equivalent of a text file that is stored as a whole,
# and returns a dictionary with the following keys, which are the names of the arguments of the writeTextBundle()
# function:
#  - \em text: a string with the complete contents of the text file.
#  - \em textType: the type of text recorded in the bundle, normally one of "free form", "XML", or "log"; or
#    "free form" if the bundle does not record the type.
#
# The function raises a user error if the bundle path does not refer to a bundle, or if the bundle does not hold a
# dataset \c text with a single string.
def readTextBundle(bundlePath):
    with _existingBundle(bundlePath) as bundle:
        dataset = _dataset(bundle, "text")
        if dataset.shape != () or not _holdsString(dataset):
            raise ut.UserError("The {} does not hold a single string".format(_where(dataset)))
        text = dataset[()]
        if isinstance(text, bytes): text = text.decode("utf-8", errors="replace")
        textType = _stringAttribute(bundle, "type")
        return dict(text=str(text), textType=textType if textType is not None else "free form")

# -----------------------------------------------------------------

## This context manager creates a new, empty bundle with the specified bundle path, replacing an existing bundle
# with the same path, and adds the attributes that describe the origin of the bundle. It yields the HDF5 group that
# represents the bundle, so that the caller can add attributes and datasets. If an exception is raised before the
# bundle is complete, the bundle is removed again. The HDF5 file is closed when the context is left. The \em kind
# argument describes the type of bundle in the message that is logged when the bundle has been written.
@contextlib.contextmanager
def _newBundle(bundlePath, kind):
    filePath, groupPath = splitBundlePath(bundlePath)
    if not groupPath:
        raise ut.UserError("Path does not specify the name of a bundle in the HDF5 file: '{}'".format(bundlePath))
    filePath = ut.absPath(filePath)
    bundlePath = joinBundlePath(filePath, groupPath)

    # open the file, creating it if it does not exist
    try:
        file = h5py.File(filePath, "a")
    except OSError as error:
        raise ut.UserError("Cannot open or create HDF5 file '{}': {}".format(filePath, error)) from None

    with file:
        _removeBundle(file, groupPath, bundlePath)
        bundle = file.create_group(groupPath)
        try:
            created = datetime.datetime.now().isoformat(timespec="milliseconds")
            _setAttributes(bundle, producer=_producer, format=_format, created=created)
            yield bundle
        except BaseException:
            del file[groupPath]
            raise

    logging.info("Created {} bundle: {}".format(kind, bundlePath))

## This helper function prepares the specified open HDF5 file for creating a bundle with the specified path within
# the file: it removes an existing bundle with that path, and it raises a user error if the bundle cannot be created
# without damaging other contents of the file. The last argument is the bundle path used in error messages.
def _removeBundle(file, groupPath, bundlePath):
    # each of the groups leading up to the bundle must be a suite; as a sanity check that remains fast for a suite
    # holding many bundles, we verify just the first member of the group, which must not be a dataset
    segments = groupPath.split("/")
    for depth in range(1, len(segments)):
        suitePath = "/".join(segments[:depth])
        if suitePath in file:
            suite = file[suitePath]
            first = next(iter(suite), None) if isinstance(suite, h5py.Group) else None
            if not isinstance(suite, h5py.Group) or (first is not None and suite.get(first, getclass=True) \
                                                                           is h5py.Dataset):
                raise ut.UserError("Cannot create bundle {} because '{}' is not a suite" \
                                   .format(bundlePath, suitePath))

    # an existing object with the same path must be a bundle, i.e. a group that does not hold other groups
    if groupPath in file:
        existing = file[groupPath]
        if not isinstance(existing, h5py.Group):
            raise ut.UserError("Cannot create bundle {} because there is a dataset with that path" \
                               .format(bundlePath))
        if any([ existing.get(name, getclass=True) is h5py.Group for name in existing ]):
            raise ut.UserError("Cannot create bundle {} because there is a suite with that path" \
                               .format(bundlePath))
        del file[groupPath]

## This helper function adds a dataset with the specified name, data, and data type to the specified bundle,
# compressed at the specified level (None for no compression), and sets the attributes of the dataset as described
# for the _setAttributes() function.
def _addDataset(bundle, name, data, dtype, compression, **attributes):
    data = np.ascontiguousarray(data, dtype=dtype)
    options = {}
    if compression is not None and data.size > 0:
        options = dict(compression="gzip", compression_opts=compression, shuffle=True)
    dataset = bundle.create_dataset(name, data=data, dtype=dtype, **options)
    _setAttributes(dataset, **attributes)

## This helper function sets the attributes of the specified bundle or dataset to the values of the corresponding
# keyword arguments. An attribute with value None is not set. The attributes are stored with the data types
# expected by SKIRT: Booleans in the form used by \c h5py, integers as 32-bit integers, floating point numbers
# as 64-bit numbers, and other values as strings of variable length in UTF-8 encoding.
def _setAttributes(target, **attributes):
    for name, value in attributes.items():
        if value is None: continue
        if isinstance(value, (bool, np.bool_)): value = np.bool_(value)
        elif isinstance(value, (int, np.integer)): value = np.int32(value)
        elif isinstance(value, (float, np.floating)): value = np.float64(value)
        else: value = str(value)
        target.attrs[name] = value

## This helper function removes leading, trailing, and repeated forward slashes from a path within an HDF5 file.
def _cleanGroupPath(groupPath):
    return "/".join([ segment for segment in groupPath.split("/") if segment ])

## This helper function verifies the specified compression level and returns it as an integer from 1 to 9, or as
# None if the data must not be compressed.
def _verifiedCompression(compression):
    if compression is None or compression == 0: return None
    if isinstance(compression, (bool, np.bool_)) or not isinstance(compression, (int, np.integer)) \
                    or not 1 <= compression <= 9:
        raise ut.UserError("The compression level must be an integer from 1 to 9, or zero for no compression")
    return int(compression)

## This helper function verifies the columns, column descriptions, and unit strings specified for a table, and
# returns them as a tuple with three lists of the same length: a one-dimensional array of 64-bit floating point
# numbers, a description string, and a unit string or None for each column.
def _verifiedColumns(columns, names, units):
    columns = [ np.asarray(column, dtype=np.float64) for column in columns ]
    if len(columns) < 1:
        raise ut.UserError("A table must have at least one column")
    if any([ column.ndim != 1 or len(column) != len(columns[0]) for column in columns ]):
        raise ut.UserError("The columns of a table must be one-dimensional arrays of the same length")

    names = [names] if isinstance(names, str) else list(names)
    if len(names) != len(columns):
        raise ut.UserError("A table with {} column(s) cannot have {} column description(s)" \
                           .format(len(columns), len(names)))

    if units is None: units = [None] * len(columns)
    units = [units] if isinstance(units, str) else list(units)
    if len(units) != len(columns):
        raise ut.UserError("A table with {} column(s) cannot have {} unit(s)".format(len(columns), len(units)))
    return columns, names, units

## This helper function returns True if the specified scale string of a stored table is 'log', False if it
# is 'lin', and raises a user error otherwise.
def _isLogScale(scale):
    if scale == "log": return True
    if scale == "lin": return False
    raise ut.UserError("The scale of an axis or quantity in a stored table must be 'lin' or 'log' rather than '{}'" \
                       .format(scale))

# -----------------------------------------------------------------

## This context manager opens the HDF5 file with the specified absolute path for reading, and yields the open file.
# The file is closed when the context is left. A user error is raised if the file does not exist or cannot be opened.
@contextlib.contextmanager
def _openFile(filePath):
    if not filePath.is_file():
        raise ut.UserError("HDF5 file does not exist: '{}'".format(filePath))
    try:
        file = h5py.File(filePath, "r")
    except OSError as error:
        raise ut.UserError("Cannot open HDF5 file '{}': {}".format(filePath, error)) from None
    with file:
        yield file

## This context manager opens the existing bundle with the specified bundle path for reading, and yields the HDF5
# group that represents the bundle. The HDF5 file is closed when the context is left. A user error is raised if
# the HDF5 file cannot be opened, or if the path does not refer to a bundle in that file.
@contextlib.contextmanager
def _existingBundle(bundlePath):
    filePath, groupPath = splitBundlePath(bundlePath)
    if not groupPath:
        raise ut.UserError("Path does not specify the name of a bundle in the HDF5 file: '{}'".format(bundlePath))
    filePath = ut.absPath(filePath)

    with _openFile(filePath) as file:
        bundle = _groupAt(file, groupPath)
        if bundle is None:
            raise ut.UserError("There is no bundle '{}' in HDF5 file {}".format(groupPath, filePath))
        if _isSuite(bundle):
            raise ut.UserError("'{}' in HDF5 file {} is a suite rather than a bundle".format(groupPath, filePath))
        yield bundle

## This helper function returns the group with the specified path within the specified open HDF5 file, or the root
# group of the file for an empty path. It returns None if there is no object with that path, and raises a user
# error if the object is not a group.
def _groupAt(file, groupPath):
    if not groupPath: return file
    try:
        target = file.get(groupPath)
    except (KeyError, TypeError, ValueError):
        target = None
    if target is None: return None
    if not isinstance(target, h5py.Group):
        raise ut.UserError("'{}' in HDF5 file {} is a dataset rather than a bundle or a suite" \
                           .format(groupPath, file.filename))
    return target

## This helper function returns True if the specified group is a suite, i.e. if it holds at least one other group.
def _isSuite(group):
    return any([ group.get(name, getclass=True) is h5py.Group for name in group ])

## This helper function returns a sorted list with the names of the bundles in the specified suite, i.e. of the
# groups in the suite that are not suites themselves.
def _bundleNames(suite):
    return sorted([ name for name in suite if suite.get(name, getclass=True) is h5py.Group \
                                              and not _isSuite(suite[name]) ])

## This helper function returns the bundle path of the specified group in an open HDF5 file, for use in messages.
def _label(bundle):
    return joinBundlePath(bundle.file.filename, bundle.name)

## This helper function returns a phrase that identifies the specified bundle or dataset, for use in messages.
def _where(target):
    if isinstance(target, h5py.Dataset):
        return "dataset '{}' in bundle {}".format(target.name.split("/")[-1], _label(target.parent))
    return "bundle {}".format(_label(target))

## This helper function returns the dataset with the specified name in the specified bundle, or raises a user error
# if the bundle does not hold such a dataset.
def _dataset(bundle, name):
    if bundle.get(name, getclass=True) is not h5py.Dataset:
        raise ut.UserError("Bundle {} does not hold a dataset '{}'".format(_label(bundle), name))
    return bundle[name]

## This helper function returns True if the specified dataset holds strings.
def _holdsString(dataset):
    return h5py.check_string_dtype(dataset.dtype) is not None or dataset.dtype.kind in "SU"

## This helper function returns the contents of the specified dataset as an array of 64-bit floating point numbers,
# or raises a user error if the dataset does not hold numbers.
def _floatData(dataset):
    if dataset.dtype.kind not in "fiu":
        raise ut.UserError("The {} does not hold numbers".format(_where(dataset)))
    return np.asarray(dataset[()], dtype=np.float64)

## This helper function returns the contents of the specified dataset as an array of 64-bit integers, or raises a
# user error if the dataset does not hold integers.
def _integerData(dataset):
    if dataset.dtype.kind not in "iu":
        raise ut.UserError("The {} does not hold integers".format(_where(dataset)))
    return np.asarray(dataset[()], dtype=np.int64)

## This helper function returns the contents of the specified dataset as an array of Booleans, or raises a user
# error if the dataset does not hold Booleans. A dataset that holds integers is accepted as well, with zero
# representing False and any other value True.
def _booleanData(dataset):
    if dataset.dtype.kind == "b": return np.asarray(dataset[()], dtype=bool)
    if dataset.dtype.kind in "iu": return np.asarray(dataset[()]) != 0
    raise ut.UserError("The {} does not hold Booleans".format(_where(dataset)))

## This helper function reads the columns of a text column file bundle or of an AMR text file bundle, i.e. the
# datasets with a "column" attribute, ignoring the datasets with one of the specified reserved names. It returns a
# tuple with three lists of the same length: an array of 64-bit floating point numbers, a description string, and a
# unit string or None for each column, in the order of the columns.
def _readColumns(bundle, reserved=()):
    # locate the dataset for each column
    datasets = {}
    for name in bundle:
        if name in reserved or bundle.get(name, getclass=True) is not h5py.Dataset: continue
        index = _integerAttribute(bundle[name], "column")
        if index is None: continue
        if index in datasets:
            raise ut.UserError("Two datasets in bundle {} have column index {}".format(_label(bundle), index))
        datasets[index] = bundle[name]
    if len(datasets) == 0:
        raise ut.UserError("Bundle {} does not hold any datasets with a column index".format(_label(bundle)))
    if sorted(datasets) != list(range(len(datasets))):
        raise ut.UserError("The column indices of the datasets in bundle {} are not consecutive integers " \
                           "starting at zero".format(_label(bundle)))
    datasets = [ datasets[index] for index in range(len(datasets)) ]

    # read the columns
    columns = [ _floatData(dataset) for dataset in datasets ]
    for dataset, column in zip(datasets, columns):
        if column.ndim != 1:
            raise ut.UserError("The {} is not one-dimensional".format(_where(dataset)))
        if len(column) != len(columns[0]):
            raise ut.UserError("The {} has {} rows, while the first column has {} rows" \
                               .format(_where(dataset), len(column), len(columns[0])))
    names = _descriptions([ dataset.name.split("/")[-1] for dataset in datasets ], reserved)
    units = [ _stringAttribute(dataset, "unit") for dataset in datasets ]
    return columns, names, units

## This helper function converts a list of dataset names, in the order of the columns, axes, or quantities they
# represent, to a list of the corresponding descriptions, by undoing the adjustments made by the datasetNames()
# function: a vertical bar is replaced by a forward slash, and a sequence number between parentheses is removed
# if the name without that number occurs earlier in the list or is one of the specified reserved names. A generic
# name that replaces an empty description cannot be recognized as such, and is left as it is.
def _descriptions(names, reserved=()):
    descriptions = []
    used = set(reserved)
    for name in names:
        match = _sequenceNumberSyntax.fullmatch(name)
        description = match.group(1) if match and match.group(1) in used else name
        used.add(name)
        descriptions.append(description.replace("|", "/"))
    return descriptions

## This helper function returns the value of the attribute with the specified name of the specified bundle or
# dataset, or None if there is no such attribute. A string is returned as a Python string, and an array with a
# single item as that item.
def _attribute(target, name):
    if name not in target.attrs: return None
    value = target.attrs[name]
    if isinstance(value, np.ndarray) and value.size == 1: value = value.reshape(-1)[0]
    if isinstance(value, bytes): value = value.decode("utf-8", errors="replace")
    return value

## This helper function returns the value of the specified attribute as a string, or None if there is no such
# attribute. It raises a user error if the value is not a string.
def _stringAttribute(target, name):
    value = _attribute(target, name)
    if value is None: return None
    if isinstance(value, str): return str(value)
    raise ut.UserError("The attribute '{}' of the {} is not a string".format(name, _where(target)))

## This helper function returns the value of the specified attribute as an integer, or None if there is no such
# attribute. It raises a user error if the value is not an integer.
def _integerAttribute(target, name):
    value = _attribute(target, name)
    if value is None: return None
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)): return int(value)
    if isinstance(value, (float, np.floating)) and float(value).is_integer(): return int(value)
    raise ut.UserError("The attribute '{}' of the {} is not an integer".format(name, _where(target)))

## This helper function returns the value of the specified attribute as a floating point number, or None if there
# is no such attribute. It raises a user error if the value is not a number.
def _numberAttribute(target, name):
    value = _attribute(target, name)
    if value is None: return None
    if isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, (bool, np.bool_)):
        return float(value)
    raise ut.UserError("The attribute '{}' of the {} is not a number".format(name, _where(target)))

## This helper function returns the value of the specified attribute as a Boolean, or None if there is no such
# attribute. It raises a user error if the value is not a Boolean; the integers zero and one are accepted as well.
def _booleanAttribute(target, name):
    value = _attribute(target, name)
    if value is None: return None
    if isinstance(value, (bool, np.bool_)): return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1): return bool(value)
    raise ut.UserError("The attribute '{}' of the {} is not a Boolean".format(name, _where(target)))

## This helper function returns the value of the specified attribute, obtained through the specified one of the
# preceding functions, or raises a user error if there is no such attribute.
def _requiredAttribute(target, name, getter):
    value = getter(target, name)
    if value is None:
        raise ut.UserError("The {} has no '{}' attribute".format(_where(target), name))
    return value

# -----------------------------------------------------------------
