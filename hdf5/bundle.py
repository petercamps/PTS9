#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.hdf5.bundle Writing SKIRT bundles to an HDF5 file
#
# The functions in this module write data held in memory to an HDF5 file, in the form expected by SKIRT when it
# reads its input from an HDF5 file, and produced by SKIRT when it writes its output to an HDF5 file. They can be
# used, for example, by a script that prepares the input for a large number of SKIRT simulations.
#
# Bundles and suites
# ------------------
#
# The contents of a SKIRT input or output file are represented in an HDF5 file by a \em bundle: an HDF5 group
# holding one or more datasets, with dataset names and attributes defined by SKIRT for each type of file. There is a
# function in this module for each type of bundle:
#
# | Function | Bundle type | Corresponding regular file
# |----------|-------------|---------------------------
# | writeColumnBundle() | Text column file bundle | Column text file (input and output)
# | writeAdaptiveMeshBundle() | AMR text file bundle | Column text file describing an adaptive mesh (input)
# | writeStoredTableBundle() | Stored table bundle | SKIRT stored table file (input)
# | writeFitsInputBundle() | FITS file bundle for input | FITS file holding just an image or a data cube (input)
# | writeFitsOutputBundle() | FITS file bundle for output | FITS file holding a frame or data cube and its axes (output)
# | writeGridPlotBundle() | Spatial grid plot file bundle | Text file with data to plot a spatial grid (output)
# | writeTextBundle() | Unstructured text file bundle | Log file, parameter file, or convergence report (output)
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
# so that they are valid and unique within the bundle. See the datasetNames() function for the rules.
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
# \note These low-level functions do \em not support astropy quantities; they expect plain numpy arrays (or
#       sequences) with values in the units declared by the corresponding unit strings.
#

# -----------------------------------------------------------------

import contextlib
import datetime
import logging
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
#  - \em image: a two-dimensional array with shape (ny, nx) or a three-dimensional array with shape (nz, ny, nx).
#    This is the index order of the array obtained when reading a FITS file with the \c astropy.io.fits package.
#  - \em compression: the optional compression level, as described in the header of this module.
#
# The bundle holds a single dataset named \c image, without attributes.
def writeFitsInputBundle(bundlePath, image, *, compression=None):
    compression = _verifiedCompression(compression)
    image = np.asarray(image, dtype=np.float64)
    if image.ndim != 2 and image.ndim != 3:
        raise ut.UserError("An image must have two or three dimensions rather than {}".format(image.ndim))

    with _newBundle(bundlePath, "FITS file (input)") as bundle:
        _addDataset(bundle, "image", image, np.float64, compression)

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
