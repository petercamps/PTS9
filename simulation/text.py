#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.text Handling SKIRT text input/output files
#
# This module offers functions for reading information from SKIRT text output files (including column text files,
# log file, convergence reports, and spatial grid plot files) and for writing SKIRT text input files (column text
# files).
#
# The functions that read an output file also accept the corresponding bundle in an HDF5 file, specified as a
# pts.utils.BundlePath instance or as a string with a bundle path (see the pts.utils.isBundlePath() function).
# They then return the same information as for the regular file. This requires the \c h5py package.
#

# -----------------------------------------------------------------

import numpy as np
import pts.utils as ut
from .units import unit as smunit

# -----------------------------------------------------------------

## This function extracts a numeric value from a text file such as a SKIRT log file or a file produced by
# SKIRT's ConvergenceInfoProbe class. The value is returned as an astropy quantity with appropriate units,
# assuming that a valid unit string has been found in the text file.
#
# The file path is interpreted as described for the pts.utils.absPath() function. The path can also refer to an
# unstructured text file bundle in an HDF5 file, such as the bundle that holds the convergence report of a simulation
# with HDF5 output.
#
# The text line containing the value is located by the \em trigger (a text string that must occur on a line
# before the one containing the value) and the \em header (a text string that must occur on the line containing
# the field). In fact, the trigger can consist of multiple subtriggers, separated by forward slashes, that
# must be triggered in sequence before the function starts looking for the header.
# If no conforming line is found, the function raises an error.
#
# If the selected line has a section between parentheses at the end, this section (including parentheses) is removed.
# Then the line is split in segments on white space. If the last segment is a valid representation of a floating point
# number, this number is returned as a dimensionless quantity. Otherwise, the last two segments represent the value and
# the unit string, respectively. If there are problems converting these representations, the function raises an error.
#
# For example, the following call will return the total gridded dust mass from a spatial grid convergence file:
#
#     getQuantityFromFile("xxx_convergence.dat", "Dust/Total mass", "Gridded")
#
# and the following call will return the total metallic mass for the first particle medium from a log file:
#
#     getQuantityFromFile("xxx_log.txt", "ParticleMedium", "Total metallic mass")
#
def getQuantityFromFile(path, trigger, header):
    triggers = trigger.split("/")
    triggered = 0
    for line in _textLines(path):
        # select the line
        if triggered < len(triggers) and triggers[triggered] in line: triggered += 1
        elif triggered==len(triggers) and header in line:
            line = line.strip()
            # remove section between parentheses and/or trailing "."
            if line.endswith(')') or line.endswith(').'):
                line = line[:line.rfind('(')]
            if line.endswith('.'):
                line = line[:-1]
            segments = line.split()
            if len(segments) >= 2:
                try:
                    return float(segments[-1]) << smunit("")
                except ValueError: pass
                return float(segments[-2]) << smunit(segments[-1])

    raise ValueError("Quantity '{}' not found in text file".format(header))

## This helper function returns an iterable over the lines of the specified text file, or of the text held by the
# specified unstructured text file bundle.
def _textLines(path):
    if ut.isBundlePath(path):
        import pts.hdf5 as h5
        yield from h5.readTextBundle(str(path))["text"].splitlines()
    else:
        with open(ut.absPath(path)) as infile:
            yield from infile

# -----------------------------------------------------------------

## This function reads the header of a SKIRT column text file, and returns a list of the column descriptions,
# in order of occurrence. If the specified file does not have a column header in SKIRT format, the function
# returns an empty list. The file path is interpreted as described for the pts.utils.absPath() function.
# The path can also refer to a text column file bundle in an HDF5 file.
def getColumnDescriptions(path):
    return [ description for description, unit in _columnHeader(path) ]

## This function loads data from a SKIRT column text file, using the metadata in the file's header to assign
# appropriate units to each column, and optionally, to locate columns by name.
#
# \note The specified file \em must have a column header in SKIRT format, because the information in the header
# is used to determine the appropriate units for each column. If the file does not have a conforming column header,
# use the numpy.loadtxt() function directly.
#
# The file path is interpreted as described for the pts.utils.absPath() function. The path can also refer to a text
# column file bundle in an HDF5 file, which holds the same metadata.
#
# The columns to be loaded can be specified in one of the following ways:
#  - A sequence (list or tuple) of integers specifying zero-based column indices.
#  - A single string including a comma-separated list of one-based column indices and/or description fragments;
#    each fragment must be contained in exactly one column description in the file header.
#  - None or omitted: causes all columns to be loaded.
#
# The function always returns a list of astropy quantities, even if only one column has been loaded. The units for
# each item are derived from the metadata in the file header. If the file contains a single row, each item in the
# returned list is a scalar astropy quantity. Otherwise each item in the list is an astropy quantity array, and
# all arrays have the same length.
#
def loadColumns(path, columns=None):
    # get the column metadata and, for a bundle, the column data
    if ut.isBundlePath(path):
        header, bundleColumns = _readBundleColumns(path)
    else:
        header = _columnHeader(path)

    # construct a list of zero-based indices for columns to be loaded
    if columns is None:
        # None --> all columns
        usecols = list(range(len(header)))
    elif isinstance(columns, str):
        # single string --> comma-separated list of one-based column indices and/or description fragments
        usecols = []
        for col in columns.split(","):
            # string representing an integer -> one-based index
            if _representsInteger(col):
                col = int(col)
                if col<1 or col>len(header):
                    raise ValueError("one-based column index is out of range: {}".format(col))
                usecols.append(col-1)
            # other string --> match column description
            else:
                usecols.append(_indexForDescriptionInHeader(col, header))
    else:
        # sequence --> zero-based column indices
        usecols = []
        for col in columns:
            col = int(col)
            if col<0 or col>=len(header):
                raise ValueError("zero-based column index is out of range: {}".format(col))
            usecols.append(col)

    # load the data and assign units
    if ut.isBundlePath(path):
        data = [ bundleColumns[colindex] for colindex in usecols ]
    else:
        data = np.loadtxt(ut.absPath(path), usecols=usecols, ndmin=2, unpack=True)
    return [coldata << smunit(header[colindex][1]) for colindex, coldata in zip(usecols, data)]

## This helper function returns a list with a tuple (description, unit) for each column of the specified text column
# file or text column file bundle, in order of occurrence. The unit is the empty string for a column without a unit.
def _columnHeader(path):
    # get the column metadata from a bundle
    if ut.isBundlePath(path):
        return _readBundleColumns(path)[0]

    # parse the header of a regular file
    header = []
    with open(ut.absPath(path)) as infile:
        for line in infile:
            # skip empty lines
            line = line.strip()
            if len(line) > 0:
                # handle end-of-header
                if not line.startswith("#"): break;
                # remove hash character and skip non-column header lines
                line = line[1:].strip()
                if line.lower().startswith("column") and ":" in line:
                    # extract the description and the unit
                    colon = line.find(":")
                    if line.endswith(')'):
                        left = line.rfind("(")
                        description = line[colon+1:left].strip()
                        unit = line[left+1:-1].strip()
                    else:
                        description = line[colon+1:].strip()
                        unit = ""
                    # add the column metadata to the internal list
                    header.append( (description, unit) )
    return header

## This helper function reads the specified text column file bundle, and returns a tuple with a list of the column
# metadata, as described for the _columnHeader() function, and a list with an array of values for each column.
# For a bundle of another type, such as an unstructured text file bundle, it returns two empty lists, just as a
# regular text file without column information lines has no columns.
def _readBundleColumns(path):
    import pts.hdf5 as h5
    try:
        table = h5.readColumnBundle(str(path))
    except ut.UserError:
        if h5.inferBundleType(str(path)) not in ("columns", "amr"): return [], []
        raise
    header = [ (name, unit if unit is not None else "") for name, unit in zip(table["names"], table["units"]) ]
    return header, table["columns"]


## This helper function returns true is the specified value can be converted to an integer, and False otherwise.
def _representsInteger(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

## This helper function searches the header metadata for a column with a description containing the specified string.
# If there is exactly one such column, the function returns its zero-based index. Otherwise, it raises an error.
def _indexForDescriptionInHeader(spec, header):
    count = 0
    result = -1
    for index,(description,unit) in enumerate(header):
        if spec.strip().lower() in description.lower():
            count += 1
            result = index

    if count == 0:
        raise ValueError("Column specification does not match header descriptions: {}".format(spec))
    if count > 1:
        raise ValueError("Column specification matches multiple header descriptions: {}".format(spec))
    return result

# -----------------------------------------------------------------

## This function loads the polylines from a spatial grid plot file produced by SKIRT, i.e. a file generated by
# a probe of the spatial grid with the data to plot the grid, and returns a tuple (points, moveto) with the
# following items:
#  - \em points: a two-dimensional numpy array with shape (numPoints, 2) or (numPoints, 3), listing the coordinates
#    of the consecutive points of all polylines, in the output length unit of the simulation. The plot file holds
#    two coordinates for each point if it represents the intersection of the grid with a coordinate plane, and three
#    coordinates if it describes the grid in three dimensions.
#  - \em moveto: a one-dimensional numpy array of Booleans with an item for each point: True for the first point
#    of a polyline, and False for a point that is connected to the preceding point.
#
# The file path is interpreted as described for the pts.utils.absPath() function. The path can also refer to a
# spatial grid plot file bundle in an HDF5 file. A regular plot file has a line with the coordinates of each point,
# separated by white space, and an empty line after each polyline. For an empty plot file, the function returns an
# array of points with shape (0, 0), or for a bundle, the shape of its empty dataset.
#
def loadGridPlot(path):
    # read a bundle
    if ut.isBundlePath(path):
        import pts.hdf5 as h5
        plot = h5.readGridPlotBundle(str(path))
        return plot["points"], plot["moveto"]

    # parse a regular file
    points = []
    moveto = []
    newPolyline = True
    with open(ut.absPath(path)) as infile:
        for line in infile:
            coords = line.split()
            if len(coords) == 0:
                newPolyline = True
            else:
                if len(points) > 0 and len(coords) != len(points[0]):
                    raise ValueError("Spatial grid plot file has points with {} and {} coordinates: {}" \
                                     .format(len(points[0]), len(coords), path))
                points.append([ float(coord) for coord in coords ])
                moveto.append(newPolyline)
                newPolyline = False
    if len(points) == 0: return np.zeros((0, 0)), np.zeros(0, dtype=bool)
    return np.array(points), np.array(moveto)

# -----------------------------------------------------------------

## This function writes data to a SKIRT column text file, generating the metadata in the file's header from the
# information provided to the function, and automatically converting the data to the specified units.
#
# The function expects the following arguments:
#  - \em path: the file path to the output file, interpreted as described for the pts.utils.absPath() function.
#  - \em quantities: a sequence of astropy quantity arrays of the same length, one for each column
#    (if the output has only a single row, each item in sequence is a scalar astropy quantity).
#  - \em units: a string containing a comma-separated list of valid SKIRT unit strings.
#  - \em descriptions: a string containing a comma-separated list of quantity descriptions.
#  - \em title: an optional title string written as the first line of the file; default is no title
#  - \em fmt: an optional format string; the default format is "%1.9e".
#
# The arguments must specify the same number of quantities, units, and descriptions. The units and descriptions
# will appear in the header of the output file. The quantities are converted to the requested units before being
# written to the file.
#
def saveColumns(path, quantities, units, descriptions, *, title=None, fmt="%1.9e"):

    # split descriptions and units into segments
    descriptions = [ s.strip() for s in descriptions.split(",") ]
    units = [ s.strip() for s in units.split(",") ]

    # verify length of input sequences
    if len(quantities) != len(units) or len(quantities) != len(descriptions):
        raise ValueError("Number of units or descriptions does not match number of quantities")

    # convert the quantities to the requested units
    quantities = [quantity.to(smunit(unit)) for quantity, unit in zip(quantities, units)]

    # open the file
    path = ut.absPath(path)
    with open(path, 'wt') as outfile:

        # write the header
        if title:
            outfile.write("# {}\n".format(title))
        for col, (description, unit) in enumerate(zip(descriptions, units)):
            outfile.write("# column {}: {} ({})\n".format(col+1, description, unit))

        # write the data
        np.savetxt(outfile, np.stack(quantities).T, fmt=fmt)

# -----------------------------------------------------------------
