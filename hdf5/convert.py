#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.hdf5.convert Converting SKIRT input and output files to bundles in an HDF5 file
#
# The functions in this module convert the regular files that SKIRT reads as input or writes as output to the
# corresponding bundle in an HDF5 file, as described in the pts.hdf5.bundle module. This allows, for example,
# migrating the input files of an existing simulation to an HDF5 file, or collecting the output files of a completed
# simulation in an HDF5 file. The source files are never changed or removed.
#
# File types
# ----------
#
# There is a conversion function for each type of file, listed in the table below. The convertFile() function
# converts a single file after deriving its type from the file name, and the convertFiles() function does the same
# for all files in a directory or matching a wildcard pattern.
#
# | Type | Function | Kind of file
# |------|----------|-------------
# | columns | convertTextColumnFile() | Column text file
# | amr | convertAdaptiveMeshFile() | Column text file describing an adaptive mesh
# | storedcolumns | convertStoredColumnsFile() | SKIRT stored columns file
# | storedtable | convertStoredTableFile() | SKIRT stored table file
# | fitsinput | convertFitsInputFile() | FITS file from which SKIRT reads just the pixel values
# | fitsoutput | convertFitsOutputFile() | FITS file written by SKIRT, with its axes and units
# | gridplot | convertGridPlotFile() | Text file with the data to plot a spatial grid
# | text | convertTextFile() | Text file stored as a whole (log file, parameter file, and some probe output)
#
# The two types of FITS file differ in what the bundle records next to the pixel values: nothing for the
# \c fitsinput type, and the axes, the units, and the observer information for the \c fitsoutput type. Both types
# of bundle hold the pixel values in the same way, so that SKIRT can read either of them as input.
#
# The type of a file is derived from its filename extension, ignoring case. Where the extension does not suffice
# (for example, \c .dat is used for several types), the type is derived from the file names that SKIRT uses for its
# output files, or from the contents of the file, as described for the inferFileType() function. A file of another
# type, or with an unusual name, can be converted by calling the appropriate function directly, or by specifying
# its type explicitly.
#
# Output location
# ---------------
#
# The conversion functions have an optional \em outPath argument that specifies where the bundle is written, with
# the same semantics as the \c -o command line option of SKIRT. Its value has the form
#
#     <dir>/<hdf>:<suite>
#
# where \c \<dir\> is an existing directory, \c \<hdf\> is the name of an HDF5 file in that directory, with the
# \c .hdf5 filename extension, and \c \<suite\> is the path of a suite within that file. The HDF5 file is created
# if it does not exist; otherwise the bundle is added to it, replacing an existing bundle of the same name. If the
# suite is omitted, the bundle is placed in the root group of the file. If the HDF5 file is omitted as well, or if
# the argument is omitted altogether, the bundle is written to an HDF5 file with the name of the source file and
# the \c .hdf5 filename extension, in the specified directory or next to the source file. The convertFiles()
# function uses the name of the source directory instead, so that all its bundles end up in the same file.
#
# Bundle name
# -----------
#
# By default, the bundle is named for the source file, including its filename extension. SKIRT then finds the
# bundle under the file name already used in the ski file, provided that the \c -i command line option specifies
# the HDF5 file, and the suite if there is one. The conversion functions have an optional \em bundleName argument
# to specify another name.
#
# \note A file in SKIRT stored columns format is the exception: the \c .scol filename extension is replaced by
# \c .txt in the name of the bundle, so that the file \c stars.scol becomes the bundle \c stars.txt. SKIRT 10 no
# longer supports the stored columns format, and rejects an input file name with the \c .scol extension. The
# bundle therefore gets the name of the column text file with the same contents, and the ski file must be adjusted
# to refer to that name, here \c stars.txt instead of \c stars.scol.
#
# If two source files would become a bundle with the same name, as for the files \c stars.scol and \c stars.txt in
# the same directory, the convertFiles() function converts neither of them and reports the conflict. Convert one
# of them, or convert them separately with different bundle names.
#

# -----------------------------------------------------------------

import glob
import logging
import os.path
import pathlib
import re
import struct
import warnings
import astropy.io.fits as fits
import numpy as np
import pts.storedtable as stab
import pts.utils as ut
from .bundle import resolveBundlePath, writeColumnBundle, writeAdaptiveMeshBundle, \
                    writeStoredTableBundle, writeFitsInputBundle, writeFitsOutputBundle, writeGridPlotBundle, \
                    writeTextBundle

# -----------------------------------------------------------------

# the supported file types, in the order of the table in the header of this module
_fileTypes = ("columns", "amr", "storedcolumns", "storedtable", "fitsinput", "fitsoutput", "gridplot", "text")

# the syntax of a header line with column information in a column text file, as defined by SKIRT; the groups are
# the optional one-based column index, the description, and the unit string
_columnInfoSyntax = re.compile(r"#\s*column\s*(\d*)\s*:\s*([^()]*)\(\s*([a-zA-Z0-9/]*)\s*\)\s*", re.IGNORECASE)

# the syntax of a floating point number in a column text file
_floatSyntax = re.compile(r"[+-]?((\d+\.?\d*|\.\d+)([eE][+-]?\d+)?|nan|inf|infinity)", re.IGNORECASE)

# the file names used by SKIRT for its log files (possibly with the rank of a process), its convergence reports,
# its grain population tables (with the index of a medium component), and its spatial grid plot files, in lowercase
_logFileName = re.compile(r"_log(p\d+)?\.txt$")
_convergenceFileName = re.compile(r"_convergence\.dat$")
_grainPopulationsFileName = re.compile(r"_grainpops_\d+\.dat$")
_gridPlotFileName = re.compile(r"_grid_(xy|xz|yz|xyz)\.dat$")

# the header cards written by SKIRT to describe the observer, with the corresponding bundle attribute
_observerCards = (("CROTA1", "inclination"), ("CROTA2", "azimuth"), ("CROTA3", "roll"), ("REDSHIFT", "redshift"),
                  ("DISTLUMI", "luminosity_distance"), ("DISTANGD", "angular_diameter_distance"),
                  ("DISTUNIT", "distance_unit"))

# the syntax of the element in a SKIRT parameter file that specifies the unit system of the simulation, and
# the unit of length for each of the unit systems
_unitSystemSyntax = re.compile(r"<units\b[^>]*>\s*<(\w+)")
_lengthUnits = dict(SIUnits="m", StellarUnits="AU", ExtragalacticUnits="pc")

# the number of text lines converted to numbers at a time when reading a file that cannot be read as a whole
_chunkSize = 100000

# -----------------------------------------------------------------

## This function returns a tuple with the names of the supported file types, as listed in the header of this module.
def fileTypes():
    return _fileTypes

## This function returns the type of the specified SKIRT input or output file, as one of the type names listed in
# the header of this module, or None if the type of the file is not recognized. The type is derived from the
# filename extension, ignoring case, as follows:
#
#  - \c .scol: a stored columns file; \c .stab: a stored table file.
#  - \c .xml: unstructured text.
#  - \c .fits: a FITS file written by SKIRT if the ORIGIN card in its primary header mentions SKIRT, which is the
#    case for the files written by SKIRT and by PTS; otherwise a FITS file that provides just pixel values.
#  - \c .txt or \c .dat: unstructured text if the file name ends in \c _log.txt, in \c _convergence.dat, or in
#    \c _grainpops_N.dat, where N is a number; data to plot a spatial grid if the name ends in \c _grid_xy.dat,
#    \c _grid_xz.dat, \c _grid_yz.dat, or \c _grid_xyz.dat; a column text file describing an adaptive mesh if the
#    first line that is not empty and not a comment starts with an exclamation mark; and a regular column text
#    file otherwise.
#
# The files with a name ending in \c _grainpops_N.dat are the grain population tables written by the
# DustGrainPopulationsProbe. Although such a table has columns, one of them holds the name of a grain material
# rather than a number, so that the table cannot be represented as a text column file bundle. It is stored as
# unstructured text instead.
#
# The file path is interpreted as described for the pts.utils.absPath() function. The file must exist if its
# contents are needed to determine the type.
def inferFileType(filePath):
    path = ut.absPath(filePath)
    name = path.name.lower()
    suffix = path.suffix.lower()

    if suffix == ".scol": return "storedcolumns"
    if suffix == ".stab": return "storedtable"
    if suffix == ".xml": return "text"
    if suffix == ".fits":
        return "fitsoutput" if _isWrittenBySkirt(path) else "fitsinput"
    if suffix == ".txt" or suffix == ".dat":
        if _logFileName.search(name) or _convergenceFileName.search(name): return "text"
        if _grainPopulationsFileName.search(name): return "text"
        if _gridPlotFileName.search(name): return "gridplot"
        return "amr" if _scanColumnFile(path).startsWithNonLeaf else "columns"
    return None

# -----------------------------------------------------------------

## This function converts the specified SKIRT input or output file to the corresponding bundle in an HDF5 file,
# and returns the bundle path. It derives the type of the file from its name, as described for the
# inferFileType() function, and calls the conversion function for that type.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em fileType: the type of the file, as one of the type names listed in the header of this module, overriding
#    the type derived from the file name.
#  - \em names, \em units: the column descriptions and units for a column text file; see convertTextColumnFile().
#    For data to plot a spatial grid, \em units specifies the unit of the coordinates.
#  - \em compression: the compression level for the datasets, from 1 to 9; by default there is no compression.
#
# The function raises a user error if the type of the file is not recognized, or if the file cannot be converted.
def convertFile(inFilePath, outPath=None, *, bundleName=None, fileType=None, names=None, units=None,
                compression=None):
    inpath = _existingFile(inFilePath)

    # determine the type of the file
    if not fileType:
        fileType = inferFileType(inpath)
        if fileType is None:
            raise ut.UserError("Cannot derive the type of file '{}' from its name; specify the type explicitly" \
                               .format(inpath))
    elif fileType not in _fileTypes:
        raise ut.UserError("Unsupported file type '{}'; the file type must be one of: {}" \
                           .format(fileType, ", ".join(_fileTypes)))

    # verify that the options are relevant for this type
    if names is not None and fileType not in ("columns", "amr"):
        raise ut.UserError("Column descriptions cannot be specified for a file of type '{}'".format(fileType))
    if units is not None and fileType not in ("columns", "amr", "gridplot"):
        raise ut.UserError("Units cannot be specified for a file of type '{}'".format(fileType))

    # call the appropriate conversion function
    if fileType == "columns":
        return convertTextColumnFile(inpath, outPath, bundleName=bundleName, names=names, units=units,
                                     compression=compression)
    if fileType == "amr":
        return convertAdaptiveMeshFile(inpath, outPath, bundleName=bundleName, names=names, units=units,
                                       compression=compression)
    if fileType == "storedcolumns":
        return convertStoredColumnsFile(inpath, outPath, bundleName=bundleName, compression=compression)
    if fileType == "storedtable":
        return convertStoredTableFile(inpath, outPath, bundleName=bundleName, compression=compression)
    if fileType == "fitsinput":
        return convertFitsInputFile(inpath, outPath, bundleName=bundleName, compression=compression)
    if fileType == "fitsoutput":
        return convertFitsOutputFile(inpath, outPath, bundleName=bundleName, compression=compression)
    if fileType == "gridplot":
        return convertGridPlotFile(inpath, outPath, bundleName=bundleName, unit=units, compression=compression)
    return convertTextFile(inpath, outPath, bundleName=bundleName)

# -----------------------------------------------------------------

## This function converts a set of SKIRT input or output files to bundles in a single HDF5 file, and returns a list
# with the bundle paths of the bundles that were written. It calls the convertFile() function for each file.
#
# The \em inPath argument specifies the files to be converted in one of the following ways:
#  - the path of a single file;
#  - the path of a directory: all files in that directory, but not in its subdirectories;
#  - a path that includes wildcards (as supported by the standard \c glob module): all files matching the pattern.
#
# Hidden files and HDF5 files are never converted. A file with a type that is not recognized is skipped, unless
# the \em fileType argument specifies the type of all files. If a file cannot be converted, the function logs an
# error message and continues with the next file. This includes the files that would become a bundle with the
# same name as another file, such as a stored columns file and a column text file with the same base name (see
# the header of this module): none of these is converted.
#
# The other arguments are the same as for the convertFile() function, except that there is no argument to specify
# the name of the bundles, and that the default HDF5 file is different. If the output location does not name an
# HDF5 file, the bundles for a single source file are written to an HDF5 file with the name of that file, and
# the bundles for a directory or a wildcard pattern to an HDF5 file with the name of the directory that holds the
# source files; in both cases with the \c .hdf5 filename extension, and in the directory specified by the output
# location or else in the directory that holds the source files.
def convertFiles(inPath, outPath=None, *, fileType=None, names=None, units=None, compression=None):
    # get the files to be converted and the default HDF5 file
    inpath = ut.absPath(inPath)
    if inpath.is_file():
        filepaths = [ inpath ]
        defaultpath = inpath.with_suffix(".hdf5")
    else:
        if inpath.is_dir():
            filepaths = inpath.iterdir()
        else:
            filepaths = [ pathlib.Path(path) for path in glob.glob(str(inpath)) ]
        filepaths = sorted([ path for path in filepaths if path.is_file() and not path.name.startswith(".") \
                                                           and path.suffix.lower() != ".hdf5" ])
        if len(filepaths) == 0:
            raise ut.UserError("There are no files to be converted for '{}'".format(inpath))
        dirpath = pathlib.Path(os.path.commonpath([ path.parent for path in filepaths ]))
        defaultpath = dirpath / ((dirpath.name or "skirt") + ".hdf5")

    # resolve the output location once, so that all bundles are written to the same HDF5 file and suite
    location = resolveBundlePath(outPath, "", defaultpath)

    # select the files with a recognized type
    selected = []
    numSkipped = 0
    numFailed = 0
    for filepath in filepaths:
        try:
            if fileType or inferFileType(filepath) is not None:
                selected.append(filepath)
            else:
                logging.info("Skipping file with unrecognized type: {}".format(filepath))
                numSkipped += 1
        except ut.UserError as error:
            logging.error("Cannot convert {}: {}".format(filepath, error.message))
            numFailed += 1
        except OSError as error:
            logging.error("Cannot convert {}: {}".format(filepath, error))
            numFailed += 1

    # convert each of the files, except those that would become a bundle with the same name as another file
    bundlePaths = []
    bundleNames = [ _defaultBundleName(filepath) for filepath in selected ]
    for filepath, bundleName in zip(selected, bundleNames):
        try:
            if bundleNames.count(bundleName) > 1:
                raise ut.UserError("{} files would become a bundle named '{}'; convert one of them, or convert " \
                                   "them separately with different bundle names" \
                                   .format(bundleNames.count(bundleName), bundleName))
            bundlePaths.append(convertFile(filepath, location, fileType=fileType, names=names, units=units,
                                           compression=compression))
        except ut.UserError as error:
            logging.error("Cannot convert {}: {}".format(filepath, error.message))
            numFailed += 1
        except OSError as error:
            logging.error("Cannot convert {}: {}".format(filepath, error))
            numFailed += 1

    # log a summary
    message = "Converted {} file(s) to bundles in {}".format(len(bundlePaths), location.rstrip(":"))
    if numSkipped > 0: message += "; skipped {} file(s)".format(numSkipped)
    if numFailed > 0: message += "; failed to convert {} file(s)".format(numFailed)
    logging.info(message)
    return bundlePaths

# -----------------------------------------------------------------

## This function converts a SKIRT column text file to a text column file bundle, and returns the bundle path.
# The file is interpreted in the same way as SKIRT does:
#
#  - Empty lines and lines starting with a hash character are ignored, except for the header lines with column
#    information. Such a line has the form "# column 1: description (unit)", where the column number is optional,
#    the description cannot contain parentheses, and the unit string can be empty. The column information lines
#    must precede the data and must be in the order of the columns.
#  - Each other line holds the values for one row, separated by white space. If the header contains column
#    information, the number of information lines determines the number of columns, and any further text on a
#    data line is ignored. If not, the number of columns is the number of values on the first data line.
#
# The bundle receives a dataset for each column, named for the column description and with the column's unit, as
# described for the pts.hdf5.bundle.writeColumnBundle() function. The comment lines at the start of the file, up to
# the first line with column information or with data, become the description of the bundle; SKIRT writes a
# single such line in its output files. For a file without column information, the datasets are named "column0",
# "column1", and so on, and they carry no unit, so that SKIRT uses its default units, just as it does for the
# text file.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em names: the column descriptions, as a sequence of strings or as a single string with a comma-separated
#    list, replacing the descriptions in the header of the file.
#  - \em units: the column units, in the same form, replacing the units in the header of the file.
#  - \em compression: the compression level for the datasets, from 1 to 9; by default there is no compression.
#
# If the file has no column information, the \em names or \em units argument also determines the number of columns.
def convertTextColumnFile(inFilePath, outPath=None, *, bundleName=None, names=None, units=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    # get the column information and load the data
    scan = _scanColumnFile(inpath)
    names, units = _columnInfo(inpath, scan, names, units)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    # numpy warns about a file without data
            rows = np.loadtxt(inpath, dtype=np.float64, comments="#", usecols=range(len(names)), ndmin=2,
                              encoding="latin-1")
    except ValueError as error:
        _verifyRows(inpath, len(names), False)
        raise ut.UserError("{} cannot be read as a column text file: {}".format(inpath, error)) from None

    writeColumnBundle(bundlePath, rows.T, names, units, description=scan.description, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a SKIRT column text file that describes an adaptive mesh to an AMR text file bundle,
# and returns the bundle path. The file lists the nodes of the mesh in the order described for the
# pts.hdf5.bundle.writeAdaptiveMeshBundle() function. A line starting with an exclamation mark represents a nonleaf
# node, and specifies the number of subdivisions of that node in each of the three spatial directions. Any other
# line represents a leaf node, and holds the values of its properties, as in a regular column text file. The
# header of the file, including any column information, is interpreted as described for the
# convertTextColumnFile() function, and the arguments of this function are the same as well.
def convertAdaptiveMeshFile(inFilePath, outPath=None, *, bundleName=None, names=None, units=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    # get the column information
    scan = _scanColumnFile(inpath)
    names, units = _columnInfo(inpath, scan, names, units)

    # read the nodes; the lines for the leaf nodes are converted to numbers in chunks to limit memory use
    isLeaf = []
    topology = []
    chunks = []
    lines = []
    try:
        with open(inpath, encoding="latin-1") as infile, warnings.catch_warnings():
            warnings.simplefilter("ignore")    # numpy warns about a chunk without data
            for line in infile:
                stripped = line.lstrip()
                if not stripped or stripped[0] == "#": continue
                if stripped[0] == "!":
                    isLeaf.append(False)
                    nx, ny, nz = [ int(token) for token in stripped[1:].split()[:3] ]
                    topology.append((nx, ny, nz))
                else:
                    isLeaf.append(True)
                    lines.append(line)
                    if len(lines) == _chunkSize:
                        chunks.append(np.loadtxt(lines, dtype=np.float64, usecols=range(len(names)), ndmin=2))
                        lines = []
            chunks.append(np.loadtxt(lines, dtype=np.float64, usecols=range(len(names)), ndmin=2))
    except ValueError as error:
        _verifyRows(inpath, len(names), True)
        raise ut.UserError("{} cannot be read as an adaptive mesh file: {}".format(inpath, error)) from None

    rows = np.concatenate(chunks)
    writeAdaptiveMeshBundle(bundlePath, np.array(isLeaf, dtype=bool), np.array(topology, dtype=np.int64),
                            rows.T, names, units, description=scan.description, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a file in SKIRT stored columns format, a binary alternative to a column text file that
# is no longer supported by SKIRT 10, to a text column file bundle, and returns the bundle path. The datasets are
# named for the column names in the file and carry the corresponding units.
#
# \note By default, the bundle is named for the file with the \c .scol filename extension replaced by \c .txt;
# for example, the file \c stars.scol becomes the bundle \c stars.txt. SKIRT 10 rejects an input file name with
# the \c .scol extension, so the ski file must be adjusted to refer to the new name. The function logs the name
# to be used.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the name of the file with the \c .scol filename extension
#    replaced by \c .txt, as described above.
#  - \em compression: the compression level for the datasets, from 1 to 9; by default there is no compression.
def convertStoredColumnsFile(inFilePath, outPath=None, *, bundleName=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    # read the file; refer to the pts.storedtable.io module for a description of the format
    try:
        with open(inpath, 'rb') as infile:
            if stab.io.stringFromFile(infile) != "SKIRT X" or stab.io.intFromFile(infile) != 0x010203040A0BFEFF \
                            or stab.io.intFromFile(infile) != 0:
                raise ut.UserError("File does not have SKIRT stored columns format: {}".format(inpath))
            numRows = stab.io.intFromFile(infile)
            numColumns = stab.io.intFromFile(infile)
            names = [ stab.io.stringFromFile(infile) for i in range(numColumns) ]
            units = [ stab.io.stringFromFile(infile) for i in range(numColumns) ]
            columns = stab.io.arrayFromFile(infile, (numColumns, numRows))
            if stab.io.stringFromFile(infile) != "SCOLEND":
                raise ut.UserError("File does not have the proper trailing tag: {}".format(inpath))
    except (ValueError, TypeError, MemoryError, OverflowError, struct.error):
        raise ut.UserError("File is not a complete SKIRT stored columns file: {}".format(inpath)) from None

    writeColumnBundle(bundlePath, columns, names, units, compression=compression)
    if not bundleName:
        logging.info("  The ski file must refer to '{}' instead of '{}'" \
                     .format(_defaultBundleName(inpath), inpath.name))
    return bundlePath

# -----------------------------------------------------------------

## This function converts a file in SKIRT stored table format to a stored table bundle, and returns the bundle path.
# The bundle receives a dataset for each axis and for each quantity in the table, as described for the
# pts.hdf5.bundle.writeStoredTableBundle() function.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em compression: the compression level for the datasets, from 1 to 9; by default there is no compression.
def convertStoredTableFile(inFilePath, outPath=None, *, bundleName=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    # read the file; refer to the pts.storedtable.io module for a description of the format
    try:
        with open(inpath, 'rb') as infile:
            if stab.io.stringFromFile(infile) != "SKIRT X" or stab.io.intFromFile(infile) != 0x010203040A0BFEFF:
                raise ut.UserError("File does not have SKIRT stored table format: {}".format(inpath))
            numAxes = stab.io.intFromFile(infile)
            if numAxes < 1 or numAxes > 9:
                raise ut.UserError("File does not have SKIRT stored table format: {}".format(inpath))
            axisNames = [ stab.io.stringFromFile(infile) for i in range(numAxes) ]
            axisUnits = [ stab.io.stringFromFile(infile) for i in range(numAxes) ]
            axisScales = [ stab.io.stringFromFile(infile) for i in range(numAxes) ]
            axisGrids = [ stab.io.arrayFromFile(infile, (stab.io.intFromFile(infile),)) for i in range(numAxes) ]
            numQuantities = stab.io.intFromFile(infile)
            quantityNames = [ stab.io.stringFromFile(infile) for i in range(numQuantities) ]
            quantityUnits = [ stab.io.stringFromFile(infile) for i in range(numQuantities) ]
            quantityScales = [ stab.io.stringFromFile(infile) for i in range(numQuantities) ]
            shape = tuple([ numQuantities ] + [ len(axisGrid) for axisGrid in axisGrids ])
            values = stab.io.arrayFromFile(infile, shape)
            if stab.io.stringFromFile(infile) != "STABEND":
                raise ut.UserError("File does not have the proper trailing tag: {}".format(inpath))
    except (ValueError, TypeError, MemoryError, OverflowError, struct.error):
        raise ut.UserError("File is not a complete SKIRT stored table file: {}".format(inpath)) from None

    writeStoredTableBundle(bundlePath, axisNames, axisUnits, axisScales, axisGrids,
                           quantityNames, quantityUnits, quantityScales, values, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a FITS file with an image or a data cube to a FITS file bundle for input, which holds
# just the pixel values, and returns the bundle path. This is the form in which SKIRT reads a FITS file as input,
# for example to define a geometry; none of the information in the header of the FITS file is used.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em hdu: the zero-based index or the name (EXTNAME) of the header-data unit that holds the image. By default,
#    the function uses the first unit that holds data, usually the primary unit, as SKIRT does. A bundle holds a
#    single image; use a separate bundle for each unit in a FITS file that is needed.
#  - \em compression: the compression level for the dataset, from 1 to 9; by default there is no compression.
#
# The image must have two or three dimensions; extra dimensions with a single element are removed.
def convertFitsInputFile(inFilePath, outPath=None, *, bundleName=None, hdu=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    with _openFits(inpath) as hdul:
        # locate the header-data unit
        if hdu is None:
            candidates = [ candidate for candidate in hdul if candidate.data is not None and candidate.data.size > 0 ]
            if len(candidates) == 0: raise ut.UserError("FITS file does not contain any data: {}".format(inpath))
            selected = candidates[0]
        else:
            try:
                selected = hdul[hdu]
            except (KeyError, IndexError):
                raise ut.UserError("FITS file {} has no header-data unit '{}'".format(inpath, hdu)) from None

        # get the image, removing extra dimensions with a single element
        if selected.data is None or selected.data.dtype.fields is not None:
            raise ut.UserError("Header-data unit in FITS file {} does not hold an image".format(inpath))
        image = np.asarray(selected.data, dtype=np.float64)
        if image.ndim > 3 and all([ size == 1 for size in image.shape[:-3] ]):
            image = image.reshape(image.shape[-3:])
        if image.ndim != 2 and image.ndim != 3:
            raise ut.UserError("Image in FITS file {} has {} dimension(s) rather than two or three" \
                               .format(inpath, image.ndim))

    writeFitsInputBundle(bundlePath, image, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a FITS file written by SKIRT to a FITS file bundle for output, which holds a frame or a
# data cube together with the grid points of its axes and with the units of the data and the axes, and returns
# the bundle path. SKIRT can read the resulting bundle as input as well, just as it can read the FITS file. The
# function obtains the information from the FITS file in the form written by SKIRT (and by
# the pts.simulation.fits.writeFits() function):
#
#  - the data from the primary header-data unit, with its unit from the BUNIT card;
#  - the grid points of each axis from the table extension for that axis ("X-axis coordinate values", and so on),
#    if present, or else from the CRPIX, CRVAL, and CDELT cards for a regular linear grid; the unit from the
#    CUNIT card;
#  - the inclination, azimuth, and roll angles of the observer, the redshift, and the distances, if present,
#    from the CROTA1, CROTA2, CROTA3, REDSHIFT, DISTLUMI, DISTANGD, and DISTUNIT cards.
#
# A FITS file written by SKIRT does not record the physical quantity represented by the data and by each axis,
# nor a description of its contents. These items can be provided through the corresponding arguments.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em quantities: a dictionary with the name of the physical quantity for "data", "x", "y" and/or "z".
#  - \em description: a free-form string that describes the contents of the bundle.
#  - \em compression: the compression level for the data, from 1 to 9; by default there is no compression.
def convertFitsOutputFile(inFilePath, outPath=None, *, bundleName=None, quantities=None, description=None,
                          compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    with _openFits(inpath) as hdul:
        # get the data
        header = hdul[0].header
        if hdul[0].data is None or hdul[0].data.ndim not in (2, 3):
            raise ut.UserError("Primary header-data unit in FITS file {} does not hold a frame or a data cube" \
                               .format(inpath))
        data = np.asarray(hdul[0].data, dtype=np.float64)
        dataUnit = _fitsString(header, "BUNIT")

        # get the axes
        x, xUnit = _fitsAxis(inpath, hdul, 1, data.shape[-1])
        y, yUnit = _fitsAxis(inpath, hdul, 2, data.shape[-2])
        z, zUnit = _fitsAxis(inpath, hdul, 3, data.shape[0]) if data.ndim == 3 else (None, None)

        # get the observer information
        observer = {}
        for card, key in _observerCards:
            if card in header: observer[key] = header[card]

    writeFitsOutputBundle(bundlePath, data, dataUnit, x, xUnit, y, yUnit, z, zUnit, quantities=quantities,
                          description=description, observer=observer, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a text file with the data to plot a spatial grid, as written by the SpatialGridPlotProbe
# in SKIRT, to a spatial grid plot file bundle, and returns the bundle path. Each line in the file holds the two
# or three coordinates of a point, separated by white space. Consecutive lines form a polyline; an empty line
# ends the polyline.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em unit: the unit string for the coordinates. The file itself does not record this unit, which is the unit
#    of length in the unit system of the simulation that wrote the file. If the argument is omitted, the function
#    derives the unit from the parameter file written by the same simulation (with a name ending in
#    \c _parameters.xml), if that file resides in the same directory. Otherwise the bundle does not record the
#    unit, and the function logs a warning.
#  - \em compression: the compression level for the datasets, from 1 to 9; by default there is no compression.
def convertGridPlotFile(inFilePath, outPath=None, *, bundleName=None, unit=None, compression=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)
    if unit is not None and not isinstance(unit, str):
        raise ut.UserError("The unit for the coordinates of a spatial grid plot must be a single string")

    # read the points; the lines are converted to numbers in chunks to limit memory use
    moveto = []
    chunks = []
    lines = []
    dimension = 0
    start = True
    try:
        with open(inpath, encoding="latin-1") as infile, warnings.catch_warnings():
            warnings.simplefilter("ignore")    # numpy warns about a chunk without data
            for line in infile:
                stripped = line.strip()
                if not stripped:
                    start = True
                elif stripped[0] != "#":
                    if dimension == 0: dimension = len(stripped.split())
                    moveto.append(start)
                    start = False
                    lines.append(line)
                    if len(lines) == _chunkSize:
                        chunks.append(np.loadtxt(lines, dtype=np.float64, ndmin=2))
                        lines = []
            chunks.append(np.loadtxt(lines, dtype=np.float64, ndmin=2))
        points = np.concatenate([ chunk for chunk in chunks if chunk.size > 0 ]) if len(moveto) > 0 else None
    except ValueError as error:
        raise ut.UserError("{} cannot be read as a spatial grid plot file: {}".format(inpath, error)) from None

    # a file without points still has a dimension, derived from its name
    if points is None:
        dimension = 3 if inpath.name.lower().endswith("_xyz.dat") else 2
        points = np.empty((0, dimension), dtype=np.float64)
    if dimension not in (2, 3) or points.shape[1] != dimension:
        raise ut.UserError("{} cannot be read as a spatial grid plot file: it has {} values on a line rather " \
                           "than two or three".format(inpath, dimension))

    if unit is None:
        unit = _lengthUnitFor(inpath)
        if unit is None:
            logging.warning("The unit of the coordinates in {} is not known; specify it to record it in the " \
                            "bundle".format(inpath))
    writeGridPlotBundle(bundlePath, points, np.array(moveto, dtype=bool), unit=unit, compression=compression)
    return bundlePath

# -----------------------------------------------------------------

## This function converts a text file to an unstructured text file bundle, which stores the contents of the file
# as a whole, without interpreting it, and returns the bundle path. This form serves for the log file, for the
# parameter file in XML format, for the convergence report written by the ConvergenceInfoProbe, and for the grain
# population table written by the DustGrainPopulationsProbe.
#
# The function accepts the following arguments:
#  - \em inFilePath: the path of the file to be converted, interpreted as described for pts.utils.absPath().
#  - \em outPath: the output location, as described in the header of this module.
#  - \em bundleName: the name of the bundle; by default the bundle is named for the file.
#  - \em textType: the type of text; one of "free form", "XML", or "log". By default, the type is "XML" for a file
#    with the \c .xml filename extension, "log" for a file with a name ending in \c _log.txt, and "free form"
#    for any other file.
#
# The file must be encoded in UTF-8, as all SKIRT output files are; a byte that cannot be decoded is replaced by the
# Unicode replacement character, and the function logs a warning.
def convertTextFile(inFilePath, outPath=None, *, bundleName=None, textType=None):
    inpath = _existingFile(inFilePath)
    bundlePath = _bundlePathFor(inpath, outPath, bundleName)

    # determine the type of text
    if textType is None:
        if inpath.suffix.lower() == ".xml": textType = "XML"
        elif _logFileName.search(inpath.name.lower()): textType = "log"
        else: textType = "free form"

    # read the text without changing its line endings
    contents = inpath.read_bytes()
    try:
        text = contents.decode("utf-8")
    except UnicodeDecodeError:
        text = contents.decode("utf-8", errors="replace")
        logging.warning("{} is not encoded in UTF-8; the offending characters have been replaced".format(inpath))
    if "\0" in text:
        raise ut.UserError("{} is not a text file".format(inpath))

    writeTextBundle(bundlePath, text, textType)
    return bundlePath

# -----------------------------------------------------------------

## This helper class holds the information obtained from the header of a column text file by the _scanColumnFile()
# function: the description assembled from the comment lines that precede the column information and the data (or
# None), the descriptions and units listed in the column information lines, the number of values on the first data
# line, and whether the first line that is not empty and not a comment specifies a nonleaf node of an adaptive mesh.
class _ColumnFileScan:
    def __init__(self):
        self.description = None
        self.titles = []
        self.units = []
        self.numValues = 0
        self.startsWithNonLeaf = False

## This helper function reads the header of the specified column text file, up to and including the first data
# line, and returns the information it finds as a _ColumnFileScan instance.
def _scanColumnFile(inpath):
    scan = _ColumnFileScan()
    comments = []
    first = True
    with open(inpath, encoding="utf-8", errors="replace") as infile:
        for line in infile:
            line = line.strip()
            if not line: continue

            # header line: column information or other comment
            if line[0] == "#":
                match = _columnInfoSyntax.fullmatch(line)
                if match:
                    scan.titles.append(" ".join(match.group(2).split()))
                    scan.units.append(match.group(3))
                    if match.group(1) and int(match.group(1)) != len(scan.titles):
                        raise ut.UserError("The header of {} specifies column {} where column {} is expected; a " \
                                           "line with column information is recognized only if it has the form " \
                                           "'# column N: description (unit)', where the unit consists of letters, " \
                                           "digits, and slashes".format(inpath, match.group(1), len(scan.titles)))
                elif len(scan.titles) == 0:
                    comment = line.lstrip("#").strip()
                    if comment: comments.append(comment)

            # line for a nonleaf node of an adaptive mesh
            elif line[0] == "!":
                if first: scan.startsWithNonLeaf = True
                first = False

            # first data line: count the values, i.e. the leading tokens formatted as a number
            else:
                for token in line.split():
                    if not _floatSyntax.fullmatch(token): break
                    scan.numValues += 1
                break

    if len(comments) > 0: scan.description = "\n".join(comments)
    return scan

## This helper function determines the descriptions and the units of the columns to be converted for the specified
# column text file, given the information scanned from its header and the descriptions and units specified by the
# caller, if any. It returns a tuple with two lists of the same length, the number of columns. An item in the list
# of units is None if the unit of the column is not known.
def _columnInfo(inpath, scan, names, units):
    names = _stringList(names)
    units = _stringList(units)

    # determine the number of columns
    if len(scan.titles) > 0: numColumns = len(scan.titles)
    elif names is not None: numColumns = len(names)
    elif units is not None: numColumns = len(units)
    else: numColumns = scan.numValues
    if numColumns == 0:
        raise ut.UserError("Cannot determine the number of columns in {}: it has no column information and " \
                           "no data".format(inpath))

    # verify the lists specified by the caller
    if names is not None and len(names) != numColumns:
        raise ut.UserError("{} column description(s) were specified for the {} column(s) in {}" \
                           .format(len(names), numColumns, inpath))
    if units is not None and len(units) != numColumns:
        raise ut.UserError("{} unit(s) were specified for the {} column(s) in {}" \
                           .format(len(units), numColumns, inpath))

    # use the information in the header for the items that were not specified by the caller
    if len(scan.titles) > 0:
        if names is None: names = scan.titles
        if units is None: units = scan.units
    else:
        if names is None and units is None:
            logging.warning("{} has no column information: the datasets get generic names and no units" \
                            .format(inpath))
        elif units is None:
            logging.warning("{} has no column information: the datasets get no units".format(inpath))
        if names is None: names = [""] * numColumns
        if units is None: units = [None] * numColumns
    return names, units

## This helper function verifies the data lines of a column text file that could not be loaded, and raises a user
# error that identifies the first offending line, if there is one. Each data line must start with the specified
# number of values. If the last argument is True, lines for the nonleaf nodes of an adaptive mesh are verified
# as well; if it is False, such lines are data lines.
def _verifyRows(inpath, numColumns, adaptiveMesh):
    with open(inpath, encoding="latin-1") as infile:
        for number, line in enumerate(infile, 1):
            line = line.strip()
            if not line or line[0] == "#": continue
            tokens = line.split()
            if adaptiveMesh and line[0] == "!":
                tokens = line[1:].split()
                if len(tokens) < 3 or not all([ re.fullmatch(r"[+-]?\d+", token) for token in tokens[:3] ]):
                    raise ut.UserError("Line {} of {} does not specify three subdivision counts" \
                                       .format(number, inpath))
                continue
            if len(tokens) < numColumns:
                raise ut.UserError("Line {} of {} has {} value(s) rather than {}" \
                                   .format(number, inpath, len(tokens), numColumns))
            for token in tokens[:numColumns]:
                if not _floatSyntax.fullmatch(token):
                    raise ut.UserError("'{}' on line {} of {} is not formatted as a floating point number" \
                                       .format(token, number, inpath))

## This helper function converts a list of strings specified by the caller, either as a sequence of strings or as
# a single string with comma-separated items, to a list of strings. It returns None for None or an empty string.
def _stringList(value):
    if value is None: return None
    if isinstance(value, str):
        if not value.strip(): return None
        return [ item.strip() for item in value.split(",") ]
    return [ str(item) for item in value ]

## This helper function returns the absolute path of the specified file, or raises a user error if it does not exist.
def _existingFile(filePath):
    path = ut.absPath(filePath)
    if not path.is_file():
        raise ut.UserError("File to be converted does not exist: '{}'".format(path))
    return path

## This helper function returns the default bundle name for the specified file: the name of the file, with the
# filename extension replaced by ".txt" for a stored columns file.
def _defaultBundleName(inpath):
    return inpath.with_suffix(".txt").name if inpath.suffix.lower() == ".scol" else inpath.name

## This helper function returns the bundle path for converting the specified file, given the output location and
# the bundle name specified by the caller, if any.
def _bundlePathFor(inpath, outPath, bundleName):
    return resolveBundlePath(outPath, bundleName if bundleName else _defaultBundleName(inpath),
                             inpath.with_suffix(".hdf5"))

## This helper function returns the unit of length used in the specified output file of a SKIRT simulation, or
# None if it cannot be determined. The unit follows from the unit system configured for the simulation, which
# the function obtains from the parameter file written by that simulation, if that file resides in the same
# directory. The parameter file is recognized by its name: the prefix of the simulation followed by
# "_parameters.xml", where the name of the specified file starts with the same prefix and an underscore.
def _lengthUnitFor(inpath):
    # locate the parameter file with the longest matching prefix
    prefixes = [ path.name[:-len("parameters.xml")] for path in inpath.parent.glob("*_parameters.xml") ]
    prefixes = [ prefix for prefix in prefixes if inpath.name.startswith(prefix) ]
    if len(prefixes) == 0: return None
    paramspath = inpath.parent / (max(prefixes, key=len) + "parameters.xml")

    # get the type of the unit system, which determines the unit of length
    try:
        match = _unitSystemSyntax.search(paramspath.read_text(encoding="utf-8", errors="replace"))
    except OSError:
        return None
    return _lengthUnits.get(match.group(1)) if match else None

## This helper function opens the specified FITS file and returns the list of its header-data units, or raises a
# user error if the file cannot be opened.
def _openFits(inpath):
    try:
        return fits.open(inpath)
    except OSError as error:
        raise ut.UserError("Cannot open FITS file {}: {}".format(inpath, error)) from None

## This helper function returns True if the primary header of the specified FITS file names SKIRT as its origin.
def _isWrittenBySkirt(inpath):
    with _openFits(inpath) as hdul:
        return "SKIRT" in str(hdul[0].header.get("ORIGIN", ""))

## This helper function returns the value of the specified card in the specified FITS header as a string without
# leading or trailing white space, or None if the header does not have the card.
def _fitsString(header, card):
    return str(header[card]).strip() if card in header else None

## This helper function returns a tuple with the grid points and the unit string for the axis with the specified
# one-based number in the specified FITS file written by SKIRT, which must have the specified number of grid
# points. An item in the tuple is None if the file does not provide it.
def _fitsAxis(inpath, hdul, number, size):
    header = hdul[0].header
    unit = _fitsString(header, "CUNIT{}".format(number))
    grid = None

    # the grid points are listed in a table extension, or they form a regular linear grid
    extension = "{}-axis coordinate values".format("XYZ"[number-1])
    cards = [ "{}{}".format(card, number) for card in ("CRPIX", "CRVAL", "CDELT") ]
    if extension in hdul:
        table = hdul[extension]
        try:
            grid = np.asarray(table.data["GRID_POINTS"], dtype=np.float64)
        except (KeyError, TypeError, ValueError):
            raise ut.UserError("The table with the grid points for axis {} in FITS file {} cannot be read" \
                               .format(number, inpath)) from None
        if unit is None: unit = _fitsString(table.header, "TUNIT1")
    elif all([ card in header for card in cards ]):
        crpix, crval, cdelt = [ float(header[card]) for card in cards ]
        grid = crval + cdelt * (np.arange(size) + 1 - crpix)

    if grid is not None and len(grid) != size:
        raise ut.UserError("Axis {} in FITS file {} has {} grid point(s) for {} pixel(s)" \
                           .format(number, inpath, len(grid), size))
    return grid, unit

# -----------------------------------------------------------------
