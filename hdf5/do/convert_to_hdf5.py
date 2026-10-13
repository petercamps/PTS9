#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.hdf5.do.convert_to_hdf5 Convert SKIRT input and output files to bundles in an HDF5 file
#
# This script converts one or more regular SKIRT input or output files to the corresponding bundles in an HDF5
# file, so that SKIRT can read its input from that file, or so that the output of a simulation is collected in
# a single file. The source files are not changed or removed.
#
# The script supports column text files (including those describing an adaptive mesh), SKIRT stored columns files,
# SKIRT stored table files, FITS files, spatial grid plot files, log files, parameter files, and convergence
# reports. The type of each file is derived from its filename extension and, where needed, from the file names used
# by SKIRT or from the contents of the file. See the pts.hdf5.convert module for more information.
#
# The script takes the following arguments:
#  - \em path (positional string argument): the file(s) to be converted, in one of the following forms:
#     - the name or path of a single file;
#     - the path of a directory: all files in that directory with a recognized type are converted;
#     - a path with wildcards: all files matching the pattern are converted. The pattern must be enclosed in
#       quotes, to prevent the shell from expanding it.
#  - \em out (string): the output location, with the same semantics as the \c -o command line option of SKIRT:
#    a directory, optionally followed by the name of an HDF5 file in that directory (with the \c .hdf5 filename
#    extension), optionally followed by a colon and the path of a suite within that file. If the HDF5 file exists,
#    the bundles are added to it, replacing any bundles of the same name; otherwise the file is created. If the
#    argument does not name an HDF5 file, or is omitted, the bundles are written to an HDF5 file named for the
#    source file, or for the source directory when converting several files, with the \c .hdf5 filename extension.
#    That file is placed in the specified directory, or else next to the source file(s).
#  - \em name (string): the name of the bundle, when converting a single file. By default, each bundle is named
#    for its source file, including the filename extension (except for the \c .scol extension), so that SKIRT
#    finds it under the name already used in the ski file.
#  - \em type (string): the type of the file(s), overriding the type derived from the file name. The type is one of
#    \c columns, \c amr, \c storedcolumns, \c storedtable, \c fitsinput, \c fitsoutput, \c gridplot, or \c text.
#  - \em names (string): a comma-separated list of column descriptions for a column text file, replacing those in
#    the header of the file. This is intended for a file without column information.
#  - \em units (string): a comma-separated list of column units for a column text file, replacing those in
#    the header of the file. For a spatial grid plot file, the argument provides the unit of the coordinates,
#    in case the parameter file of the simulation, from which this unit is normally derived, is not in the same
#    directory.
#  - \em compress (integer): the level at which the datasets are compressed, from 1 (fastest) to 9 (smallest);
#    by default the datasets are not compressed.
#
# For example, to convert all input files of a simulation to bundles in a suite of a new or existing HDF5 file:
#
#     pts convert_to_hdf5 in --out in/data.hdf5:galaxy042
#
# after which SKIRT can be run with the command line option <tt>-i in/data.hdf5:galaxy042</tt>. Or to collect all
# output files of a simulation with prefix "galaxy" in a new or existing HDF5 file:
#
#     pts convert_to_hdf5 "out/galaxy_*" --out out/results.hdf5
#

# -----------------------------------------------------------------

def do( path : (str,"file, directory, or quoted wildcard pattern specifying the file(s) to be converted"),
        out : (str,"output location in the form <dir>/<hdf>:<suite>, as for the SKIRT -o option") = "",
        name : (str,"name of the bundle when converting a single file") = "",
        type : (str,"type of the file(s), overriding the type derived from the file name: columns, amr, "
                    "storedcolumns, storedtable, fitsinput, fitsoutput, gridplot, or text") = "",
        names : (str,"comma-separated column descriptions, replacing those in the file header") = "",
        units : (str,"comma-separated column units, replacing those in the file header") = "",
        compress : (int,"compression level for the datasets, from 1 to 9, or 0 for no compression") = 0,
        ) -> "convert SKIRT input and output files to bundles in an HDF5 file":

    import pts.hdf5 as hdf
    import pts.utils as ut

    # an empty string means that the option was not specified
    names = names if names else None
    units = units if units else None

    # a single file can be given another bundle name; errors are reported as user errors
    if ut.absPath(path).is_file():
        hdf.convertFile(path, out, bundleName=name, fileType=type, names=names, units=units, compression=compress)

    # for several files, a file that cannot be converted is reported, and the other files are still converted
    else:
        if name:
            raise ut.UserError("The --name option can be used only when converting a single file")
        hdf.convertFiles(path, out, fileType=type, names=names, units=units, compression=compress)

# -----------------------------------------------------------------
