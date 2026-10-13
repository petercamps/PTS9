#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.hdf5.do.convert_from_hdf5 Convert bundles in an HDF5 file to SKIRT input and output files
#
# This script converts one or more bundles in an HDF5 file to the corresponding regular SKIRT input or output
# files, so that SKIRT, or a tool that does not read HDF5 files, can read them as regular files. It performs the
# opposite conversion of the \c convert_to_hdf5 script. The HDF5 file is not changed.
#
# The script supports all bundles defined by SKIRT: a text column file bundle becomes a column text file, an AMR
# text file bundle a column text file describing an adaptive mesh, a stored table bundle a SKIRT stored table file,
# a FITS file bundle a FITS file, a spatial grid plot file bundle a spatial grid plot file, and an unstructured text
# file bundle a text file. The type of each bundle is derived from the datasets and attributes it holds. See the
# pts.hdf5.convert module for more information, including a description of the contents of the resulting files.
#
# The script takes the following arguments:
#  - \em path (positional string argument): the bundle(s) to be converted, in one of the following forms:
#     - the input location of a SKIRT simulation, with the same semantics as the \c -i command line option of
#       SKIRT: a directory, followed by the name of an HDF5 file in that directory (with the \c .hdf5 filename
#       extension), optionally followed by a colon and the path of a suite within that file. All bundles in the
#       root group of the file or in that suite are converted, but not those in nested suites;
#     - the same, followed by a forward slash and the name of a bundle: just that bundle is converted;
#     - the same, with wildcards in the name of the bundle: all bundles in the root group of the file or in the
#       suite with a name matching the pattern are converted. The pattern must be enclosed in quotes, to prevent
#       the shell from expanding it.
#  - \em out (string): the output location: the existing directory in which the files are placed, as for the \c -o
#    command line option of SKIRT when it does not name an HDF5 file. If the argument is omitted, the files are
#    placed in the directory that holds the HDF5 file, which is where SKIRT looks for a regular input file before
#    it looks for the bundle.
#  - \em name (string): the name of the file, when converting a single bundle. By default, each file is named
#    for its bundle, which is the name under which SKIRT looks for the regular file.
#  - \em names (string): a comma-separated list of column descriptions for a text column file bundle or an AMR
#    text file bundle, replacing the descriptions derived from the names of the datasets, when converting a single
#    bundle.
#  - \em units (string): a comma-separated list of column units for a text column file bundle or an AMR text file
#    bundle, replacing the units recorded in the bundle, when converting a single bundle. This is intended for a
#    bundle that does not record the units of its columns, because the header of a column text file cannot
#    describe a column without specifying its unit.
#
# The script never replaces an existing file without confirmation. If one or more of the files to be written exist,
# the script lists them and asks to enter "YES" before it converts anything; with any other answer, no bundles
# are converted.
#
# For example, to convert all bundles in a suite of an HDF5 file to regular files in an existing directory:
#
#     pts convert_from_hdf5 in/data.hdf5:galaxy042 --out galaxy042
#
# after which SKIRT can be run with the command line option <tt>-i galaxy042</tt>. Or to obtain the FITS files
# for the output of a simulation with prefix "galaxy" that was written to an HDF5 file, next to that file:
#
#     pts convert_from_hdf5 "out/results.hdf5:galaxy_*.fits"
#
# Because a file in SKIRT stored columns format becomes a text column file bundle, the two scripts together
# convert such a file, which SKIRT 10 no longer supports, to a column text file:
#
#     pts convert_to_hdf5 in/stars.scol
#     pts convert_from_hdf5 in/stars.hdf5:stars.txt
#

# -----------------------------------------------------------------

def do( path : (str,"HDF5 file, suite, bundle, or quoted wildcard pattern, in the form <dir>/<hdf>:<suite>/<bundle>"),
        out : (str,"output directory, as for the SKIRT -o option without HDF5 file") = "",
        name : (str,"name of the file when converting a single bundle") = "",
        names : (str,"comma-separated column descriptions, replacing those in a single bundle") = "",
        units : (str,"comma-separated column units, replacing those in a single bundle") = "",
        ) -> "convert bundles in an HDF5 file to SKIRT input and output files":

    import logging
    import pts.hdf5 as hdf
    import pts.utils as ut

    # an empty string means that the option was not specified
    names = names if names else None
    units = units if units else None

    # only a single bundle can be given another file name and other column information
    bundlePaths = hdf.findBundles(path)
    single = len(bundlePaths) == 1
    if not single and (name or names or units):
        raise ut.UserError("The --name, --names, and --units options can be used only when converting " \
                           "a single bundle")

    # ask for confirmation before replacing existing files
    existing = hdf.filesToBeReplaced(path, out, fileName=name)
    if len(existing) > 0:
        listed = ", ".join([ filepath.name for filepath in existing[:10] ]) + (", ..." if len(existing) > 10 else "")
        logging.warning("** This will replace {} existing file(s) in {}: {} **" \
                        .format(len(existing), existing[0].parent, listed))
        try:
            response = input("Enter 'YES' to replace the file(s) ({}): ".format(len(existing)))
        except EOFError:
            response = ""
        if response != "YES":
            logging.info("No bundles were converted")
            return

    # for a single bundle, errors are reported as user errors
    if single:
        hdf.convertBundle(bundlePaths[0], out, fileName=name, names=names, units=units, overwrite=True)

    # for several bundles, a bundle that cannot be converted is reported, and the other bundles are still converted
    else:
        hdf.convertBundles(path, out, overwrite=True)

# -----------------------------------------------------------------
