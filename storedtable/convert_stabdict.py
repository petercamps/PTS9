#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_stabdict Function to bundle stored table files into a stabdict archive
#
# The function in this module bundles a wildcarded set of files that already have the correct SKIRT
# stored table format into a single uncompressed tar archive, i.e. a "stabdict" file, without altering
# the contents of the individual files.

# -----------------------------------------------------------------

import tarfile
import pathlib
import glob

# -----------------------------------------------------------------

## This function expects a single input file path containing a "*" wildcard and a single output file path.
# It bundles all input files matching the wildcarded path into a single uncompressed tar archive at the
# output path, using the bare input filenames (i.e. without any directory prefix) as the archive member
# names. This is equivalent to invoking "COPYFILE_DISABLE=1 tar -cf outFilePath *.stab" from within the
# directory containing the input files, except that, because this function reads the file contents
# directly instead of shelling out to the platform tar utility, there is no risk of extraneous metadata
# (such as macOS AppleDouble "._filename" entries) ending up in the archive.
#
def convertToStabDict(inFilePaths, outFilePaths):
    with tarfile.open(outFilePaths[0], "w") as tar:
        for inFilePath in sorted(glob.glob(inFilePaths[0])):
            tar.add(inFilePath, arcname=pathlib.Path(inFilePath).name)

# -----------------------------------------------------------------
