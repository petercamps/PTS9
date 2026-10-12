#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.do.create_skirt_resources Create all or part of the SKIRT resources from original data
#
# This script creates all or part of the SKIRT resources (in SKIRT stored table format) from the original data
# residing in the \c SKIRT/Resources9/OriginalData directory hierarchy, according to the instructions embedded
# in that same hierarchy.
# See the pts.storedtable.conversionspec module for more information.
#
# Provide the name of a (possibly nested) subdirectory as the single argument to restrict execution to that
# portion of the directory tree, or the string "." to create all resources in the directory hierarchy.
#

# -----------------------------------------------------------------

def do( subDirectory : (str,"name of subdirectory to process or . for all"),
        ) -> "create SKIRT resources from original data":

    import pts.storedtable as stab
    import pts.utils as ut

    # get the paths to the top-level SKIRT resources input and output directories
    parentPath = ut.projectParentPath()
    inputPath = parentPath / "Resources9" / "OriginalData"
    outputPath = parentPath / "Resources9" / "StoredTables"

    # create the conversion suite
    specs = stab.createConversionSpecs(str(inputPath), str(outputPath), subDirectory if subDirectory!="." else None)

    # perform the conversion(s)
    specs.perform()

# -----------------------------------------------------------------
