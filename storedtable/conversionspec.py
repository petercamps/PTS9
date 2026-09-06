#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.conversionspec Facilities for representing SKIRT resource conversion specifications
#
# An instance of the ConversionSpec class in this module represents a particular specification for converting
# some original data (in some given format) to the SKIRT stored table format, and it provides facilities to
# actually perform the specified conversion by calling upon the appropriate functionality outside this module.
# A conversion specification is derived from the information in a particular subdirectory of the
# original data directory hierarchy.
#
# An instance of the ConversionSpecs class in this module aggregates a set of ConversionSpec instances with
# facilities to perform them as a group.
#
# The createConversionSpecs function creates the ConversionSpec instances corresponding to all specifications
# residing in a given nested directory hierarchy.

# -----------------------------------------------------------------

import logging
import os
import pts.utils as ut

# import the modules that implement the conversions
from .convert_band import *
from .convert_copy import *
from .convert_enthalpies import *
from .convert_kernel import *
from .convert_opticalprops import *
from .convert_sed import *
from .convert_singlegrain import *
from .convert_stabdict import *

# -----------------------------------------------------------------
#  createConversionSpecs function
# -----------------------------------------------------------------

## This function creates and returns the ConversionSpec instances corresponding to all specifications
# residing in the nested directory hierarchy below the specified sub-directory within the overall
# hierarchy rooted at the given input path. The output root path is passed on to each conversion spec.
def createConversionSpecs(inputPath, outputPath, subDirectory):

    # expand the directory paths
    inputPath = os.path.realpath(os.path.expanduser(inputPath))
    outputPath = os.path.realpath(os.path.expanduser(outputPath))
    subDirPath = findSubDirectory(inputPath, subDirectory)
    setName = os.path.basename(subDirPath)

    # recursively look for and create conversion specs in the nested hierarchy
    specs = ConversionSpecs(setName)
    for dirpath, dirnames, filenames in os.walk(subDirPath):
        if "ConversionSpec.txt" in filenames:
            spec = ConversionSpec(dirpath, outputPath)
            specs.add(spec)

    # return the aggregated set of specs
    return specs

## This helper function locates the specified sub-directory in the nested hierarchy specified by the root path.
def findSubDirectory(rootPath, subDirectory):
    # if the name of the subdirectory is empty, simply return the root directory
    if not subDirectory:
        return rootPath

    # recursively walk all subdirectories of the parent directory;
    # the path of the first subdirectory with a name matching the specified name is returned
    for dirpath, dirnames, filenames in os.walk(rootPath):
        for dirname in dirnames:
            # compare the directory name itself, and the directory name prefixed with the name of its parent directory,
            # to the target name (to allow search strings like "SED/Sun")
            if (dirname.lower() == subDirectory.lower()) or  \
               (os.path.join(os.path.basename(dirpath), dirname).lower() == subDirectory.lower()):
                   return os.path.join(dirpath, dirname)

    # if no match is found, raise an error
    raise ut.UserError("{} not found in {}".format(subDirectory, rootPath))


# -----------------------------------------------------------------
#  ConversionSpecs class
# -----------------------------------------------------------------

## An instance of this class aggregates a set of ConversionSpec instances with
# facilities to perform them as a group.
class ConversionSpecs:

    ## The constructor remembers the given name for the set and initializes an empty list of conversion specs
    def __init__(self, name):
        self._name = name
        self._specs = [ ]

    ## This function adds the given conversion spec to the list aggregated by this instance
    def add(self, spec):
        self._specs.append(spec)

    ## This function performs all of the conversion specs currently in the list
    def perform(self):
        logging.info("Starting set {} of {} conversion specs".format(self._name, len(self._specs)))
        for spec in self._specs:
            spec.perform()
        logging.info("Finished set of conversion specs")


# -----------------------------------------------------------------
#  ConversionSpec class
# -----------------------------------------------------------------

## An instance of this class represents a particular specification for converting
# some original data (in some given format) to the SKIRT stored table format, and it provides facilities to
# actually perform the specified conversion by calling upon the appropriate functionality outside this module.
#
# The conversion specification is read from the "ConversionSpec.txt" text file in a particular subdirectory
# of the SKIRTresources directory hierarchy. The specification file must contain one or more sequences of
# exactly three lines:
#  - line 1: the name of the conversion function to be called (function must be imported into this module)
#  - line 2: the input file path(s), relative to the location of the specification file, separated by spaces
#  - line 3: the output file path(s), relative to the global root of the output hierarchy, separated by spaces
#
# Because the filename arguments are separated by spaces, individual file names cannot contain spaces.
#
class ConversionSpec:

    ## The constructor remembers the specified input and output paths. The input path specifies the directory in
    # which the "ConversionSpec.txt" file resides. The output path specifies the global root of the output hierarchy.
    def __init__(self, inputPath, outputPath):
        self._name = os.path.basename(inputPath)
        self._inputPath = inputPath
        self._outputPath = outputPath

    ## This function performs the conversion specified for this instance.
    def perform(self):
        logging.info("Performing conversion spec {}".format(self._name))

        # open the specification file
        with open(os.path.join(self._inputPath, "ConversionSpec.txt")) as specfile:
            # read specifications until the end of the file has been reached
            while True:
                functionspec = specfile.readline().split()
                if len(functionspec) != 1: break
                function = functionspec[0]
                inputpaths = specfile.readline().split()
                outputpaths = specfile.readline().split()

                # add the appropriate prefixes to the input and output paths
                inputpaths = [ os.path.join(self._inputPath, path) for path in inputpaths ]
                outputpaths = [ os.path.join(self._outputPath, path) for path in outputpaths ]

                # create the output directories if needed
                for path in outputpaths:
                    dirpath = os.path.dirname(path)
                    if not os.path.exists(dirpath):
                        os.makedirs(dirpath)

                # perform the conversion (using any of the functions imported in the header of this module)
                globals()[function](inputpaths, outputpaths)


# -----------------------------------------------------------------
