#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.storedtable.convert_band Functions for converting broadband transmission curves to stored table format
#
# The functions in this module convert broadband response or transmission curves from their original formats
# to normalized transmission curves in SKIRT stored table format.

# -----------------------------------------------------------------

import glob
import pathlib
import xml.etree.ElementTree
import numpy as np
import astropy.constants as const
from .io import writeStoredTable

# -----------------------------------------------------------------

# NumPy 2.0 renamed numpy.trapz() to numpy.trapezoid(), and later versions no longer offer the old name
_trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# common private helper function to normalize a transmission curve to unity and output it as a stored table
def _writeNormalized(outFilePath, wavelengths, transmissions):
    transmissions /= _trapezoid(x=wavelengths, y=transmissions)
    writeStoredTable(outFilePath, ['lambda'], ['m'], ['lin'], [wavelengths],
                                  ['T'], ['1/m'], ['lin'], [transmissions])

# -----------------------------------------------------------------

## This function converts a set of filter response curves given in SVO XML format to stored table files.
# The function expects a single input file path template with wildcards expanding to the set of input file names,
# and a single output file path template in which the asterisk will be replaced by the appropriate band name.
# This band name is derived from the input file name without the filename extension by
# removing suffixes "mu" or "band", replacing dots by underscores, and converting to uppercase.
def convertResponseCurveSVO(inFilePaths, outFilePaths):

    # helper function to remove a suffix from a string
    def rchop(s, suffix):
        if suffix and s.endswith(suffix):
            return s[:-len(suffix)]
        return s

    # loop over the list of input files
    for inFilePath in glob.glob(inFilePaths[0]):

        # load the XML tree
        with open(inFilePath, 'r') as bandfile: bandtree = xml.etree.ElementTree.parse(bandfile)

        # get the response curve, converting wavelengths from Angstrom to meter
        values = np.array([ float(el.text) for el in bandtree.findall(".//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD") ])
        if len(values) < 4: raise ValueError("Transmission table not found in '{}'".format(inFilePath))
        wavelengths, transmissions = np.reshape(values, (-1, 2)).T
        wavelengths *= 1e-10

        # convert from response curve to transmission curve, with arbitrary scaling
        transmissions *= wavelengths

        # construct the output file path based on the input file name
        name = pathlib.Path(inFilePath).stem
        name = rchop(name, "mu")
        name = rchop(name, "band")
        name = name.replace(".", "_").upper()
        outFilePath = outFilePaths[0].replace("*",name)

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts filter transmission curves given in SVO XML format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
def convertTransmissionCurveSVO(inFilePaths, outFilePaths):

    # loop over the list of input/ouput filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # load the XML tree
        with open(inFilePath, 'r') as bandfile: bandtree = xml.etree.ElementTree.parse(bandfile)

        # get the transmission curve, converting wavelengths from Angstrom to meter
        values = np.array([ float(el.text) for el in bandtree.findall(".//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD") ])
        if len(values) < 4: raise ValueError("Transmission table not found in '{}'".format(inFilePath))
        wavelengths, transmissions = np.reshape(values, (-1, 2)).T
        wavelengths *= 1e-10

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts filter transmission curves given in JCMT text column format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
def convertTransmissionCurveJCMT(inFilePaths, outFilePaths):

    # loop over the list of input/output filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # get the transmission curve, converting wavelengths from micron to meter
        wavelengths, transmissions = np.loadtxt(inFilePath, unpack=True)
        wavelengths *= 1e-6

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts filter transmission curves given in Planck HFI text column format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
def convertTransmissionCurvePlanckHFI(inFilePaths, outFilePaths):

    # loop over the list of input/output filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # load the data
        wavenumbers, transmissions, uncertainties = np.loadtxt(inFilePath, skiprows=3, usecols=(0, 1, 2), unpack=True)

        # only keep the rows where the transmission value is higher than 10 x the uncertainty
        mask = transmissions > 10. * uncertainties
        wavenumbers = wavenumbers[mask]
        transmissions = transmissions[mask]

        # only keep the rows where the transmission is above 1/5000 of the peak transmission
        peak = np.max(transmissions)
        mask = transmissions > peak / 5000.
        wavenumbers = wavenumbers[mask]
        transmissions = transmissions[mask]

        # convert to wavelengths
        wavenumbers *= 1e2  # from 1/cm to 1/m
        wavelengths = 1. / wavenumbers

        # reverse order
        wavelengths = np.flipud(wavelengths)
        transmissions = np.flipud(transmissions)

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts filter transmission curves given in Planck LFI text column format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
def convertTransmissionCurvePlanckLFI(inFilePaths, outFilePaths):

    # loop over the list of input/output filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # load the data
        frequencies, transmissions = np.loadtxt(inFilePath, skiprows=1, usecols=(0, 1), unpack=True)

        # convert to wavelengths
        frequencies *= 1e9  # from GHz to Hz
        wavelengths = const.c.si.value / frequencies

        # reverse order
        wavelengths = np.flipud(wavelengths)
        transmissions = np.flipud(transmissions)

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts response curves given in Euclid/Rubin text column format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
def convertResponseCurveEuclidRubin(inFilePaths, outFilePaths):

    # loop over the list of input/output filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # load text columns and convert from Angstrom to meter
        wavelengths, transmissions = np.loadtxt(inFilePath, unpack=True)
        wavelengths *= 1e-10

        # remove leading and trailing values that are close to zero
        zeroT = transmissions.max() / 500.
        first = max(0, np.argmax(transmissions > zeroT) - 1)
        last = len(transmissions) - np.argmax(transmissions[::-1] > zeroT) + 1
        wavelengths = wavelengths[first:last]
        transmissions = transmissions[first:last]

        # convert from response curve to transmission curve, with arbitrary scaling
        transmissions *= wavelengths

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------

## This function converts transmission curves given in ALMA text column format to stored table files.
# The function expects one or more input file paths and the same number of corresponding output file paths.
# The appropriate wavelength range for each ALMA band is derived from the band number in the output file name.
def convertTransmissionCurveALMA(inFilePaths, outFilePaths):

    # define wavelength ranges, in micron, indexed on ALMA band number
    # (band 0 does not exist and bands 1 and 2 are not supported)
    wavelengthRanges = [ None, None, None,
       (2590, 3570), (1840, 2400), (1420, 1840), (1090, 1420), (800, 1090), (600, 780), (420, 500), (320, 380) ]

    # loop over the list of input/output filenames
    for inFilePath, outFilePath in zip(inFilePaths, outFilePaths):

        # get the band and thus the wavelength range from the filename
        band = int(outFilePath.split("_")[-2])
        wmin, wmax = wavelengthRanges[band]
        wmin *= 1e-6  # from micron to m
        wmax *= 1e-6

        # load text column data
        frequencies, transmissions = np.loadtxt(inFilePath, usecols=(0, 1), unpack=True)
        frequencies *= 1e9  # from GHz to Hz
        wavelengths = const.c.si.value / frequencies

        # only keep the rows inside the wavelength range for this band
        mask = (wmin < wavelengths) * (wavelengths < wmax)
        wavelengths = wavelengths[mask]
        transmissions = transmissions[mask]

        # reverse order
        wavelengths = np.flipud(wavelengths)
        transmissions = np.flipud(transmissions)

        # write normalized stored table
        _writeNormalized(outFilePath, wavelengths, transmissions)

# -----------------------------------------------------------------
