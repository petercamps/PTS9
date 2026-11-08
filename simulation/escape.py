#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.escape Loading the output of escape instruments
#
# This module offers functions for loading the output of a SKIRT EscapeInstrument, which records the luminosity
# escaping from the model in bins of escape direction. For each bin, the instrument outputs the flux density at the
# distance of the instrument, divided by the sky fraction of the bin, so that an isotropic source yields the same
# flux in every bin as an SEDInstrument at the same distance. Depending on the output style of the instrument, the
# fluxes are written to an SED file for each bin ("*_bin<N>_sed.dat", with the statistics in "*_bin<N>_sedstats.dat")
# or to a data cube for each flux component, with the bins along the spatial axes ("*_total.fits" and so on, with the
# statistics in "*_stats<k>.fits"). In both cases, the instrument also writes a bin list ("*_bins.dat"), which lists
# for each bin its index, its sky fraction, and the directions it covers. The functions in this module hide these
# differences: they load the fluxes of all bins of an instrument as a single array, whatever the output style.
#
# The bin list describes the binning completely:
#  - AllDirectionsEscapeBinning: a single bin, without further columns.
#  - InclinationEscapeBinning: the minimum and maximum inclination of each bin, plus the minimum and maximum
#    azimuth if there is more than one azimuth bin. Without folding, the inclination bins cover 0° to 180°; with
#    folding, they cover 0° to 90°, and each bin also includes the mirror directions in the other hemisphere.
#    The bins are numbered with the azimuth index running fastest.
#  - HEALPixEscapeBinning: the ring index, the pixel-in-ring index, and the center direction of each HEALPix pixel,
#    with the pixels numbered in the ring scheme.
#
# In a data cube, a bin of a HEALPixEscapeBinning occupies the pixel with x equal to its pixel-in-ring index and y
# equal to its ring index, as in the output of the HEALPixSkyInstrument; for the other binnings, the bins form a
# single row, with x equal to the bin index.
#
# The functions accept the output of an escape instrument as regular files or as bundles in an HDF5 file.
#

# -----------------------------------------------------------------

import re
import astropy.units as u
import numpy as np
from .fits import loadFits, getFitsAxes
from .text import getColumnDescriptions, loadColumns

# -----------------------------------------------------------------

## An instance of the EscapeBins class represents the bins of an escape instrument, as described by the bin list
# written by the instrument. The constructor loads the bin list; the other functions return information about
# the bins, in the order of the bin index.
class EscapeBins:

    ## The constructor accepts an Instrument instance representing an EscapeInstrument in a simulation, and loads
    # its bin list. If the instrument has no bin list, or the bin list has an unexpected format, the constructor
    # raises an error.
    def __init__(self, instrument):
        if instrument.type() != "EscapeInstrument":
            raise ValueError("Instrument '{}' is not an EscapeInstrument".format(instrument.name()))
        paths = instrument.outFilePaths("bins.dat")
        if len(paths) != 1:
            raise ValueError("Cannot find the bin list of escape instrument '{}'".format(instrument.name()))

        # load the columns by their description
        descriptions = [ description.split(";")[0].strip() for description in getColumnDescriptions(paths[0]) ]
        columns = [ np.atleast_1d(column) for column in loadColumns(paths[0]) ]
        self._columns = dict(zip(descriptions, columns))
        if "bin index" not in self._columns or "sky fraction" not in self._columns:
            raise ValueError("The bin list of escape instrument '{}' has an unexpected format".format(
                             instrument.name()))
        self._index = self._columns["bin index"].value.astype(int)
        if not np.array_equal(self._index, np.arange(len(self._index))):
            raise ValueError("The bin list of escape instrument '{}' does not list the bins in order".format(
                             instrument.name()))

        # determine the type of binning
        if "ring index" in self._columns:
            self._type = "HEALPixEscapeBinning"
        elif "minimum inclination" in self._columns:
            self._type = "InclinationEscapeBinning"
        else:
            self._type = "AllDirectionsEscapeBinning"

    ## This function returns the type of the escape binning (the corresponding SKIRT class name) as a string:
    # "AllDirectionsEscapeBinning", "InclinationEscapeBinning", or "HEALPixEscapeBinning".
    def binningType(self):
        return self._type

    ## This function returns the number of bins.
    def numBins(self):
        return len(self._index)

    ## This function returns a numpy array with the sky fraction of each bin, i.e. the fraction of the full sphere
    # of directions covered by the bin.
    def skyFractions(self):
        return self._columns["sky fraction"].value

    ## For an inclination binning, this function returns True if the hemispheres are folded, i.e. if each bin also
    # includes the mirror directions in the other hemisphere, and False otherwise. For other binnings, it returns
    # False.
    def isFolded(self):
        if self._type != "InclinationEscapeBinning": return False
        return self._columns["maximum inclination"].max().to_value(u.deg) < 179.

    ## For an inclination binning, this function returns a tuple with the number of inclination bins and the number
    # of azimuth bins; for other binnings, it raises an error.
    def inclinationAzimuthShape(self):
        if self._type != "InclinationEscapeBinning":
            raise ValueError("Escape binning is not an inclination binning")
        numAzimuth = len(np.unique(self._columns["minimum azimuth"])) if "minimum azimuth" in self._columns else 1
        return self.numBins() // numAzimuth, numAzimuth

    ## For an inclination binning, this function returns a tuple with two astropy quantity arrays: the minimum and
    # maximum inclination of each bin. For other binnings, it raises an error.
    def inclinationRanges(self):
        if self._type != "InclinationEscapeBinning":
            raise ValueError("Escape binning is not an inclination binning")
        return self._columns["minimum inclination"], self._columns["maximum inclination"]

    ## For an inclination binning with more than one azimuth bin, this function returns a tuple with two astropy
    # quantity arrays: the minimum and maximum azimuth of each bin. For an inclination binning with a single azimuth
    # bin, both arrays hold 0° and 360°. For other binnings, the function raises an error.
    def azimuthRanges(self):
        if self._type != "InclinationEscapeBinning":
            raise ValueError("Escape binning is not an inclination binning")
        if "minimum azimuth" not in self._columns:
            return np.zeros(self.numBins()) << u.deg, np.full(self.numBins(), 360.) << u.deg
        return self._columns["minimum azimuth"], self._columns["maximum azimuth"]

    ## For a HEALPix binning, this function returns its order, i.e. the base-2 logarithm of the number of pixels
    # along the side of each base pixel. For other binnings, it raises an error.
    def healpixOrder(self):
        if self._type != "HEALPixEscapeBinning":
            raise ValueError("Escape binning is not a HEALPix binning")
        return int(round(np.log2(np.sqrt(self.numBins() / 12))))

    ## For a HEALPix binning, this function returns a tuple with two integer numpy arrays: the ring index and the
    # pixel-in-ring index of each bin. For other binnings, it raises an error.
    def healpixPixels(self):
        if self._type != "HEALPixEscapeBinning":
            raise ValueError("Escape binning is not a HEALPix binning")
        return self._columns["ring index"].value.astype(int), self._columns["pixel-in-ring index"].value.astype(int)

    ## For a HEALPix binning, this function returns a tuple with two astropy quantity arrays: the inclination and
    # azimuth of the center of each bin. For other binnings, it raises an error.
    def centerDirections(self):
        if self._type != "HEALPixEscapeBinning":
            raise ValueError("Escape binning is not a HEALPix binning")
        return self._columns["center inclination"], self._columns["center azimuth"]

    ## This function returns a tuple with two integer numpy arrays: the x and y indices of the pixel occupied by
    # each bin in the data cubes written by the instrument with the DataCubes output style.
    def framePixels(self):
        if self._type == "HEALPixEscapeBinning":
            rings, pixels = self.healpixPixels()
            return pixels, rings
        return self._index.copy(), np.zeros(self.numBins(), dtype=int)

    ## This function returns a list with a short text describing the directions covered by each bin, for use in
    # plot legends, e.g. "all directions", "incl 45.57°–53.13°", "incl 0°–36.87°, azi 90°–180°" (with "folded"
    # added for folded hemispheres), or "ring 1 pixel 0 at incl 48.19°, azi 22.5°".
    def labels(self):
        if self._type == "AllDirectionsEscapeBinning":
            return [ "all directions" ]
        if self._type == "HEALPixEscapeBinning":
            rings, pixels = self.healpixPixels()
            theta, phi = self.centerDirections()
            return [ "ring {} pixel {} at incl {}, azi {}".format(j, i, _deg(t), _deg(p))
                     for j, i, t, p in zip(rings, pixels, theta, phi) ]
        tmin, tmax = self.inclinationRanges()
        labels = [ "incl {}–{}".format(_deg(a), _deg(b)) for a, b in zip(tmin, tmax) ]
        if "minimum azimuth" in self._columns:
            pmin, pmax = self.azimuthRanges()
            labels = [ "{}, azi {}–{}".format(label, _deg(a), _deg(b)) for label, a, b in zip(labels, pmin, pmax) ]
        if self.isFolded():
            labels = [ label + " folded" for label in labels ]
        return labels

## This helper function returns a string with the specified angle in degrees with up to four significant digits.
def _deg(angle):
    return "{:.4g}°".format(angle.to_value(u.deg))

# -----------------------------------------------------------------

## This function returns an EscapeBins instance representing the bins of the specified escape instrument,
# as described for the constructor of that class.
def loadEscapeBins(instrument):
    return EscapeBins(instrument)

## This function returns the zero-based bin index of the escape bin for the specified SED or SED statistics file
# path written by an escape instrument, i.e. the number following "_bin" in a file name ending with "_sed.dat" or
# "_sedstats.dat". If the file name does not have this form, the function returns None.
def escapeBinIndex(path):
    match = re.search(r"_bin(\d+)_sed(stats)?\.dat$", str(path))
    return int(match.group(1)) if match else None

## This function loads the fluxes of all bins of the specified escape instrument for the specified flux component,
# from the SED files or the data cubes written by the instrument, whichever are present. It returns a tuple with two
# astropy quantity arrays: the characteristic wavelengths of the instrument's wavelength grid (or frequencies or
# energies, for the corresponding output style, in the order of the output files), and the fluxes with shape
# (number of bins, number of wavelengths), in the flavor and units of the output files.
#
# The \em component argument specifies the flux component by the name used for the data cube files: "total",
# "transparent", "primarydirect", "primaryscattered", "secondarytransparent", "secondarydirect",
# "secondaryscattered", or "primaryscatteredN" with N the number of scattering events. If the instrument did not
# record the requested component, the function raises an error.
def loadEscapeFluxes(instrument, component="total"):
    bins = EscapeBins(instrument)

    # load from the SED files, if present
    sedPaths = _binFilePaths(instrument, bins, "sed.dat")
    if sedPaths is not None:
        description = _columnDescription(component)
        waves = None
        fluxes = []
        for path in sedPaths:
            descriptions = [ description.split(";")[0].strip() for description in getColumnDescriptions(path) ]
            if description not in descriptions:
                raise ValueError("Escape instrument '{}' did not record the {} flux".format(instrument.name(),
                                                                                         component))
            wave, flux = loadColumns(path, (0, descriptions.index(description)))
            if waves is None: waves = np.atleast_1d(wave)
            fluxes.append(np.atleast_1d(flux.to(fluxes[0].unit if fluxes else flux.unit)))
        return waves, np.stack([ flux.value for flux in fluxes ]) << fluxes[0].unit

    # load from the data cube
    path = _cubeFilePath(instrument, component)
    x, y = bins.framePixels()
    return getFitsAxes(path)[2], loadFits(path)[x, y, :]

## This function loads the relative statistical error of the flux in each bin and at each wavelength for the
# specified escape instrument, from the SED statistics files or the statistics data cubes written by the instrument,
# whichever are present. If the instrument did not record statistics, the function raises an error.
#
# The relative error is calculated as \f$R = \sqrt{\sum w_i^2/(\sum w_i)^2 - 1/N}\f$, where \f$w_i\f$ is the
# contribution of the i-th photon packet and \f$N\f$ the number of photon packets, if the \em numPackets argument
# is given, and as \f$R = \sqrt{\sum w_i^2}/\sum w_i\f$ otherwise, which is slightly larger. The function returns a
# tuple with two astropy quantity arrays: the characteristic wavelengths of the instrument's wavelength grid (as for
# the loadEscapeFluxes() function) and the dimensionless relative errors with shape (number of bins, number of
# wavelengths). The relative error is NaN where no photon packet contributed.
def loadEscapeRelativeErrors(instrument, numPackets=None):
    bins = EscapeBins(instrument)

    # load the sums of the first and second powers of the contributions
    statPaths = _binFilePaths(instrument, bins, "sedstats.dat")
    if statPaths is not None:
        columns = [ loadColumns(path, (0, 2, 3)) for path in statPaths ]
        waves = np.atleast_1d(columns[0][0])
        w1 = np.stack([ np.atleast_1d(column[1]).value for column in columns ])
        w2 = np.stack([ np.atleast_1d(column[2]).value for column in columns ])
    else:
        x, y = bins.framePixels()
        path1 = _cubeFilePath(instrument, "stats1")
        path2 = _cubeFilePath(instrument, "stats2")
        waves = getFitsAxes(path1)[2]
        w1 = loadFits(path1).value[x, y, :]
        w2 = loadFits(path2).value[x, y, :]

    # calculate the relative error
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = w2 / w1**2
        if numPackets is not None: variance -= 1 / numPackets
        errors = np.where(w1 > 0, np.sqrt(np.maximum(variance, 0)), np.nan)
    return waves, errors << u.dimensionless_unscaled

## This function returns the luminosities escaping in the bins of the specified escape instrument, calculated from
# the total fluxes loaded with the loadEscapeFluxes() function and the distance of the instrument. It returns a tuple
# with three astropy quantity arrays:
#  - the characteristic wavelengths of the instrument's wavelength grid (as for the loadEscapeFluxes() function);
#  - the isotropic-equivalent luminosity \f$4\pi d^2 F_b\f$ of each bin \f$b\f$, with shape (number of bins, number
#    of wavelengths), i.e. the luminosity that an isotropic source would need to emit to produce the flux of the
#    bin; dividing by \f$4\pi\f$ yields the luminosity per unit solid angle escaping in the directions of the bin;
#  - the total luminosity escaping from the model \f$\sum_b 4\pi d^2 f_b F_b\f$, with \f$f_b\f$ the sky fraction of
#    bin \f$b\f$, which equals the luminosity of the single bin of an all-directions binning.
#
# The luminosities have the flavor of the fluxes, with the area unit removed (e.g. W/micron for fluxes in
# W/m2/micron). The function supports models in the local universe only; for a redshifted model, it raises an error.
def loadEscapeLuminosities(instrument):
    if _redshift(instrument.simulation()) > 0:
        raise ValueError("Escape luminosities are not supported for a redshifted model")
    distance = instrument.distance()
    if distance <= 0:
        raise ValueError("Escape instrument '{}' has no distance".format(instrument.name()))
    bins = EscapeBins(instrument)
    waves, fluxes = loadEscapeFluxes(instrument)
    luminosities = (4 * np.pi * distance**2 * fluxes).to(fluxes.unit * u.m**2)
    total = (bins.skyFractions()[:, np.newaxis] * luminosities).sum(axis=0)
    return waves, luminosities, total

# -----------------------------------------------------------------

## This helper function returns a list with the paths of the files with the specified final segment ("sed.dat" or
# "sedstats.dat") for each bin of the specified escape instrument, in the order of the bin index, or None if the
# instrument has no such files. If only some of the bins have such a file, the function raises an error.
def _binFilePaths(instrument, bins, fileType):
    paths = {}
    for path in instrument.outFilePaths(fileType):
        index = escapeBinIndex(path)
        if index is not None and path.name.endswith("_" + fileType): paths[index] = path
    if len(paths) == 0:
        return None
    if sorted(paths.keys()) != list(range(bins.numBins())):
        raise ValueError("Escape instrument '{}' lacks some of its *_{} files".format(instrument.name(), fileType))
    return [ paths[index] for index in range(bins.numBins()) ]

## This helper function returns the path of the data cube file of the specified escape instrument for the
# specified flux component or statistics (e.g. "total" or "stats1"), or raises an error if there is no such file.
def _cubeFilePath(instrument, name):
    fileName = "{}_{}_{}.fits".format(instrument.prefix(), instrument.name(), name)
    paths = [ path for path in instrument.outFilePaths(".fits") if path.name == fileName ]
    if len(paths) != 1:
        raise ValueError("Escape instrument '{}' has no SED files and no {} data cube".format(instrument.name(), name))
    return paths[0]

## This helper function returns the column description in an SED file for the specified flux component name.
def _columnDescription(component):
    descriptions = {
        "total": "total flux",
        "transparent": "transparent flux",
        "primarydirect": "direct primary flux",
        "primaryscattered": "scattered primary flux",
        "secondarytransparent": "transparent secondary flux",
        "secondarydirect": "direct secondary flux",
        "secondaryscattered": "scattered secondary flux",
    }
    if component in descriptions:
        return descriptions[component]
    match = re.fullmatch(r"primaryscattered(\d+)", component)
    if match:
        return "{}-times scattered primary flux".format(match.group(1))
    raise ValueError("Unknown flux component: '{}'".format(component))

## This helper function returns the redshift of the specified simulation, or zero for a model in the local universe.
def _redshift(simulation):
    try:
        cosmology = simulation.getStringAttribute("//MonteCarloSimulation/cosmology/*[1]", None)
    except ValueError:
        return 0.
    if cosmology == "LocalUniverseCosmology":
        return 0.
    return simulation.getFloatAttribute("//MonteCarloSimulation/cosmology/*[1]", "redshift")

# -----------------------------------------------------------------
