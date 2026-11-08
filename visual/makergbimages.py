#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.makergbimages Creating RGB images for SKIRT simulation surface brightness output
#
# The function in this module creates RGB images for surface brightness maps produced by the instruments
# of a SKIRT simulation. The caller can specify the wavelengths corresponding to the R,G,B frames in the image.

# -----------------------------------------------------------------

import logging
import numpy as np
import pts.simulation as sm
import pts.utils as ut
from .rgbimage import RGBImage

# -----------------------------------------------------------------

## This function creates RGB images for surface brightness maps produced by the instruments of a SKIRT simulation.
# Specifically, the function accepts a single Simulation or Instrument instance, and creates one or more RGB images
# for each ".fits" file actually produced by the instruments of the specified simulation or the specified instrument.
#
# If the \em fileType argument is missing, the function by default handles the "*_total.fits" output files. Otherwise,
# the \em fileType argument must be a string indicating the final segment of the name of the instrument output files
# to be handled, excluding the ".fits" suffix. For example, "primaryscattered".
#
# If the \em wavelengthTuples argument is missing, a single image is created for each ".fits" file using
# the frames loaded by default by the RGBImage constructor. Otherwise, the \em wavelengthTuples argument must
# contain a sequence or a dictionary of 3-tuples with (R,G,B) wavelengths. Each of these tuples causes an image
# to be created, loading the frames most closely corresponding to the specified wavelengths.
#
# The \em fromPercentile and \em toPercentile arguments take the percentile values, in range [0,100], used
# to clip the luminosity values loaded from the fits file.
#
# The images are saved in PNG format and have the same name as the corresponding ".fits" file but with a ".pdf"
# filename extension. If there are multiple images per ".fits" file, a serial number (for a sequence of wavelength
# tuples) or the tuple name (for a dictionary of wavelength tuples) is added.
# By default the images are placed next to the corresponding ".fits" files in the same directory, but another output
# directory can be specified if needed.
#
# \note The function uses the same pixel range scaling for all RGB images created for a particular wavelength tuple,
# across multiple instruments. To achieve this, the function makes a separate preprocessing pass over all files.
#
# The function skips escape instruments, whose data cubes hold escape bins rather than images, and logs a message
# for each of them.
#
def makeRGBImages(simulation, wavelengthTuples=None, *, fileType="total", fromPercentile=30, toPercentile=100,
                  outDirPath=None):

    # get the (instrument, output file path) tuples to be handled
    instr_paths = _withoutEscapeInstruments(sm.instrumentOutFilePaths(simulation, fileType+".fits"))
    if len(instr_paths) < 1:
        return

    # convert the wavelength tuples argument into a list of names and wavelength tuples
    if wavelengthTuples is None:
        name_waves = [ ("", None) ]
    elif isinstance(wavelengthTuples, (tuple,list)):
        name_waves = [ ("_" + str(index+1), wavetuple) for index, wavetuple in enumerate(wavelengthTuples) ]
    elif isinstance(wavelengthTuples, dict):
        name_waves = [ ("" if not name else "_" + name, wavetuple) for name, wavetuple in wavelengthTuples.items() ]
    else:
        raise ValueError("Invalid wavelengthTuples type")

    # loop over all wavelength tuples
    for name, wavetuple in name_waves:

        # determine the appropriate pixel range for ALL images
        ranges = []
        for instrument, filepath in instr_paths:
            image = RGBImage(filepath, frameIndices=instrument.wavelengthIndices(wavetuple))
            ranges += list(image.percentilePixelRange(fromPercentile, toPercentile))
        rmin = min(ranges)
        rmax = max(ranges)

        # create an RGB file for each output file
        for instrument, filepath in instr_paths:
            image = RGBImage(filepath, frameIndices=instrument.wavelengthIndices(wavetuple))
            image.setRange(rmin, rmax)
            image.applyLog()
            image.applyCurve()

            # determine output file path
            saveFilePath = ut.savePath(filepath.with_name(filepath.stem+name+".png"), ".png", outDirPath=outDirPath)
            image.saveTo(saveFilePath)
            logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------

## This function creates broadband-convolved RGB images for surface brightness data cubes produced by the instruments
# of a SKIRT simulation. The color makeup of the image is defined by a list of broadbands, where each band can
# contribute with some given weight to each of the R,G,B channels. This allows arbitrary mapping of wavelength ranges
# to RGB color.
#
# After convolving the SKIRT output data cubes with each of the specified broadbands, the function converts the
# surface brightness values for each band to per-frequency flavor (MJy/sr) before performing the RGB color mixing.
#
# Each value in a SKIRT data cube is the mean surface brightness over a wavelength bin. To convolve the cubes in a way
# that conserves the surface brightness in each bin, which matters for emission lines and noisy output (see the
# BroadBand.convolve() function), the function needs the borders of the wavelength bins of each instrument. It loads
# these borders from the output of an InstrumentWavelengthGridProbe in the simulation, and passes them to that function.
# If the simulation has no such probe output for an instrument, or if the wavelength grid of the instrument does not
# consist of disjoint bins (for example, because it consists of bands), the function logs a warning and interpolates
# the surface brightness between the characteristic wavelengths of the bins instead.
#
# The final images are saved in PNG format and have the same name as the corresponding ".fits" file, possibly
# extended with a custom name segment. By default the images are placed next to the corresponding ".fits" files in
# the same directory, but another output directory can be specified if needed.
#
# The function takes the following arguments:
#  - \em simulation: a single Simulation or Instrument instance for which to create RGB images. When given a
#        simulation, the function uses the same pixel range scaling for the RGB images across all instruments.
#  - \em contributions: a sequence of 4-tuples (band, r, g, b), each containing a BroadBand instance and three
#        weights that specify the contribution of the surface brightness in this band to each RGB channel.
#  - \em name: a string that will be added to the image file name; defaults to the empty string.
#  - \em fileType: a string indicating the final segment of the name of the instrument output files
#        to be handled, excluding the ".fits" suffix (for example: "primaryscattered"); if this argument is
#        missing, the function by default handles the "*_total.fits" output files.
#  - \em decades: the dynamic range for the surface brightness values in the images, expressed in order of magnitudes
#        (decades); the default value is 3. Pixels with a value less than \em decades order of magnitudes below
#        the maximum value (in any of the channels/images) are clipped to zero.
#  - \em fmax: the largest surface brightness value that will be shown in the images without clipping;
#        this must be an astropy quantity with per-frequency surface brightness units (e.g. MJy/sr).
#        If this argument is None or missing, the function uses the largest value found in any of the channels/images.
#  - \em fmin: the smallest surface brightness value that will be shown in the image without clipping;
#        this must be an astropy quantity with per-frequency surface brightness units (e.g. MJy/sr).
#        If this argument is None or missing, the function determines this value from \em fmax and \em decades.
#        using the formula \f$f_\mathrm{range}=\log_{10}(f_\mathrm{max}/f_\mathrm{min})\f$.
#        If \em fmin is specified, the value of \em decades is ignored.
#  - \em imageDirPath: path to the PNG output directory; default is next to the fits files.
#
#  The function returns a tuple (fmin, fmax) specifying the surface brightness range actually used in the images,
#  as an astropy quantity with per-frequency surface brightness units (MJy/sr).
#
#  The function skips escape instruments, whose data cubes hold the flux density of escape bins rather than surface
#  brightness images, and logs a message for each of them.
#
def makeConvolvedRGBImages(simulation, contributions, name="", *, fileType="total", decades=3,
                           fmax=None, fmin=None, outDirPath=None):

    # get the (instrument, output file path) tuples to be handled
    instr_paths = _withoutEscapeInstruments(sm.instrumentOutFilePaths(simulation, fileType+".fits"))
    if len(instr_paths) > 0:
        name = "" if not name else "_" + name
        sbunit = sm.unit("MJy/sr")

        # construct an image object with integrated color channels for each output file
        # and keep track of the largest surface brightness value (in MJy/sr) in each of the images
        images = []
        fmaxes = []
        bordersPerInstrument = {}
        for instrument, filepath in instr_paths:
            logging.info("Convolving for RGB image {}{}".format(filepath.stem, name))

            # get the data cube in its intrinsic units
            cube = sm.loadFits(filepath)

            # get the wavelengths and, once for each instrument, the wavelength bin borders (or None)
            wavelengths = instrument.wavelengths()
            if instrument.name() not in bordersPerInstrument:
                bordersPerInstrument[instrument.name()] = _loadBinBorders(instrument)
            left, right = bordersPerInstrument[instrument.name()]

            # initialize an RGB frame
            dataRGB = np.zeros( (cube.shape[0], cube.shape[1], 3) ) << sbunit

            # add color for each filter
            for band,w0,w1,w2 in contributions:
                # convolve and convert to per frequency units
                data = band.convolve(wavelengths, cube, flavor=sbunit, leftBorders=left, rightBorders=right)
                dataRGB[:,:,0] += w0*data
                dataRGB[:,:,1] += w1*data
                dataRGB[:,:,2] += w2*data

            # construct the image object
            images.append(RGBImage(dataRGB.value))
            fmaxes.append(dataRGB.max())

        # determine the appropriate pixel range for all output images
        fmax = max(fmaxes) if fmax is None else fmax << sbunit
        fmin = fmax/10**decades if fmin is None else fmin << sbunit

        # create an RGB file for each output file
        for (instrument, filepath), image in zip(instr_paths,images):
            image.setRange(fmin.value, fmax.value)
            image.applyLog()
            image.applyCurve()

            # determine output file path
            saveFilePath = ut.savePath(filepath.with_name(filepath.stem+name+".png"), ".png", outDirPath=outDirPath)
            image.saveTo(saveFilePath)
            logging.info("Created convolved RGB image file {}".format(saveFilePath))

        # return the surface brightness range used for these images
        return fmin,fmax

## This private helper function returns the specified list of (instrument, output file path) tuples without the
# tuples for escape instruments, whose data cubes hold escape bins rather than images, and logs a message for each
# escape instrument that it leaves out.
def _withoutEscapeInstruments(instr_paths):
    escapeNames = { instrument.name() for instrument, _ in instr_paths if instrument.type() == "EscapeInstrument" }
    for name in sorted(escapeNames):
        logging.info("Skipping escape instrument {}: its data cubes hold escape bins rather than images".format(name))
    return [ (instrument, filepath) for instrument, filepath in instr_paths if instrument.type() != "EscapeInstrument" ]

## This private helper function returns a tuple with the left and right borders of the wavelength bins of the
# specified instrument, loaded from the output of an InstrumentWavelengthGridProbe, if the simulation has such output
# and if the bins are disjoint, i.e. if the effective width of each bin equals the distance between its borders.
# Otherwise, the function logs a warning and returns a tuple with two None values, so that the band convolution
# interpolates between the characteristic wavelengths of the bins instead.
def _loadBinBorders(instrument):
    try:
        left, right = instrument.wavelengthBinBorders()
        widths = instrument.wavelengthBinWidths()
    except ValueError:
        logging.warning("No wavelength bin borders for instrument {}; add an InstrumentWavelengthGridProbe to the "
                        "simulation to convolve over the wavelength bins; interpolating between the characteristic "
                        "wavelengths instead".format(instrument.name()))
        return None, None
    # allow for the rounding of the values written to the probe output file
    if not np.allclose(widths.to_value(left.unit), np.abs(right - left).to_value(left.unit), rtol=1e-3, atol=0):
        logging.warning("The wavelength grid of instrument {} does not consist of disjoint bins; interpolating between "
                        "the characteristic wavelengths instead".format(instrument.name()))
        return None, None
    return left, right

# -----------------------------------------------------------------
