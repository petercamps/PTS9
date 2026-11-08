#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.do.convert_healpix Create a projected image based on a raw HEALPix output cube
#
# This script creates a projected image based on the raw output of a HEALPixSkyInstrument, or on a data cube written
# by an EscapeInstrument with a HEALPixEscapeBinning and the DataCubes output style, which uses the same layout of the
# HEALPix pixels. The input can be a FITS file or a FITS file bundle in an HDF5 file (see the pts.utils.isBundlePath()
# function). The output is a FITS file with the projected images, in the units of the input, and the wavelength axis
# of the input; it has no spatial axes, because the projection is not linear.
#
# For the output of an escape instrument, the center of the image corresponds to inclination 90° and azimuth 180°;
# the inclination 0° (the positive z-axis) is at the top of the image, and the azimuth increases to the left. The
# angles \em thetaCenter and \em phiCenter then rotate the image relative to this center.
#
# The function takes the following arguments:
#  - \em inputFileName (string): the name of the file containing the raw HEALPix data cube
#  - \em outputFileName (string): the name of the generated projected image
#  - \em projection (string): the projection to use. Accepted values are Mollweide (for the Mollweide projection)
#    and HammerAitoff (for the Hammer-Aitoff projection); default is Mollweide.
#  - \em nPixelY (int): Number of pixels in the vertical direction for the generated image. The number of pixels
#    in the horizontal direction is twice this value. Default is 250.
#  - \em thetaCenter (float): Zenith angle of the central position relative to the original crosshair of the
#    HEALPixSkyInstrument (in degrees). Corresponds to a vertical rotation that goes downward in the right half
#    of the image. Default is 0.
#  - \em phiCenter (float): Azimuth angle of the central position relative to the original crosshair of the
#    HEALPixSkyInstrument (in degrees). Corresponds to a horizontal rotation to the right. Default is 0.
#

# -----------------------------------------------------------------


def do(
    inputFileName: (
        str,
        "HEALPixSkyInstrument or HEALPix EscapeInstrument output cube containing the raw HEALPix pixel data",
    ),
    outputFileName: (
        str,
        "output image containing an all sky projection of the input HEALPix grid",
    ),
    projection: (str, "projection to use (Mollweide/HammerAitoff)") = "Mollweide",
    nPixelY: (
        int,
        "number of pixels to use in the vertical direction for the output image",
    ) = 250,
    thetaCenter: (
        float,
        "zenith angle of the central position relative to the original crosshair (in degrees),"
        " corresponding to a vertical rotation that goes downward in the right half of the image",
    ) = 0.0,
    phiCenter: (
        float,
        "azimuth angle of the central position relative to the original crosshair (in degrees),"
        " corresponding to a horizontal rotation to the right",
    ) = 0.0,
) -> "create a projected image based on a raw HEALPix output cube":

    import pts.simulation as sm
    import pts.simulation.healpix as healpix
    import pts.utils as ut
    import pathlib
    from pts.utils.error import UserError
    import logging

    if not ut.isBundlePath(inputFileName):
        inputPath = pathlib.Path(inputFileName)
        if not inputPath.exists():
            raise UserError("Input file does not exist!")
        if not inputPath.suffix == ".fits":
            raise UserError("Input file is not a FITS file!")

    outputPath = pathlib.Path(outputFileName)
    if outputPath.exists():
        raise UserError("Output file already exists!")
    if not outputPath.suffix == ".fits":
        raise UserError("Output file is not a FITS file!")

    if nPixelY <= 0:
        raise UserError("Invalid number of pixels: {0}!".format(nPixelY))

    # load the data cube, with indices in the order [x, y, z], and its axes
    data = sm.loadFits(inputFileName)
    axes = sm.getFitsAxes(inputFileName)
    logging.info("Opened input file {0}".format(inputFileName))

    # project the data cube, which the HEALPixGrid class expects with indices in the order [z, y, x]
    hData = healpix.HEALPixGrid(data.value.T)
    hData.printInfo()
    outputData = hData.getProjectionMap(nPixelY, projection, thetaCenter, phiCenter)

    # write the projected images with the wavelength axis of the input, but without the spatial axes of the raw
    # HEALPix layout, which do not apply to the projection
    sm.writeFits(outputPath, outputData.T << data.unit, zaxis=axes[2] if len(axes) > 2 else None)
    logging.info("Wrote output file {0}".format(outputPath))
