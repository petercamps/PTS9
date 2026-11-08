#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.makewavelengthmovie Creating a movie that runs through all wavelengths in the SKIRT simulation output.
#
# The function in this module creates a movie for the output of a SKIRT simulation. The movie combines
# the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
# running through all wavelengths in the simulation.

# -----------------------------------------------------------------

import logging
import numpy as np
import pts.simulation as sm
import pts.utils as ut
from .moviefile import MovieFile
from .rgbimage import RGBImage

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

# -----------------------------------------------------------------

## This function creates a movie for the output of the specified simulation. The movie combines the SEDs (bottom panel)
# and the pixel frames (top panel, from left to right) for up to three instruments, running through all wavelengths
# in the simulation.
#
# The function takes the following arguments:
#  - simulation:    a sequence of Simulation and/or Instrument instances (or a single instance of one of these types).
#                   The sequence is expanded into a list of up to three instruments; any extra instruments are ignored.
#                   Escape instruments are ignored. All instruments in the final list must include an integrated SED
#                   ("*_sed.dat" file) as well as total flux pixel frames ("*_total.fits" file), which may be written
#                   by two instruments with the same name, such as an SEDInstrument and a FrameInstrument; otherwise
#                   the function does nothing. All instruments must use the same wavelength grid, and the instruments
#                   frames must have the same pixel size. These requirements are not verified.
#                   The fluxes and surface brightness values are plotted in the units of the SKIRT output.
#  - maxPercentile: the percentile value in range [0,100] used to determine the maximum surface brightness in the
#                   images. The default value is 100, so that the largest surface brightness value in the frame(s)
#                   determines the maximum value. A smaller percentile value such as 99 could be used to exclude the
#                   very largest surface brightness values (for example, the direct light from a point source).
#  - minPercentile: the percentile value in range [0,100] used to determine the minimum surface brightness in the
#                   images. The default value is 10. This argument is ignored if \em decades is specified.
#  - decades:       the number of decades in the surface brightness range in the images. This is an alternative way
#                   to determine the minimum value once the maximum has been obtained using  \em maxPercentile.
#                   By default, the \em minPercentile mechanism is used instead.
#  - renormalize:   a Boolean flag selecting one of two image scaling options. When False (the default value), the
#                   surface brightness range is determine for the complete data cube and kept constant for all frames.
#                   As a result, image frames might be completely dark at wavelengths with very low overall luminosity.
#                   When True, the surface brightness range is determined for each frame separately. This allows to see
#                   the spatial structure of the surface brightness even at wavelengths with very low luminosity.
#  - outDirPath:    when specified, overrides the default output directory path.
#  - outFileName:   when specified, overrides the filename portion of the default output path.
#  - outFilePath:   when specified, overrides the the complete output file path.
#  - rate:          the frame rate of the movie, in frames per second. The default value is 7.
#
# By default, the movie is saved in the output directory of the first instrument, using a name starting with the
# corresponding simulation prefix and ending with ".mp4". This can be overridden with the out* arguments as described
# for the pts.utils.savePath() function.
#
def makeWavelengthMovie(simulation, *, maxPercentile=100, minPercentile=10, decades=None, renormalize=False,
                        outDirPath=None, outFileName=None, outFilePath=None, rate=7):

    # get the list of instruments and corresponding output file paths, ignoring escape instruments
    sedEntries = [ (instrument, path) for instrument, path in sm.instrumentOutFilePaths(simulation, "sed.dat")
                   if instrument.type() != "EscapeInstrument" ][:3]
    cubEntries = [ (instrument, path) for instrument, path in sm.instrumentOutFilePaths(simulation, "total.fits")
                   if instrument.type() != "EscapeInstrument" ][:3]
    if len(sedEntries) < 1 or len(sedEntries) != len(cubEntries) \
                           or any([ a.name() != b.name() for (a,_),(b,_) in zip(sedEntries,cubEntries) ]):
        return
    instrA, sedPaths = zip(*sedEntries)
    instrB, cubPaths = zip(*cubEntries)
    instruments = instrA

    # get the wavelength grid for the first instrument (assumed to be the same for all instruments)
    wavelengths = instruments[0].wavelengths()
    if len(wavelengths) < 3:
        return
    nlambda = len(wavelengths)

    # load the data
    logging.info("Creating movie for {} ({} wavelengths and {} instruments)..." \
                 .format(instruments[0].prefix(), nlambda, len(instruments)))
    sedData = [ sm.loadColumns(sedPath, "total flux")[0] for sedPath in sedPaths ]
    cubData = [ sm.loadFits(cubPath).value for cubPath in cubPaths ]

    # determine the shape (assume that frames in all fits files have the same shape)
    imgShape = cubData[0].shape[:2]
    sedShape = (len(cubData)*imgShape[0], max(imgShape[1]//2, 300))
    totalShape = (sedShape[0], imgShape[1]+sedShape[1])

    # determine the global surface brightness range
    fmax = max([ np.percentile(np.unique(cube), maxPercentile) for cube in cubData ])
    if decades is None:
        fmin = min([ np.percentile(np.unique(cube), minPercentile) for cube in cubData ])
    else:
        fmin = fmax / 10**decades

    # determine the global integrated flux range
    Fmax = max([ sed.max() for sed in sedData ])
    Fmin = Fmax * fmin / fmax

    # open the movie file
    defSaveFilePath = sedPaths[0].with_name(instruments[0].prefix() + "_wavemovie.mp4")
    saveFilePath = ut.savePath(defSaveFilePath, (".mp4",),
                               outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
    movie = MovieFile(saveFilePath, shape=totalShape, rate=rate)

    # for each wavelength, construct and add a movie frame
    for frame in range(nlambda):
        logging.info("  adding frame " + str(frame+1)+"/"+str(nlambda) + "...")

        # determine the surface brightness range for this frame, if needed
        if renormalize:
            fmax = max([np.percentile(np.unique(cube[:,:,frame]), maxPercentile) for cube in cubData])
            if decades is None:
                fmin = min([np.percentile(np.unique(cube[:,:,frame]), minPercentile) for cube in cubData])
            else:
                fmin = fmax / 10 ** decades

        # assemble the top panel
        image = None
        for cube in cubData:
            im = RGBImage(np.dstack(( cube[:,:,frame],cube[:,:,frame],cube[:,:,frame] )))
            im.setRange(fmin,fmax)
            im.applyLog()
            im.applyColorMap("gnuplot2")
            if image==None: image = im
            else: image.addRight(im)

        # plot the seds in the bottom panel
        dpi = 100
        figure = Figure(dpi=dpi, figsize=(sedShape[0]/dpi, sedShape[1]/dpi), facecolor='w', edgecolor='w')
        canvas = FigureCanvasAgg(figure)
        ax = figure.add_subplot(111)
        colors = ('r','g','b')
        for sed, instrument, color in zip(sedData, instruments, colors):
             ax.loglog(wavelengths.value, sed.value, color=color, label=instrument.name())
        ax.axvline(wavelengths[frame].value, color='m')
        ax.set_ylabel(sm.latexForSpectralFlux(sedData[0]) + sm.latexForUnit(sedData[0]))
        ax.set_ylim(Fmin.value/1.1, Fmax.value*1.1)
        ax.legend(loc='lower right', title=sm.latexForWavelength(wavelengths) \
                    + r"$={0:.4g}\,$".format(wavelengths[frame].value)+sm.latexForUnit(wavelengths))
        canvas.draw()
        im = RGBImage(figure)
        image.addBelow(im)

        # add the frame to the movie
        movie.addFrame(image)

    # close the movie file
    movie.close()
    logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------
