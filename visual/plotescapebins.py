#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotescapebins Plot the flux escaping from the model as a function of direction
#
# This module offers a function to plot the output of SKIRT escape instruments as a function of the escape direction.
#

# -----------------------------------------------------------------

import logging
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pts.simulation as sm
import pts.simulation.healpix as healpix
import pts.utils as ut
from .makergbimages import _loadBinBorders

# -----------------------------------------------------------------

## This function plots the flux recorded by one or more escape instruments as a function of the escape direction,
# at a single wavelength or convolved with a broadband. It accepts a sequence of Simulation and/or Instrument
# instances (or a single instance of either of these types), and creates a plot for each EscapeInstrument among them,
# except for instruments with an all-directions binning, which have a single bin. The plot depends on the binning:
#  - inclination bins: a step plot of the flux versus inclination, with error bars if the instrument recorded
#    statistics (only when plotting at a single wavelength); with folded hemispheres, the inclination runs from 0° to
#    90°, and each bin also includes the mirror directions;
#  - inclination and azimuth bins: a map of the flux versus azimuth (horizontal) and inclination (vertical, with 0°
#    at the top), with a logarithmic color scale;
#  - HEALPix bins: a Mollweide or Hammer-Aitoff projection of the sky of escape directions, with a logarithmic color
#    scale; the inclination 0° (the positive z-axis) is at the top, the azimuth is 180° at the center of the map and
#    increases to the left, as for the pts.simulation.do.convert_healpix command.
#
# The function loads the fluxes with the pts.simulation.loadEscapeFluxes() function, so that it handles both output
# styles of an escape instrument, and output files written as bundles in an HDF5 file. It accepts the following
# arguments:
#  - \em simulation: the simulations and/or instruments to be handled, as described above.
#  - \em wavelength: an astropy quantity specifying the wavelength (or frequency or energy) at which to plot the flux;
#    the function uses the characteristic wavelength of the instrument's wavelength grid nearest to it.
#  - \em band: a BroadBand instance (see the pts.band package) with which to convolve the flux; the function then
#    plots the mean flux density over the band. If the simulation includes the output of an
#    InstrumentWavelengthGridProbe for the instrument, the convolution uses the borders of the wavelength bins (see
#    the BroadBand.convolve() function). If both \em wavelength and \em band are None, the function plots the flux at
#    the wavelength where the total escaping flux, summed over the bins and weighted by their sky fractions, is
#    largest. The \em wavelength and \em band arguments cannot both be specified.
#  - \em component: the flux component to be plotted, as described for the pts.simulation.loadEscapeFluxes() function;
#    the default is "total".
#  - \em decades: the number of decades in the range of the vertical axis or color scale; the default is 3.
#  - \em projection: the projection for HEALPix bins, "Mollweide" (the default) or "HammerAitoff".
#  - \em outDirPath: the directory in which the plots are saved; by default, the directory of the instrument's output.
#  - \em figSize: the size of the figure in inch; the default is (8, 6).
#  - \em interactive: see below.
#
# By default, each plot is saved as a PDF file named after the simulation prefix and the instrument name, ending in
# "_escape.pdf". In interactive mode (see the pts.utils.interactive() function), the figures are not saved and they
# are left open so that they are displayed in notebooks.
#
def plotEscapeBins(simulation, *, wavelength=None, band=None, component="total", decades=3, projection="Mollweide",
                   outDirPath=None, figSize=(8, 6), interactive=None):
    if wavelength is not None and band is not None:
        raise ValueError("The wavelength and band arguments cannot both be specified")

    # get the escape instruments
    entities = [ simulation ] if isinstance(simulation, (sm.Simulation, sm.Instrument)) else simulation
    instruments = []
    for entity in entities:
        instruments += entity.instruments() if isinstance(entity, sm.Simulation) else [ entity ]
    instruments = [ instrument for instrument in instruments if instrument.type() == "EscapeInstrument" ]

    for instrument in instruments:
        bins = sm.loadEscapeBins(instrument)
        if bins.numBins() < 2:
            logging.info("Skipping escape instrument {}: it has a single bin".format(instrument.name()))
            continue

        # get the flux in each bin and a text describing the wavelength or band
        waves, fluxes = sm.loadEscapeFluxes(instrument, component)
        errors = None
        if band is not None:
            left, right = _loadBinBorders(instrument)
            values = band.convolve(waves, fluxes, leftBorders=left, rightBorders=right)
            spectralText = "in band {}".format(band.name())
            if not isinstance(values, u.Quantity):
                logging.warning("Skipping escape instrument {}: band {} does not overlap its wavelength grid".format(
                                instrument.name(), band.name()))
                continue
        else:
            if wavelength is not None:
                index = np.argmin(np.abs(waves - wavelength.to(waves.unit, equivalencies=u.spectral())))
            else:
                index = np.argmax((bins.skyFractions()[:, np.newaxis] * fluxes.value).sum(axis=0))
            values = fluxes[:, index]
            spectralText = "at {:.4g} {}".format(waves[index].value, waves.unit.to_string())
            try:
                errors = sm.loadEscapeRelativeErrors(instrument)[1][:, index].value * values
            except ValueError:
                pass
        title = "{} {} {} {}".format(instrument.prefix(), instrument.name(), component, spectralText)

        # determine the range of the values
        vmax = values.value.max()
        if not vmax > 0:
            logging.warning("Skipping escape instrument {}: no positive flux {}".format(instrument.name(),
                                                                                     spectralText))
            continue
        vmin = vmax * 10**(-decades)
        fluxLabel = sm.latexForSpectralFluxWithUnit(values.unit)

        plt.figure(figsize=figSize)
        if bins.binningType() == "InclinationEscapeBinning" and bins.inclinationAzimuthShape()[1] == 1:
            # step plot versus inclination
            tmin, tmax = bins.inclinationRanges()
            edges = np.append(tmin.to_value(u.deg), tmax[-1].to_value(u.deg))
            plt.stairs(values.value, edges, baseline=None, color='b', lw=1.5)
            if errors is not None:
                plt.errorbar(0.5 * (edges[:-1] + edges[1:]), values.value, yerr=errors.value, fmt='none', ecolor='b',
                             elinewidth=1, capsize=3)
            plt.yscale('log')
            plt.ylim(vmin, vmax * 10**0.2)
            plt.xlim(edges[0], edges[-1])
            plt.xlabel("{}inclination (deg)".format("folded " if bins.isFolded() else ""), fontsize='large')
            plt.ylabel(fluxLabel, fontsize='large')

        elif bins.binningType() == "InclinationEscapeBinning":
            # map versus azimuth and inclination
            numInclination, numAzimuth = bins.inclinationAzimuthShape()
            tmin, tmax = bins.inclinationRanges()
            pmin, pmax = bins.azimuthRanges()
            tedges = np.append(tmin.to_value(u.deg)[::numAzimuth], tmax[-1].to_value(u.deg))
            pedges = np.append(pmin.to_value(u.deg)[:numAzimuth], pmax[numAzimuth-1].to_value(u.deg))
            grid = values.value.reshape(numInclination, numAzimuth)
            plt.pcolormesh(pedges, tedges, np.ma.masked_less_equal(grid, 0), cmap='viridis', edgecolors='face',
                           norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
            plt.gca().invert_yaxis()
            plt.xlabel("azimuth (deg)", fontsize='large')
            plt.ylabel("{}inclination (deg)".format("folded " if bins.isFolded() else ""), fontsize='large')
            plt.colorbar(label=fluxLabel)

        else:
            # projection of the HEALPix sky; place the values in the raw layout used by the HEALPixSkyInstrument
            numSide = 2**bins.healpixOrder()
            frame = np.zeros((4 * numSide - 1, 4 * numSide))
            x, y = bins.framePixels()
            frame[y, x] = values.value
            image = healpix.HEALPixGrid(frame).getProjectionMap(250, projection)
            plt.imshow(np.ma.masked_less_equal(image, 0), origin='lower', cmap='viridis', interpolation='nearest',
                       norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
            plt.xticks([])
            plt.yticks([])
            plt.xlabel("{} projection; inclination 0° at the top; azimuth 180° at the center, increasing to the left"
                       .format(projection), fontsize='small')
            plt.colorbar(label=fluxLabel, orientation='horizontal', shrink=0.8, pad=0.08)

        plt.title(title)

        # if not in interactive mode, save the figure; otherwise leave it open
        if not ut.interactive(interactive):
            defSaveFilePath = instrument.outFilePaths("bins.dat")[0]
            defSaveFilePath = defSaveFilePath.with_name("{}_{}_escape.pdf".format(instrument.prefix(),
                                                                                 instrument.name()))
            saveFilePath = ut.savePath(defSaveFilePath, (".pdf",".png"), outDirPath=outDirPath)
            plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
            plt.close()
            logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------
