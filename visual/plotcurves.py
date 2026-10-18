#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotcurves Plot SEDs or other SKIRT output as a function of wavelength
#
# This module offers functions to plot spectral distributions or other SKIRT output as a function of wavelength.
# Each function accepts a list of Instrument or Probe instances of a certain type, and plots the corresponding output.
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import numpy as np
import pts.simulation as sm
import pts.storedtable as stab
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates a plot of the SEDs produced by one or more SKIRT simulations. Specifically, the function
# accepts a sequence of Simulation and/or Instrument instances (or a single instance of either of these types),
# and it plots the total flux density for each of the instruments that actually produced an "*_sed.dat" file.
#
# If there is only a single SED, and the corresponding output file includes the individual flux components,
# the function plots these as well.
#
# If there are multiple SEDs, the plot uses the flux flavor and units of the first SED and the data for all
# other SEDs are converted, if needed.
#
# If there is no SED in the list, the function does nothing.
#
# By default, the figure is saved in the output directory of the first instrument in the list,
# using a name starting with the corresponding simulation prefix, possibly including the instrument name,
# and ending with ".pdf". This can be overridden with the out* arguments as described for the
# pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSeds(simulation, minWavelength=None, maxWavelength=None, decades=None, *,
             outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 6), interactive=None):

    # private function to get the maximum flux within the wavelength range passed to the plotSeds function
    def maxFluxInRange(flux, wave):
        wmin = minWavelength if minWavelength is not None else wave[0]
        wmax = maxWavelength if maxWavelength is not None else wave[-1]
        mask = (wave>=wmin) & (wave<=wmax)
        if np.count_nonzero(mask) > 0:
            return flux[mask].max()
        else:
            return flux.max()

    # get the (instrument, output file path) tuples
    instr_paths = sm.instrumentOutFilePaths(simulation, "sed.dat")
    if len(instr_paths) < 1:
        return

    # setup the figure
    plt.figure(figsize=figSize)

    # if there is a single output file, and it has components, plot the components
    if len(instr_paths) == 1 and any(["transparent" in col for col in sm.getColumnDescriptions(instr_paths[0][1])]):
        instrument, filepath = instr_paths[0]
        # load the columns (we assume that all components have the same units)
        wave, tot, tra, dirpri, scapri, dirsec, scasec = sm.loadColumns(filepath, (0, 1, 2, 3, 4, 5, 6))
        waveUnit = wave.unit
        fluxUnit = tot.unit
        fluxMax = max(maxFluxInRange(tot, wave), maxFluxInRange(tra, wave))
        # plot the various components
        label = "{} {} ".format(instrument.prefix(), instrument.name())
        plt.plot(wave.value, tot.value, color='k', ls='solid', label=label + "total")
        plt.plot(wave.value, tra.value, color='b', ls='dotted', label=label + "transparent")
        plt.plot(wave.value, (dirpri + scapri).value, color='b', ls='solid', label=label + "primary")
        plt.plot(wave.value, (dirsec + scasec).value, color='r', ls='solid', label=label + "secondary")

    # otherwise loop over all SEDs
    else:
        colors = ('r', 'g', 'b', 'c', 'm', 'y')
        colorindex = 0
        first = True
        for instrument, filepath in instr_paths:
            # load the total flux; first time remember units; thereafter convert units
            wave, flux = sm.loadColumns(filepath, (0, 1))
            if first:
                waveUnit = wave.unit
                fluxUnit = flux.unit
                fluxMax = maxFluxInRange(flux, wave)
                first = False
            else:
                wave <<= waveUnit
                flux = sm.convertToFlavor(wave, flux, fluxUnit)
                fluxMax = max(fluxMax, maxFluxInRange(flux, wave))
            # plot
            plt.plot(wave.value, flux.value, color=colors[colorindex],
                     label="{} {} total".format(instrument.prefix(), instrument.name()))
            # advance color index
            colorindex = (colorindex + 1) % len(colors)

    # set axis details and add a legend
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(_adjustWavelengthRange(plt.xlim(), waveUnit, minWavelength, maxWavelength))
    if decades is not None:
        plt.ylim(fluxMax.value * 10 ** (-decades), fluxMax.value * 10 ** 0.2)
    plt.xlabel(sm.latexForWavelengthWithUnit(waveUnit), fontsize='large')
    plt.ylabel(sm.latexForSpectralFluxWithUnit(fluxUnit), fontsize='large')
    plt.legend(loc='best')

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        # use the first instrument output path; if there are multiple instruments, remove the instrument name
        defSaveFilePath = instr_paths[0][1]
        if len(instr_paths) > 1:
            defSaveFilePath = defSaveFilePath.with_name(instr_paths[0][0].prefix() + "_sed.pdf")
        saveFilePath = ut.savePath(defSaveFilePath, (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

## This function creates a plot of the behavior of the sources in a SKIRT simulations. Specifically, it plots
# the luminosity of and the number of photon packets launched by each source as a function of wavelength.
# Both curves are normalized so that their maximum value is equal to one. This allows a user to evaluate
# whether the photon packet distribution (as a function of wavelength) is appropriate for the source spectra.
#
# The function accepts a single Simulation instance and it assumes that the simulation includes exactly one
# LuminosityProbe and exactly one LaunchedPacketsProbe. If this is not the case, the function does nothing.
#
# By default, the figure is saved in the simulation output directory, using a name starting with the simulation prefix
# and ending with "sources.pdf". This can be overridden with the out* arguments as described for the
# pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSources(simulation, minWavelength=None, maxWavelength=None, decades=None, *,
                outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 6), interactive=None):

    # find the required probes
    lumiProbes = simulation.probes("LuminosityProbe")
    packProbes = simulation.probes("LaunchedPacketsProbe")
    if len(lumiProbes) != 1 or len(packProbes) != 1:
        return
    lumiProbe = lumiProbes[0]
    packProbe = packProbes[0]

    # load the luminosities
    lumiFilePath = lumiProbe.outFilePaths("luminosities.dat")[0]
    descriptions = sm.getColumnDescriptions(lumiFilePath)
    columns = list(range(len(descriptions)))
    del columns[1]  # remove the "specific luminosity column"
    lumiWave, lumiTot, *lumiFracs = sm.loadColumns(lumiFilePath, columns)

    # load the number of launched packets
    packFilePath = packProbe.outFilePaths("launchedpackets.dat")[0]
    descriptions = sm.getColumnDescriptions(packFilePath)
    columns = list(range(len(descriptions)))
    packWave, packTot, *packSplits = sm.loadColumns(packFilePath, columns)
    packWave <<= lumiWave.unit

    # setup the figure
    plt.figure(figsize=figSize)
    label = "{} ".format(simulation.prefix())

    # plot the total
    lumiMax = lumiTot.max()
    packMax = packTot.max()
    plt.plot(lumiWave.value, lumiTot/lumiMax, color='k', ls='solid', label=label + "total")
    plt.plot(packWave.value, packTot/packMax, color='k', ls='dashed')

    # loop over all sources
    colors = ('r', 'g', 'b', 'c', 'm', 'y')
    colorindex = 0
    sourceindex = 1
    for lumiFrac, packSplit in zip(lumiFracs, packSplits):
        plt.plot(lumiWave.value, lumiFrac*lumiTot/lumiMax, color=colors[colorindex], ls='solid',
                 label=label + str(sourceindex))
        plt.plot(packWave.value, packSplit/packMax, color=colors[colorindex], ls='dashed')
        # advance color and source index
        colorindex = (colorindex + 1) % len(colors)
        sourceindex += 1

    # set axis details and add a legend
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(_adjustWavelengthRange(plt.xlim(), lumiWave.unit, minWavelength, maxWavelength))
    if decades is not None:
        plt.ylim(10 ** (-decades), 10 ** 0.2)
    plt.xlabel(sm.latexForWavelengthWithUnit(lumiWave.unit), fontsize='large')
    plt.ylabel(r"Normalized $L$ and $N_\mathrm{pp}$", fontsize='large')
    plt.legend(loc='best')

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        saveFilePath = ut.savePath(simulation.outFilePath("sources.pdf"), (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

## This function creates a plot of the spectral resolution \f$R=\lambda/\Delta\lambda\f$ of a wavelength grid
# loaded from a SKIRT stored table (".stab"), a SKIRT text column file (".dat") that includes a wavelength axis,
# or a SKIRT instrument or probe data cube that includes a wavelength axis. A text column file or a data cube can
# also be specified as a bundle in an HDF5 file (see the pts.utils.isBundlePath() function).
#
# By default, the figure is saved in the current directory, using the name of the input file but
# with the ".pdf" filename extension. This can be overridden with the out* arguments as described for the
# pts.utils.savePath() function. In interactive mode (see the pts.utils.interactive() function),
# the figure is not saved and it is left open so that is displayed in notebooks.
#
def plotSpectralResolution(inFilePath, minWavelength=None, maxWavelength=None, decades=None, *, title=None,
                outDirPath=None, outFileName=None, outFilePath=None, figSize=(8, 5), interactive=None):

    # load the wavelength grid
    inFilePath = ut.BundlePath(inFilePath) if ut.isBundlePath(inFilePath) else ut.absPath(inFilePath)
    if inFilePath.suffix.lower() == ".stab":
        table = stab.readStoredTable(inFilePath)
        if "lambda" not in table:
            raise ValueError("No wavelength axis in stored table: {}".format(inFilePath))
        grid = table["lambda"]
    elif inFilePath.suffix.lower() == ".dat":
        if "wavelength" not in sm.getColumnDescriptions(inFilePath)[0].lower():
            raise ValueError("First text column is not labeled 'wavelength': {}".format(inFilePath))
        grid = sm.loadColumns(inFilePath, "1")[0]
    elif inFilePath.suffix.lower() == ".fits":
        axes = sm.getFitsAxes(inFilePath)
        if len(axes) != 3:
            raise ValueError("FITS file does not have embedded wavelength axis")
        grid = axes[2]
    else:
        raise ValueError("Filename does not have the .stab, .dat, or .fits extension: {}".format(inFilePath))

    # calculate the spectral resolution
    R = grid[:-1] / (grid[1:] - grid[:-1])
    Rmax = R.max()

    # choose wavelength units from grid
    wunit = grid.unit

    # setup the plot
    plt.figure(figsize=figSize)
    plt.xlabel(sm.latexForWavelengthWithUnit(wunit), fontsize='large')
    plt.ylabel(r"$R=\frac{\lambda}{\Delta\lambda}$", fontsize='large')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='major', axis='both', ls=":")
    plt.xlim(_adjustWavelengthRange(plt.xlim(), wunit, minWavelength, maxWavelength))
    if decades is not None:
        plt.ylim(Rmax* 10 ** (-decades), Rmax * 10 ** 0.2)

    # plot the spectral resolution
    if title is None or len(title)==0: title = inFilePath.stem
    label = "{}\n{} pts from {:g} to {:g} {}".format(title, len(grid), grid[0].to_value(wunit), grid[-1].to_value(wunit),
                                              sm.latexForUnit(wunit))
    plt.plot(grid[:-1].to_value(wunit), R, label=label)
    plt.legend()

    # if not in interactive mode, save the figure; otherwise leave it open
    if not ut.interactive(interactive):
        saveFilePath = ut.savePath(inFilePath.stem+".pdf", (".pdf",".png"),
                                   outDirPath=outDirPath, outFileName=outFileName, outFilePath=outFilePath)
        plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        logging.info("Created {}".format(saveFilePath))

# ----------------------------------------------------------------------

## This private function adjusts the given wavelength range to the given minimum and/or maximum, if present,
# and returns the result. The input and output wavelength ranges are expressed in the specified units but
# do not carry explicit unit information. The output range will always be in increasing order. The minimum
# and/or maximum values, if present, must include explicit unit information and will be converted to plot units.
# Their ordering depends on the style of the units in which they are specified (and not on the plot units).
def _adjustWavelengthRange(wrange, wunit, wmin=None, wmax=None):
    # get the incoming range in increasing order
    wrange = list(wrange)
    if (wrange[0] > wrange[1]):
        wrange = wrange[::-1]

    # adjust for given minimum and/or maximum
    if wmin is not None:
        wrange[1 if sm.hasReverseWavelengthOrder(wmin) else 0] = wmin.to_value(wunit)
    if wmax is not None:
        wrange[0 if sm.hasReverseWavelengthOrder(wmax) else 1] = wmax.to_value(wunit)

    # get the outgoing range in increasing order
    if (wrange[0] > wrange[1]):
        wrange = wrange[::-1]
    return wrange

# ----------------------------------------------------------------------
