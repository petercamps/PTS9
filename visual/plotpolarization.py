#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.plotpolarization Plot polarization maps from SKIRT simulation output
#
# The function in this module creates PDF polarization maps for the output of polarization-enabled SKIRT simulations
# that include instruments generating surface brightness frames or data cubes (*_total.fits and *_stokes*.fits files).
#

# -----------------------------------------------------------------

import logging
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.lines
import matplotlib.patches
import numpy as np
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function creates PDF polarization maps for the output of polarization-enabled SKIRT simulations
# that include instruments generating surface brightness frames or data cubes (*.fits files).
# Specifically, the function accepts a sequence of Simulation and/or Instrument instances (or a single instance of
# either of these types), and it creates polarization maps for each of the instruments that actually produced
# both *_total.fits and \*_stokes\*.fits files. Other instruments are silently ignored.
#
# By default, the figures are saved in the output directory of the corresponding instrument, using a name starting
# with the simulation prefix and instrument name, and ending with ".pdf". The output directory can be overridden.
# In interactive mode (see the pts.utils.interactive() function), the figures are not saved and are left open
# for display in notebooks.
#
# The function takes the following arguments:
#   - simulation: a sequence of Simulation and/or Instrument instances, or a single instance of either of these types.
#   - plotLinMap: whether to plot a linear polarization map (binned degree and angle); True or False.
#   - plotDegMap: whether to plot a linear polarization degree map; True or False.
#   - plotDegAvg: whether to plot the y-axis averaged linear polarization degree for all x pixels; True or False.
#   - plotCirMap: whether to plot a circular polarization map; True or False.
#   - wavelength: an astropy quantity or a sequence of astropy quantities specifying the wavelength(s) for which to
#     create the plot, or the string 'all' to create a plot for all wavelengths in the data cube, or None to use the
#     first frame in the data cube.
#   - binSize: the number of pixels in each bin, in horizontal and vertical directions.
#   - degreeScale: the maximum polarization degree to be visualized, or None for automatic scale.
#   - decades: number of surface brightness decades (dex) shown in the background image; default is 5.
#   - outDirPath: string or Pathlib.Path object that specifies (overrides) the output directory.
#   - figSize: the horizontal and vertical size of the output figure in inch; default is 6x6 inch.
#   - interactive: whether to leave figures open (True) or save them to file (False).
#
def plotPolarization(simulation, *, plotLinMap=True, plotDegMap=False, plotDegAvg=False, plotCirMap=False,
                     wavelength=None, binSize=(7,7), degreeScale=None, decades=5,
                     outDirPath=None, figSize=(8, 6), interactive=None):

    # loop over all applicable instruments
    for instrument, filepath in sm.instrumentOutFilePaths(simulation, "stokesQ.fits"):
        # form the simulation/instrument name
        insname = "{}_{}".format(instrument.prefix(), instrument.name())

        # get the file paths for the frames/data cubes
        filepathI = instrument.outFilePaths("total.fits")[0]
        filepathQ = instrument.outFilePaths("stokesQ.fits")[0]
        filepathU = instrument.outFilePaths("stokesU.fits")[0]
        filepathV = instrument.outFilePaths("stokesV.fits")[0]

        # load datacubes with shape (nx, ny, nlambda)
        Is = sm.loadFits(filepathI)
        Qs = sm.loadFits(filepathQ)
        Us = sm.loadFits(filepathU)
        Vs = sm.loadFits(filepathV)

        # load the axes grids (assuming all files have the same axes)
        xgrid, ygrid, wavegrid = sm.getFitsAxes(filepathI)
        xmin = xgrid[0].value
        xmax = xgrid[-1].value
        ymin = ygrid[0].value
        ymax = ygrid[-1].value
        extent = (xmin, xmax, ymin, ymax)

        # determine binning configuration
        binX = binSize[0]
        orLenX = Is.shape[0]
        dropX =  orLenX % binX
        startX = dropX//2
        binY = binSize[1]
        orLenY = Is.shape[1]
        dropY = orLenY % binY
        startY = dropY//2

        # construct arrays with central bin positions in pixel coordinates
        posX = np.arange(startX - 0.5 + binX / 2.0, orLenX - dropX + startX - 0.5, binX)
        posY = np.arange(startY - 0.5 + binY / 2.0, orLenY - dropY + startY - 0.5, binY)

        # determine the appropriate wavelength index or indices
        if wavelength is None:
            indices = [0]
        elif wavelength == 'all':
            indices = range(Is.shape[2])
        else:
            if not isinstance(wavelength, (list,tuple)): wavelength = [wavelength]
            indices = instrument.wavelengthIndices(wavelength)

        # loop over all requested wavelength indices
        for index in indices:
            wave = wavegrid[index]
            wavename = "{:09.4f}".format(wave.to_value(sm.unit("micron")))
            wavelatex = r"$\lambda={:.4g}$".format(wave.value) + sm.latexForUnit(wave)

            # extract the corresponding frame, and transpose to (y,x) style for compatibility with legacy code
            I = Is[:,:,index].T.value
            Q = Qs[:,:,index].T.value
            U = Us[:,:,index].T.value
            V = Vs[:,:,index].T.value

            # perform the actual binning
            binnedI = np.zeros((len(posY),len(posX)))
            binnedQ = np.zeros((len(posY),len(posX)))
            binnedU = np.zeros((len(posY),len(posX)))
            binnedV = np.zeros((len(posY),len(posX)))
            for x in range(len(posX)):
                for y in range(len(posY)):
                    binnedI[y,x] = np.sum(I[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedQ[y,x] = np.sum(Q[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedU[y,x] = np.sum(U[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                    binnedV[y,x] = np.sum(V[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])

            # -----------------------------------------------------------------

            # plot a linear polarization map
            if plotLinMap:
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel("y" + sm.latexForUnit(ygrid), fontsize='large')

                # determine intensity range for the background image, ignoring pixels with outrageously high flux
                Ib = I.copy()
                highmask = Ib > 1e6*np.nanmedian(np.unique(Ib))
                vmax = np.nanmax(Ib[~highmask])
                Ib[highmask] = vmax
                vmin = vmax / 10**decades

                # plot the background image and the corresponding color bar
                normalizer = matplotlib.colors.LogNorm(vmin, vmax)
                cmap = plt.get_cmap('PuRd').copy()
                cmap.set_under('w')
                backPlot = ax.imshow(Ib, norm=normalizer, cmap=cmap, extent=extent,
                                     aspect='equal', interpolation='bicubic', origin='lower')
                cbarlabel = sm.latexForSpectralFlux(Is) + sm.latexForUnit(Is) + " @ " +  wavelatex
                plt.colorbar(backPlot, ax=ax).set_label(cbarlabel, fontsize='large')

                # compute the linear polarization degree
                degreeLD = np.sqrt(binnedQ**2 + binnedU**2)
                degreeLD[degreeLD>0] /= binnedI[degreeLD>0]

                # determine a characteristic 'high' degree of polarization in the frame
                # (this has to be done before degreeLD contains 'np.NaN')
                charDegree = np.percentile(degreeLD, 99.0)
                if not 0<charDegree<1:
                    charDegree = np.nanmax((np.nanmax(degreeLD), 0.0001))

                # remove pixels with minuscule polarization
                degreeLD[degreeLD<charDegree/50] = np.NaN

                # determine the scaling so that the longest arrows do not to overlap with neighboring arrows
                if degreeScale is None:
                    degreeScale = _roundUp(charDegree)
                lengthScale = 2.2 * degreeScale * max(float(len(posX))/figSize[0], float(len(posY))/figSize[1])
                key = "{:.3g}%".format(100 * degreeScale)

                # compute the polarization angle
                angle = 0.5 * np.arctan2(binnedU, binnedQ) # angle from North through East while looking at the sky

                # create the polarization vector arrays
                xPolarization = -degreeLD*np.sin(angle) #For angle = 0: North & x=0, For angle = 90deg: West & x=-1
                yPolarization =  degreeLD*np.cos(angle) #For angle = 0: North & y=1, For angle = 90deg: West & y=0

                # plot the vector field (scale positions to data coordinates)
                X,Y = np.meshgrid(xmin + posX * (xmax - xmin) / orLenX, ymin + posY * (ymax - ymin) / orLenY)
                quiverPlot = ax.quiver(X,Y, xPolarization, yPolarization, pivot='middle', units='inches',
                                        angles='xy', scale=lengthScale, scale_units='inches', headwidth=0,
                                        headlength=1, headaxislength=1, minlength=0.8, width=0.02)
                ax.quiverkey(quiverPlot, 0.85, 0.02, degreeScale, key, coordinates='axes', labelpos='E')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepath, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_pollinmap.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot a linear polarization degree map
            if plotDegMap:
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel("y" + sm.latexForUnit(ygrid), fontsize='large')

                # calculate polarization degree for each pixel, in percent
                # set degree to zero for pixels with very low intensity
                cutmask = I < ( np.nanmax(I[I < 1e6*np.nanmedian(np.unique(I))]) / 10**decades )
                degreeHD = np.sqrt(Q**2 + U**2)
                degreeHD[~cutmask] /= I[~cutmask]
                degreeHD[cutmask] = 0
                degreeHD *= 100

                # plot the image and the corresponding color bar
                vmax = degreeScale if degreeScale is not None else np.percentile(degreeHD, 99)
                normalizer = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
                backPlot = ax.imshow(degreeHD, norm=normalizer, cmap='plasma', extent=extent,
                                     aspect='equal', interpolation='bicubic', origin='lower')
                plt.colorbar(backPlot, ax=ax).set_label("Linear polarization degree (%)" + " @ " +  wavelatex,
                                                        fontsize='large')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepath, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_poldegmap.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot the y-axis averaged linear polarization degree
            if plotDegAvg:
                # construct the plot
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)
                degreeHD = np.sqrt(np.average(Q, axis = 0)**2 + np.average(U, axis=0)**2)
                degreeHD /= np.average(I,axis = 0)
                ax.plot(xgrid.value, degreeHD*100)
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(0, degreeScale)
                ax.set_title("{}   {}".format(insname, wavelatex), fontsize='large')
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel('Average linear polarization degree (%)', fontsize='large')

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepathI, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_poldegavg.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

            # -----------------------------------------------------------------

            # plot a circular polarization map
            if plotCirMap:
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figSize)

                # configure the axes
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel("x" + sm.latexForUnit(xgrid), fontsize='large')
                ax.set_ylabel("y" + sm.latexForUnit(ygrid), fontsize='large')

                # determine intensity range for the background image, ignoring pixels with outrageously high flux
                Ib = I.copy()
                highmask = Ib > 1e6*np.nanmedian(np.unique(Ib))
                vmax = np.nanmax(Ib[~highmask])
                Ib[highmask] = vmax
                vmin = vmax / 10**decades

                # plot the background image and the corresponding color bar
                normalizer = matplotlib.colors.LogNorm(vmin, vmax)
                cmap = plt.get_cmap('PuRd').copy()
                cmap.set_under('w')
                backPlot = ax.imshow(Ib, norm=normalizer, cmap=cmap, extent=extent,
                                     aspect='equal', interpolation='bicubic', origin='lower')
                cbarlabel = sm.latexForSpectralFlux(Is) + sm.latexForUnit(Is) + " @ " +  wavelatex
                plt.colorbar(backPlot, ax=ax).set_label(cbarlabel, fontsize='large')

                # compute the circular polarization degree
                degreeLD = binnedV.copy()
                degreeLD[binnedI>0] /= binnedI[binnedI>0]

                # determine the scaling and add legend
                if degreeScale is None:
                    degreeScale = _roundUp(np.percentile(np.abs(degreeLD), 99))
                lengthScale = 0.7/max(len(posX),len(posY))
                _circArrow(ax, 0.84-lengthScale/2, 0.01+lengthScale/2, lengthScale)
                key = r'$+{} \%$'.format(100*degreeScale)
                ax.text(0.85, 0.01+lengthScale/2, key, transform=ax.transAxes, ha='left', va='center')

                # actual plotting
                for x in range(len(posX)):
                    for y in range(len(posY)):
                        if np.isfinite(degreeLD[y,x]) and abs(degreeLD[y,x]) > degreeScale/50:
                            _circArrow(ax, posX[x]/orLenX, posY[y]/orLenY, degreeLD[y,x]/degreeScale*lengthScale)

                # if not in interactive mode, save the figure; otherwise leave it open
                if not ut.interactive(interactive):
                    saveFilePath = ut.savePath(filepath, ".pdf", outDirPath=outDirPath,
                                               outFileName = "{}_{}_polcirmap.pdf".format(insname, wavename))
                    plt.savefig(saveFilePath, bbox_inches='tight', pad_inches=0.25)
                    plt.close()
                    logging.info("Created {}".format(saveFilePath))

# -----------------------------------------------------------------

## This helper function rounds up a given number (away from zero) to one significant digit, in a logarithmic fashion.
# The \em slack argument allows a number just larger than the nearest limit to be clamped to that limit.
def _roundUp(x, slack = 0.01):
    x *= 1.-slack
    signX = np.sign(x)
    x = np.abs(x)
    dec = 10**(np.floor(np.log10(x)))
    mant = x/dec
    if mant <= 1: return signX*1*dec
    if mant <= 2: return signX*2*dec
    if mant <= 4: return signX*4*dec
    if mant <= 5: return signX*5*dec
    return signX*10*dec

## This helper function draws a circular polarization arrow in the given matplotlib axes.
def _circArrow(ax, posx, posy, size):
    # add the arc
    lw = 1
    oa = 33
    ax.add_artist(matplotlib.patches.Arc((posx,posy), size, size, theta1=oa, theta2=-oa,
                                         lw=lw, capstyle='round',
                                         color='k', zorder = 10, clip_on=False, transform=ax.transAxes))
    # add the arrow
    x = posx + size / 2 * np.cos(oa / 180 * np.pi)
    y = posy - abs(size) / 2 * np.sin(oa / 180 * np.pi)
    ax.add_artist(matplotlib.lines.Line2D((x - 0.5 * size / 2, x, x), (y, y, y - 0.5 * abs(size) / 2),
                                          lw=lw, solid_capstyle='round', dash_capstyle='round',
                                          color='k', zorder = 10, clip_on=False, transform=ax.transAxes))

# -----------------------------------------------------------------
