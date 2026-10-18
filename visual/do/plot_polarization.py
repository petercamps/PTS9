#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_polarization Plot polarization maps from SKIRT simulation output
#
# This script creates PDF polarization maps for the output of polarization-enabled SKIRT simulations that
# include instruments generating surface brightness frames or data cubes (*_total.fits and *_stokes*.fits files).
# Other simulations and instruments are silently ignored.
#
# The script takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory.
#  - \em plot (string): the type of plot to make (default is "linear"):
#    - "linmap": a linear polarization map
#    - "degmap": a linear polarization degree map
#    - "degavg": the y-axis averaged linear polarization degree for all x pixels
#    - "cirmap": a circular polarization map
#  - \em wavelength: the wavelength in micron for which to create the plot, or zero to create a plot for each
#                    wavelength in the data cube.
#  - bin: the number of pixels in each bin in both horizontal and vertical directions.
#  - \em dex (float): the number of decades to be included in the background intensity range (color bar); default is 5.
#
# In all cases, the plot files are placed next to the simulation output file(s) being handled. The file names include
# the simulation prefix, the instrument name, a wavelength indicator, and the ".pdf" filename extension.
#

# -----------------------------------------------------------------

def do( simDirPath : (str, "SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        plot : (str,"type of plot: linmap, degmap, degavg, or cirmap") = "linmap",
        wave: (float,"wavelength of the frame to be plotted; 0 means all frames") = 0,
        bin : (int,"number of pixels in a bin, in both x and y directions") = 7,
        dex : (float,"number of decades to be included in the background intensity range (color bar)") = 5,
        ) -> "plot polarization maps for the output of one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix) > 0 else None):
        vis.plotPolarization(sim,
                             plotLinMap=plot.lower().startswith("lin"),
                             plotDegMap=plot.lower().startswith("degm"),
                             plotDegAvg=plot.lower().startswith("dega"),
                             plotCirMap=plot.lower().startswith("cir"),
                             wavelength='all' if wave==0 else wave << sm.unit("micron"),
                             binSize=(bin,bin),
                             decades=dex)

# -----------------------------------------------------------------
