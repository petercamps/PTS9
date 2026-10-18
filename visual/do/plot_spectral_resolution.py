#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_spectral_resolution Plot the spectral resolution of a wavelength axis
#
# This script creates a plot of the spectral resolution \f$R=\lambda/\Delta\lambda\f$ of a wavelength grid
# loaded from a SKIRT stored table (".stab") or a SKIRT output file (".dat") that includes a wavelength axis.
#
# Provide the name or (relative or absolute) file path of a SKIRT stored table file or a SKIRT output file
# as the single required argument; a SKIRT output file written as a bundle in an HDF5 file can be specified as
# <tt>\<filepath\>:\<suite\>/\<bundle\></tt>. Use optional arguments to further configure the plot:
#  - \em wmin (float): smallest wavelength on the horizontal axis; default is 0.1
#  - \em wmax (float): largest wavelength on the horizontal axis; default is 1000
#  - \em unit (string): unit of the wavelength limits (any spectral unit); default is micron
#  - \em dex (float): if specified, the number of decades to be plotted on the vertical axis; default is 3
#  - \em title (string): the title used in the plot legend; default is the name of the input file
#
# The plot file has the same name as the input file but with the ".pdf" filename extension and is placed
# in the current working directory.
#

# -----------------------------------------------------------------

def do( filepath : (str,"name or path of SKIRT stored table or text column file"),
        wmin : (float,"smallest wavelength on the horizontal axis") = 0.1,
        wmax : (float,"largest wavelength on the horizontal axis") = 1000.,
        unit: (str,"unit of the wavelength limits (any spectral unit)") = "micron",
        dex : (float,"number of decades to be plotted on the vertical axis") = 3,
        title : (str,"title used in plot legend") = "",
        ) -> "plot the spectral resolution of a wavelength axis":

    import pts.simulation as sm
    import pts.visual as vis

    vis.plotSpectralResolution(inFilePath=filepath, minWavelength=wmin*sm.unit(unit), maxWavelength=wmax*sm.unit(unit),
                               decades=dex, title=title)

# -----------------------------------------------------------------
