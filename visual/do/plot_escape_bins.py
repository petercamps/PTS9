#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_escape_bins Plot the flux escaping from the model as a function of direction
#
# This script plots the flux recorded by the escape instruments of one or more SKIRT simulations as a function of the
# escape direction, at a single wavelength or convolved with a broadband: a step plot versus inclination for
# inclination bins, a map versus azimuth and inclination for inclination and azimuth bins, and a projection of the
# sky for HEALPix bins (see the pts.visual.plotEscapeBins() function). It takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#  - \em instr (string): the name of the escape instrument to handle; by default handles all escape instruments
#  - \em wave (float): the wavelength at which to plot the flux; by default (or if zero), the wavelength where the
#                      total escaping flux is largest
#  - \em unit (string): the unit of the wavelength (any spectral unit); default is micron
#  - \em band (string): the name of a builtin broadband (see the pts.band package) with which to convolve the flux
#                       instead; cannot be combined with a nonzero wavelength
#  - \em component (string): the flux component to plot, such as "total" (the default), "transparent",
#                            "primaryscattered", or "secondarydirect"
#  - \em dex (float): the number of decades in the range of the vertical axis or color scale; default is 3
#  - \em projection (string): the projection for HEALPix bins, Mollweide (the default) or HammerAitoff
#
# The plot files are placed next to the output files of each escape instrument; their names end with "_escape.pdf".
#

# -----------------------------------------------------------------

def do( simDirPath : (str,"SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        instr : (str,"escape instrument name") = "",
        wave : (float,"wavelength at which to plot the flux (default: where the escaping flux is largest)") = 0,
        unit: (str,"unit of the wavelength (any spectral unit)") = "micron",
        band : (str,"name of a builtin broadband with which to convolve the flux") = "",
        component : (str,"flux component to plot") = "total",
        dex : (float,"number of decades in the range of the vertical axis or color scale") = 3,
        projection : (str,"projection for HEALPix bins (Mollweide/HammerAitoff)") = "Mollweide",
        ) -> "plot the flux escaping from the model as a function of direction":

    import pts.band as bnd
    import pts.simulation as sm
    import pts.utils as ut
    import pts.visual as vis

    if wave != 0 and len(band) > 0:
        raise ut.UserError("wave and band arguments cannot both be specified")
    wavelength = wave * sm.unit(unit) if wave != 0 else None
    try: broadband = bnd.builtinBand(band) if len(band) > 0 else None
    except ValueError as error: raise ut.UserError(str(error))

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix)>0 else None):
        instruments = [ ins for ins in sim.instruments() if ins.type() == "EscapeInstrument" ]
        if len(instr) > 0:
            instruments = [ ins for ins in instruments if ins.name() == instr ]
            if len(instruments) == 0:
                raise ut.UserError("simulation '{}' does not have escape instrument '{}'".format(sim.prefix(), instr))
        vis.plotEscapeBins(instruments, wavelength=wavelength, band=broadband, component=component, decades=dex,
                           projection=projection)

# ----------------------------------------------------------------------
