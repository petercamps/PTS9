#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_seds Plot the SEDs produced by one or more SKIRT simulations
#
# This script creates plots of the SEDs produced by one or more SKIRT simulations. It takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#  - \em instr (string): the name of the instrument to handle; by default handles all instruments in the simulation
#  - \em bins (string): comma-separated zero-based indices of the escape bins to be plotted for escape instruments
#                       with the SEDFiles output style; by default plots all bins
#  - \em wmin (float): smallest wavelength on the horizontal axis; default is 0.1
#  - \em wmax (float): largest wavelength on the horizontal axis; default is 1000
#  - \em unit (string): unit of the wavelength limits (any spectral unit); default is micron
#  - \em dex (float): if specified, the number of decades to be plotted on the vertical axis; default is 5
#
# If there is only a single instrument (because the instrument name was specified or because the simulation only has
# one instrument), and the corresponding output file includes the individual flux components, the script creates
# a plot showing all flux components. Otherwise, it creates a plot showing the total flux for all instruments.
# An escape instrument contributes an SED for each of its (selected) bins, labelled with the bin index and the
# directions covered by the bin; if there are more than six such bins, they are colored along a color map instead.
#
# In all cases, the plot file is placed next to the "*_sed.dat" file(s) being represented. The filename starts
# with the simulation prefix and ends with "sed.pdf".
#

# -----------------------------------------------------------------

def do( simDirPath : (str,"SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        instr : (str,"instrument name") = "",
        bins : (str,"comma-separated indices of the escape bins to be plotted (default: all)") = "",
        wmin : (float,"smallest wavelength on the horizontal axis") = 0.1,
        wmax : (float,"largest wavelength on the horizontal axis") = 1000.,
        unit: (str,"unit of the wavelength limits (any spectral unit)") = "micron",
        dex : (float,"number of decades to be plotted on the vertical axis") = 5,
        ) -> "plot the SEDs produced by one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.utils as ut
    import pts.visual as vis

    # parse the bin indices
    binIndices = None
    if len(bins) > 0:
        try: binIndices = [ int(index) for index in bins.split(",") ]
        except ValueError: raise ut.UserError("bins argument must be a comma-separated list of integers")

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix)>0 else None):
        if len(instr)>0:
            instrument = [ ins for ins in sim.instruments() if ins.name() == instr ]
            vis.plotSeds(instrument, minWavelength=wmin*sm.unit(unit), maxWavelength=wmax*sm.unit(unit), decades=dex,
                         bins=binIndices)
        else:
            vis.plotSeds(sim, minWavelength=wmin*sm.unit(unit), maxWavelength=wmax*sm.unit(unit), decades=dex,
                         bins=binIndices)

# ----------------------------------------------------------------------
