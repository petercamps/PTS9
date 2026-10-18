#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_convergence Plot convergence density cuts from one or more SKIRT simulations
#
# This script creates plots of the planar cuts through the input and gridded medium density produced by the
# convergence cuts probe in a SKIRT simulation. This allows a visual comparison between the medium density distribution
# of the input model and the corresponding spatially discretized version used in the simulation.
# If the simulation does not include any ConvergenceCutsProbe instances, the script does nothing.
#
# A single plot is produced for each medium (dust, electrons, gas) including panels for both input and gridded density
# (rows) and for all orientations (columns) depending on the symmetry of the model. All panels share the same color
# scale. The dynamic range of the logarithmic color scale is specified by \em decades.
#
# The script takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#  - \em dex (float): if specified, the number of decades to be included in the density range (color bar); default is 5
#
# In all cases, the plot file is placed next to the simulation output file(s) being handled. The filename includes the
# simulation prefix, the probe name, and the medium type indicator, and has the ".pdf" filename extension.
#

# -----------------------------------------------------------------

def do( simDirPath : (str, "SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        dex : (float,"number of decades to be included in the density range (color bar)") = 5,
        ) -> "plot convergence density cuts from one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix) > 0 else None):
        vis.plotConvergenceCuts(sim, decades=dex)

# ----------------------------------------------------------------------
