#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_density Plot planar density cuts or projections from one or more SKIRT simulations
#
# This script creates plots of the density cuts or projections produced by one of the relevant probes
# in a SKIRT simulation. The script supports DensityProbe, ImportedSourceDensityProbe, and/or ImportedMediumDensityProbe
# instances with an associated probe form that produces a planar cut (DefaultCutsForm, PlanarCutsForm) or planar
# projection (ParallelProjectionForm, AllSkyProjectionForm). For cuts, the plotted quantity is a mass or number density;
# for projections, the plotted quantity is a surface mass density or a column density.
# If the simulation does not include any supported probes, the script does nothing.
#
# If the output for a given medium component or medium type (dust, gas, or electrons) includes two or three files
# with names that differ only by orientation labels ("_xy", "_xz", and/or "_yz"), the density maps for these
# files are included in a single plot and share the same scale.
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
# simulation prefix, the probe name, and the medium component or type indicator, and has the ".pdf" filename extension.
#

# -----------------------------------------------------------------

def do( simDirPath : (str, "SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        dex : (float,"number of decades to be included in the density range (color bar)") = 5,
        ) -> "plot planar density cuts or projections from one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix) > 0 else None):
        vis.plotScalarCuts(sim, probeTypes=("DensityProbe", "ImportedSourceDensityProbe", "ImportedMediumDensityProbe"),
                           decades=dex)

# ----------------------------------------------------------------------
