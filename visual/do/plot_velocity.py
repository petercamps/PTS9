#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.plot_velocity Plot planar velocity cuts or projections from one or more SKIRT simulations
#
# This script creates PDF vector maps for the planar velocity cuts or projections produced by one of the relevant
# probes in a SKIRT simulation, allowing a visual evaluation of the velocity field.
#
# The script supports VelocityProbe, ImportedSourceVelocityProbe, and/or ImportedMediumVelocityProbe
# instances with an associated probe form that produces a planar cut (DefaultCutsForm, PlanarCutsForm) or planar
# projection (ParallelProjectionForm, AllSkyProjectionForm). If the simulation does not include any supported probes,
# the script does nothing.
#
# The script takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory
#  - \em bin: the number of pixels in each bin in both horizontal and vertical directions.
#
# In all cases, the plot file is placed next to the simulation output file(s) being handled. The plot file has the same
# name as the SKIRT output file from which it is derived but with the ".pdf" filename extension instead of ".fits".
#

# -----------------------------------------------------------------

def do( simDirPath : (str, "SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        bin : (int,"number of pixels in a bin, in both x and y directions") = 32,
        ) -> "plot planar velocity cuts or projections from one or more SKIRT simulations":

    import pts.simulation as sm
    import pts.visual as vis

    for sim in sm.createSimulations(simDirPath, prefix if len(prefix) > 0 else None):
        vis.plotVectorCuts(sim, ("VelocityProbe", "ImportedSourceVelocityProbe", "ImportedMediumVelocityProbe"),
                           binSize=(bin,bin))

# ----------------------------------------------------------------------
