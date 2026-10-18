#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.do.make_wavelength_movie Create a movie that runs through all wavelengths in the SKIRT simulation output
#
# This script creates a movie in MP4 format for the output of the specified simulation. The movie combines the SEDs
# (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments, running through all
# wavelengths in the simulation.
#
# The script takes the following arguments:
#  - \em simDirPath (positional string argument): the path to the SKIRT simulation output directory,
#                                                 or "." for the current directory; output bundles in an HDF5 file are
#                                                 located through the log file, or by appending the HDF5 file and suite
#                                                 to the path, as for the -o option of SKIRT.
#  - \em prefix (string): the prefix of the simulation to handle; by default handles all simulations in the directory.
#  - \em instruments (string): up to three comma-separated instrument names; by default handles the first three
#                              instruments in the simulated configuration.
#  - \em percentile (float): the percentile value in range [0,100] used to determine the maximum surface brightness in
#                   the images. The default value is 100, so that the largest surface brightness value in the frame(s)
#                   determines the maximum value. A smaller percentile value such as 99 could be used to exclude the
#                   very largest surface brightness values (for example, the direct light from a point source).
#  - \em dex (float): the number of decades in the surface brightness range in the images. The default value is 5.
#  - \em renormalize (int): a flag selecting one of two image scaling options. When zero (the default value), the
#                   surface brightness range is determined for the complete data cube and kept constant for all frames.
#                   As a result, image frames might be completely dark at wavelengths with very low overall luminosity.
#                   When nonzero, the surface brightness range is determined for each frame separately. This allows to
#                   see the spatial structure of the surface brightness even at wavelengths with very low luminosity.
#  - \em rate:      the frame rate of the movie, in frames per second. The default value is 7.
#
# By default, the movie is saved in the output directory of the first instrument, using a name starting with the
# corresponding simulation prefix and ending with ".mp4". This can be overridden with the out* arguments as described
# for the pts.utils.savePath() function.
#

# -----------------------------------------------------------------

def do( simDirPath : (str, "SKIRT simulation output directory"),
        prefix : (str,"SKIRT simulation prefix") = "",
        instruments : (str,"up to three comma-separated SKIRT instrument names") = "",
        percentile : (float,"percentile in range [0,100] used to determine maximum surface brightness") = 100,
        dex : (float,"number of decades in the surface brightness range") = 5,
        renormalize : (int,"if nonzero, the surface brightness range is determined for each frame separately") = 0,
        rate : (int,"frame rate of the movie in frames per second") = 7,
        ) -> "create a movie that runs through all wavelengths in the SKIRT simulation output":

    import pts.simulation as sm
    import pts.utils as ut
    import pts.visual as vis

    # private function to retrieve the instrument with a specified name from a simulation
    def instrumentWithName(sim, instrname):
        instrlist = [ instr for instr in sim.instruments() if instr.name()==instrname ]
        if len(instrlist) == 0:
            raise ut.UserError("simulation '{}' does not have instrument '{}'".format(sim.prefix(), instrname))
        return instrlist[0]

    # loop over the simulations to be handled
    for sim in sm.createSimulations(simDirPath, prefix if len(prefix) > 0 else None):

        # get a list of the requested instruments (or simply specify the complete simulation)
        simspec = sim
        if len(instruments) > 0:
            simspec = [ instrumentWithName(sim, instrname) for instrname in instruments.split(',') ]

        # create the movie
        vis.makeWavelengthMovie(simspec, maxPercentile=percentile, decades=dex, renormalize=renormalize!=0, rate=rate)

# ----------------------------------------------------------------------
