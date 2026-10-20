#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.do.convert_dust_import_file Convert a SKIRT 9 dust import file to SKIRT 10
#
# This script converts an import file for a dust medium or an imported geometry that used the \em importMetallicity
# or \em importTemperature options of SKIRT 9. SKIRT 10 no longer offers these options, so that the import file must
# provide the effective mass or density, without metallicity and temperature columns. The script multiplies the
# values in the mass or density column(s) by the metallicity, sets them to zero for the entities with a temperature
# above the maximum temperature (SKIRT ignores smoothed particles with zero mass), and removes the metallicity and
# temperature columns. See the pts.simulation.text.convertDustImportFile() function for more information.
#
# The script takes the following arguments:
#  - \em inFilePath (positional string argument): the file path of the SKIRT 9 import file; a column text file
#    in an HDF5 file is not supported.
#  - \em outFilePath (positional string argument): the file path of the converted import file, which must not yet
#    exist.
#  - \em massColumns (positional string argument): a comma-separated list of the one-based indices of the columns
#    holding the mass or density, usually a single index.
#  - \em metallicityColumn (integer): the one-based index of the metallicity column, or zero if there is none
#    (the default).
#  - \em temperatureColumn (integer): the one-based index of the temperature column, or zero if there is none
#    (the default).
#  - \em maxTemperature (float): the value of the \em maxTemperature option in K, or zero if there is no
#    temperature cutoff (the default).
#
# For example, the following command converts a particle file with the mass in column 5, the metallicity in column 6,
# and the temperature in column 7, for a temperature cutoff at 8000 K:
#
#     pts convert_dust_import_file particles.txt dust.txt 5 --metallicityColumn 6 --temperatureColumn 7
#         --maxTemperature 8000
#

# -----------------------------------------------------------------

def do( inFilePath : (str,"file path of the SKIRT 9 import file"),
        outFilePath : (str,"file path of the converted import file"),
        massColumns : (str,"comma-separated one-based indices of the mass or density columns"),
        metallicityColumn : (int,"one-based index of the metallicity column, or zero if there is none") = 0,
        temperatureColumn : (int,"one-based index of the temperature column, or zero if there is none") = 0,
        maxTemperature : (float,"maximum temperature in K, or zero if there is no temperature cutoff") = 0.,
        ) -> "convert a SKIRT 9 import file for a dust medium or imported geometry to SKIRT 10":

    import pts.simulation as sm
    import pts.utils as ut

    try:
        columns = [ int(column) for column in massColumns.split(",") ]
    except ValueError:
        raise ut.UserError("Mass or density column indices must be integers: {}".format(massColumns))
    sm.convertDustImportFile(inFilePath, outFilePath, columns, metallicityColumn, temperatureColumn, maxTemperature)

# -----------------------------------------------------------------
