#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skiupgrade.do.upgrade_functional_tests Upgrade ski files for all or a subsuite of functional tests
#
# This script upgrades the ski files for a selection of the standard functional test cases for SKIRT 9.
# The test case definitions are expected to reside in the \c SKIRT/Functional10 directory hierarchy.
#
# The script takes a single positional string argument, which can be one of the following:
#  - "." (a single period): upgrade all test cases in the standard suite.
#  - "testcase" (the name of a test case directory): upgrade all test cases with that name.
#  - "subsuite" (the name of an intermediate directory in the hierarchy): upgrade the test cases in all sub-suites
#    with that name.
#  - "parentsubsuite/testcase" or "parentsubsuite/subsuite": upgrade the indicated test case(s) or sub-suite(s)
#    that reside immediately inside the indicated parent sub-suite; this can disambiguate items with the same name.
#
# For each of the selected test cases, the script performs as follows:
# - if the ski file is up to date with the latest version of SKIRT 9, an informational message is logged and the
#   file is otherwise ignored.
# - if the ski file needs upgrading, it is first copied to a backup version with a filename including a time stamp,
#   it is then replaced by a version upgraded to the latest version of SKIRT 9, and a warning message is logged.
#
# See the pts.test.functional module for more information about test cases and test (sub-)suites.
# See the pts.skiupgrade.skiupgrade module for more information on the ski file upgrade process.
#

# -----------------------------------------------------------------

def do( subSuite : (str,"name of sub-suite or test case to be performed or '.' to perform all"),
        ) -> "upgrade ski files for functional tests to the latest version of SKIRT 9":

    import pts.skiupgrade
    import pts.test

    for skipath in pts.test.SkirtTestSuite(subSuite=subSuite).skiPaths():
        pts.skiupgrade.upgradeSkiFile(skipath)

# -----------------------------------------------------------------
