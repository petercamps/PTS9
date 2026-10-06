#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.test.do.endorse_functional Endorse the current output of standard functional SKIRT tests
#
# This script "endorses" the current output of all or a selection of the standard functional SKIRT test cases
# by replacing the contents of the "ref" directory by the contents of the "out" directory.

# \note This script destroys the current reference output for the specified test cases; USE WITH CARE!
#
# The test case definitions are expected to reside in the \c SKIRT/Functional10 directory hierarchy.
# See the pts.test.functional module for more information.
#
# The script takes a single positional string argument, which can be one of the following:
#  - "." (a single period): clean all test cases in the standard suite.
#  - "testcase" (the name of a test case directory): clean all test cases with that name.
#  - "subsuite" (the name of an intermediate directory in the hierarchy): clean the test cases in all sub-suites
#    with that name.
#  - "parentsubsuite/testcase" or "parentsubsuite/subsuite": clean the indicated test case(s) or sub-suite(s)
#    that reside immediately inside the indicated parent sub-suite; this can disambiguate items with the same name.
#

# -----------------------------------------------------------------

def do( subSuite : (str,"name of sub-suite or test case to be endorsed or '.' to endorse all"),
        ) -> "endorse the output of all or a selection of the standard functional SKIRT tests":

    import logging
    import pts.test

    suite = pts.test.SkirtTestSuite(subSuite=subSuite)
    logging.warning("** This will replace the reference output for {} test cases by the current test output **" \
                    .format(suite.size()))
    response = input("Enter 'YES' to overwrite test cases ({}): ".format(suite.size()))
    if response == "YES":
        logging.info("Replacing the reference output for {} test cases...".format(suite.size()))
        suite.endorse()

# -----------------------------------------------------------------
