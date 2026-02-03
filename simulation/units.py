#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.units Handling SKIRT units
#
# This module offers functions for converting SKIRT units to astropy units, and for handling conversions that
# may not be trivial to perform using pure astropy.
#

# -----------------------------------------------------------------

import astropy.constants as const
import astropy.units as u
import warnings

# -----------------------------------------------------------------

# globally enable automatic conversion between photon wavelength, frequency, and energy units
u.set_enabled_equivalencies(u.spectral())

# -----------------------------------------------------------------

## This function returns an astropy unit for any of the following "unit-like" inputs:
#  - a string representing a valid SKIRT unit: returns the corresponding astropy unit.
#  - a string representing a valid astropy unit in generic format: returns the corresponding astropy unit.
#  - an astropy quantity: returns the unit of the quantity.
#  - an astropy unit: simply returns the unit itself.
#
# While the function is mainly intended for converting SKIRT units to astropy units, supporting the other
# input types is useful when the function is invoked directly or indirectly from some other context.
#
def unit(unitlike):
    if isinstance(unitlike, str):
        # handle SKIRT-specific pressure equivalence: 1 K/m3 <==> k_B Pa
        if unitlike=="K/m3":
            unitlike = "1.3806488e-23 Pa"
        # parse string ignoring warnings about multiple divisions, as in "W/m2/sr"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=u.UnitsWarning)
            # astropy does not support multiple division for units starting with "1/" or "/", as in
            # "1/s/keV" or "/s/keV"; replace the leading "1" by an arbitrary unit and then remove it again
            if unitlike.startswith("1/"):
                return u.Unit("A" + unitlike[1:]) / u.Unit("A")
            if unitlike.startswith("/"):
                return u.Unit("A" + unitlike) / u.Unit("A")
            # in other cases, directly use the astropy parser
            return u.Unit(unitlike)
    if isinstance(unitlike, u.Quantity):
        return unitlike.unit
    if isinstance(unitlike, u.UnitBase):
        return unitlike
    raise ValueError("Unsupported unit-like type: {}".format(type(unitlike)))

# -----------------------------------------------------------------

## This function accepts a flux density, surface brightness, spectral radiance, or spectral luminosity
# (generically called "flux" for the purposes of this function), converts this flux from an arbitrary unit
# in any "flavor" (neutral, per wavelength, per frequency, or per energy) to an equivalent unit in the
# specified flavor, and returns the result.
#
# Both the wavelength and the flux must be astropy quantities in one of the following shape combinations:
#  - The wavelength is a scalar: this wavelength is used to convert all values in the flux argument, which can be a
#    scalar or an array of arbitrary shape (e.g. a 2D data frame) containing flux values at the same wavelength.
#  - The wavelength is a 1D array: the length of this array must match the last (or only) axis of the flux argument
#    which can be, for example an 1D spectrum or a 3D data cube. The wavelengths are broad-casted to the other axes.
#
# The target flavor can be specified in one of two ways:
#  - As one of the strings "neutral", "wavelength", "frequency", or "energy": the result will have the unit that
#    follows from the conversion calculation and thus depends on the incoming units of both flux and wavelength.
#    This option is useful in case the incoming flux type (density, brightness, luminosity) is not known a priori.
#  - As an astropy unit instance (such as the return value of the sm.unit() function): the flavor is derived from
#    this unit, and the result is explicitly converted to it before being returned. The unit must be compatible
#    with the incoming flux type (density, brightness, luminosity), after flavor conversion.
#
def convertToFlavor(wavelength, flux, flavor):
    # make a copy because the operations below will happen in place
    flux = flux.copy()

    # determine the conversion scheme: first digit = input flavor, second digit = output flavor
    scheme =  _flavor(flux.unit) * 10 + _flavor(flavor)

    # shortcuts for constants
    c = const.c
    hc2 = (const.h * const.c)**2

    # apply the conversion appropriate for each scheme (schemes 11, 22, 33, and 44 require no conversion)
    if scheme==12:            # neutral to wavelength
        flux /= wavelength
    elif scheme==13:          # neutral to frequency
        flux *= wavelength / c
    elif scheme==14:          # neutral to energy
        flux *= wavelength**2 / hc2

    elif scheme==21:          # wavelength to neutral
        flux *= wavelength
    elif scheme==23:          # wavelength to frequency
        flux *= wavelength**2 / c
    elif scheme==24:          # wavelength to energy
        flux *= wavelength**3 / hc2

    elif scheme==31:          # frequency to neutral
        flux *= c / wavelength
    elif scheme==32:          # frequency to wavelength
        flux *= c / wavelength**2
    elif scheme==34:          # frequency to energy
        flux *= wavelength * (c / hc2)

    elif scheme==41:          # energy to neutral
        flux *= hc2 / wavelength**2
    elif scheme==42:          # energy to wavelength
        flux *= hc2 / wavelength**3
    elif scheme==43:          # energy to frequency
        flux *= (hc2 / c) / wavelength

    # if the output flavor is specified as an explicit unit, convert to that unit
    if isinstance(flavor, u.UnitBase):
        flux <<= flavor
    return flux

## This helper function returns a numeric code corresponding to the specified flavor or unit:
#   - 1 for the string "neutral" and for units that have a neutral flavor.
#   - 2 for the string "wavelength" and for units that have a per wavelength flavor.
#   - 3 for the string "frequency" and for units that have a per frequency flavor.
#   - 4 for the string "energy" and for units that have a per energy flavor.
def _flavor(flavorunit):
    if isinstance(flavorunit, str):
        if flavorunit == "neutral": return 1
        if flavorunit == "wavelength": return 2
        if flavorunit == "frequency": return 3
        if flavorunit == "energy": return 4
        raise ValueError("Invalid flavor specification: '{}'".format(flavorunit))

    if flavorunit.is_equivalent(( unit("W"), unit("W/m2"), unit("W/m2/sr") )): return 1
    if flavorunit.is_equivalent(( unit("W/m"), unit("W/m2/m"), unit("W/m2/sr/m") )): return 2
    if flavorunit.is_equivalent(( unit("W/Hz"), unit("W/m2/Hz"), unit("W/m2/sr/Hz") )): return 3
    if flavorunit.is_equivalent(( unit("1/s/J"), unit("1/s/m2/J"), unit("1/s/m2/sr/J") )): return 4
    raise ValueError("Not a flux-like unit: '{}'".format(flavorunit))

# -----------------------------------------------------------------

## This function returns a latex-formatted string representation for the "unit-like" input,
# enclosed in square brackets and preceded by a short space. The input argument
# is interpreted as described for the unit() function in this module.
def latexForUnit(unitlike):
    latex = unit(unitlike).to_string("latex_inline")
    # work around for bug in astropy 3.1.2: put parentheses around LaTeX for arcsec if it carries an exponent
    latex = latex.replace(r"{}^{\prime\prime}^",r"({}^{\prime\prime})^")
    # use "1" for dimensionless quantities instead of the empty string
    if latex==r"$\mathrm{}$": latex = r"$\mathrm{1}$"
    return r"$\,$[" + latex + r"]"

## This function returns a latex-formatted string representation for wavelength (\f$\lambda\f$),
# frequency (\f$\nu\f$), or energy (\f$E\f$), depending on the units represented by the "unit-like" input
# (the latex string does not include the unit itself). The input argument is interpreted as described
# for the unit() function in this module. If no match is found, the function returns an empty string.
def latexForWavelength(unitlike):
    un = unit(unitlike)
    with u.set_enabled_equivalencies([]):
        if un.is_equivalent(unit("m")): return r"$\lambda$"
        if un.is_equivalent(unit("Hz")): return r"$\nu$"
        if un.is_equivalent(unit("J")): return r"$E$"
    return ""

## This function returns True if the units represented by the "unit-like" input are those for
# frequency (\f$\nu\f$) or energy (\f$E\f$), and False if not (including for wavelength (\f$\lambda\f$).
# The input argument is interpreted as described for the unit() function in this module.
def hasReverseWavelengthOrder(unitlike):
    un = unit(unitlike)
    with u.set_enabled_equivalencies([]):
        if un.is_equivalent(unit("Hz")): return True
        if un.is_equivalent(unit("J")): return True
    return False

## This function returns a latex-formatted string representation for a flux density (\f$\lambda F_\lambda\f$,
# \f$F_\lambda\f$, \f$F_\nu\f$, or \f$F_E\f$) or a surface brightness (\f$\lambda f_\lambda\f$, \f$f_\lambda\f$,
# \f$f_\nu\f$, or \f$f_E\f$) that has the units represented by the "unit-like" input (the latex string does not
# include the unit itself). The input argument is interpreted as described for the unit() function in this module.
# If no match is found, the function returns an empty string.
def latexForSpectralFlux(unitlike):
    un = unit(unitlike)
    with u.set_enabled_equivalencies([]):
        if un.is_equivalent(unit("W/m2")): return r"$\lambda\,F_\lambda$"
        if un.is_equivalent(unit("W/m2/m")): return r"$F_\lambda$"
        if un.is_equivalent(unit("W/m2/Hz")): return r"$F_\nu$"
        if un.is_equivalent(unit("1/s/m2/J")): return r"$F_E$"
        if un.is_equivalent(unit("W/m2/sr")): return r"$\lambda\,f_\lambda$"
        if un.is_equivalent(unit("W/m2/sr/m")): return r"$f_\lambda$"
        if un.is_equivalent(unit("W/m2/sr/Hz")): return r"$f_\nu$"
        if un.is_equivalent(unit("1/s/m2/sr/J")): return r"$f_E$"
    return ""

## This function returns a latex-formatted string representation for the mean intensity or spectral radiance
# (\f$\lambda J_\lambda\f$, \f$J_\lambda\f$, \f$J_\nu\f$, or \f$J_E\f$) that has the units represented by the
# "unit-like" input (the latex string does not include the unit itself). The input argument is interpreted as
# described for the unit() function in this module. If no match is found, the function returns an empty string.
def latexForSpectralRadiance(unitlike):
    un = unit(unitlike)
    with u.set_enabled_equivalencies([]):
        if un.is_equivalent(unit("W/m2/sr")): return r"$\lambda\,J_\lambda$"
        if un.is_equivalent(unit("W/m2/sr/m")): return r"$J_\lambda$"
        if un.is_equivalent(unit("W/m2/sr/Hz")): return r"$J_\nu$"
        if un.is_equivalent(unit("1/s/m2/sr/J")): return r"$J_E$"
    return ""

## This function returns a latex-formatted string representation for the spectral luminosity
# (\f$\lambda L_\lambda\f$, \f$L_\lambda\f$, \f$L_\nu\f$, or \f$L_E\f$) that has the units represented by the
# "unit-like" input (the latex string does not include the unit itself). The input argument is interpreted as
# described for the unit() function in this module. If no match is found, the function returns an empty string.
def latexForSpectralLuminosity(unitlike):
    un = unit(unitlike)
    with u.set_enabled_equivalencies([]):
        if un.is_equivalent(unit("W")): return r"$\lambda\,L_\lambda$"
        if un.is_equivalent(unit("W/m")): return r"$L_\lambda$"
        if un.is_equivalent(unit("W/Hz")): return r"$L_\nu$"
        if un.is_equivalent(unit("1/s/J")): return r"$L_E$"
    return ""

## This function returns the result of latexForWavelength() + latexForUnit() with the same argument.
def latexForWavelengthWithUnit(unitlike):
    return latexForWavelength(unitlike) + latexForUnit(unitlike)

## This function returns the result of latexForSpectralFlux() + latexForUnit() with the same argument.
def latexForSpectralFluxWithUnit(unitlike):
    return latexForSpectralFlux(unitlike) + latexForUnit(unitlike)

## This function returns the result of latexForSpectralRadiance() + latexForUnit() with the same argument.
def latexForSpectralRadianceWithUnit(unitlike):
    return latexForSpectralRadiance(unitlike) + latexForUnit(unitlike)

## This function returns the result of latexForSpectralLuminosity() + latexForUnit() with the same argument.
def latexForSpectralLuminosityWithUnit(unitlike):
    return latexForSpectralLuminosity(unitlike) + latexForUnit(unitlike)

# -----------------------------------------------------------------
