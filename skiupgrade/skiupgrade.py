#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skiupgrade.skiupgrade Contains the upgradeSkiFile function for upgrading SKIRT parameter files
#
# The upgradeSkiFile function in this module allows upgrading SKIRT parameter files (\em ski files)
# between versions of SKIRT 9.

# -----------------------------------------------------------------

import logging
import pts.simulation as sm
import pts.utils as ut

# -----------------------------------------------------------------

## This function upgrades the specified SKIRT 9 parameter file (\em ski file) so that it becomes appropriate for the
# latest SKIRT 9 version. The upgrade process supports all ski files created by the SKIRT project version 9 tools
# (including the SKIRT command line Q&A and the graphical MakeUp wizard) since SKIRT 9 was publicly released.
# Ski files created for older SKIRT versions (such as SKIRT 7 and 8) cannot be upgraded to SKIRT 9 automatically.
# The function logs an appropriate error message when an unsupported file is specified.
#
# The function accepts three arguments:
# - inpath: the absolute or relative path to the ski file to be handled; the filename extension should be ".ski"
# - backup: if the input ski file needs upgrading and \em backup is True (the default), a copy of the input ski file
#           is created with a filename including a time stamp and ending with "_backupski.xml"
# - replace: if the input ski file needs upgrading and \em replace is True (the default), the input ski file is
#            overwritten by the upgraded version; otherwise it is saved as a new file using a similar filename ending
#            with "_upgradedski.xml"
#
def upgradeSkiFile(inpath, *, backup=True, replace=True):
    # load the ski file
    inpath = ut.absPath(inpath)
    try:
        ski = sm.SkiFile(inpath)
    except SyntaxError:
        logging.error("File does not contain well-formed XML: {}".format(inpath))
        return

    # verify the ski file format version
    try:
        version = ski.getStringAttribute("/skirt-simulation-hierarchy", "format")
    except ValueError:
        logging.error("XML file does not have ski file format: {}".format(inpath))
        return
    if version!="9":
        logging.error("Ski file is older than version 9: {}".format(inpath))
        return

    # perform the upgrade in the XML tree in memory, keeping track of whether the contents has actually changed
    changed = False
    for condition,templates in _getUpgradeDefinitions():
        changed |= ski.transformIf(condition, templates)

    # save the upgraded version if needed
    if changed:
        if backup:
            inpath.rename(inpath.with_name(inpath.stem + "_" + ut.timestamp() + "_backupski.xml"))
        if replace:
            ski.saveTo(inpath)
            logging.warning("Ski file UPGRADED:  {}".format(inpath))
        else:
            outpath = inpath.with_name(inpath.stem + "_upgradedski.xml")
            ski.saveTo(outpath)
            logging.warning("Ski file UPGRADED:  {} --> {}".format(inpath, outpath.name))
    else:
        logging.info("Ski file unchanged: {}".format(inpath))

# -----------------------------------------------------------------

## This private function returns a sequence of 2-tuples, each defining the XPath condition and XSLT template
# for a single modification to the ski file format.
#
# Using XSLT is extremely flexible, but unfortunately the XSLT language is fairly obscure an thus has a steep
# learning curve. To alleviate this problem to some extent, we provide a set of functions that generate
# upgrade definitions for specific types of changes, such as, for example, changing the name of a property.
# That way, as soon as the XSLT sheet for a particular type of change has been developed, it can be more easily
# reused for other, similar changes. As a result, this function consists of a sequence of calls to definition
# generators. New generators will have to be added as the need arises.
#
# The comments before each set of generator calls identifies the change in the SKIRT 9 code
# requiring the corresponding ski file modification(s).
#
def _getUpgradeDefinitions():
    # construct dictionaries used by the upgrade definitions replacing existing probes by form probes (see below)
    propsAtPositionsForm = dict(filename=None, useColumns=None)
    propsLinearCutForm = dict(numSamples=None, startX=None, startY=None, startZ=None, endX=None, endY=None, endZ=None)
    propsMeridionalCutForm = dict(numSamples=None, radius=None, azimuth=None)
    grid = "//MediumSystem/grid/*"
    propsPlanarCutsForm = dict(minX=grid, maxX=grid, minY=grid, maxY=grid, minZ=grid, maxZ=grid,
                               positionX=None, positionY=None, positionZ=None,
                               numPixelsX=None, numPixelsY=None, numPixelsZ=None)
    propsParallelProjectionForm = dict(inclination=None, azimuth=None, roll=None,
                                       fieldOfViewX=None, numPixelsX=None, centerX=None,
                                       fieldOfViewY=None, numPixelsY=None, centerY=None)
    propsAllSkyProjectionForm = dict(numPixelsY=None, observerX=None, observerY=None, observerZ=None,
                                     crossX=None, crossY=None, crossZ=None, upX=None, upY=None, upZ=None)

    # construct the list of upgrade definitions
    return [
        # SKIRT update (15 jan 2020): replace single bulk velocity with velocity field for geometric sources and media
        # !!! We just remove the zero-velocity components without upgrading non-zero velocities
        # !!! A true upgrade is hard to implement and would be applicable to very few ski files
        _removeScalarPropertyWithValue("GeometricSource", "velocityX", "0 "),
        _removeScalarPropertyWithValue("GeometricSource", "velocityY", "0 "),
        _removeScalarPropertyWithValue("GeometricSource", "velocityZ", "0 "),
        _removeScalarPropertyWithValue("GeometricMedium", "velocityX", "0 "),
        _removeScalarPropertyWithValue("GeometricMedium", "velocityY", "0 "),
        _removeScalarPropertyWithValue("GeometricMedium", "velocityZ", "0 "),

        # SKIRT update (4 feb 2021): change base type for SED properties of Lya SED decorators to continuous SED
        _changeCompoundPropertyBaseType("LyaSEDDecorator", "sedOriginal", "SED", "ContSED"),
        _changeCompoundPropertyBaseType("LyaSEDDecorator", "sedLymanAlpha", "SED", "ContSED"),
        _changeCompoundPropertyBaseType("LyaSEDFamilyDecorator", "sedLymanAlpha", "SED", "ContSED"),

        # SKIRT update (oct/nov 2021): reorganize simulation mode and medium options to enable gas emission
        # -- simulation mode and top-level iterate options
        _adjustLyaExtinctionMode(),
        _adjustDustSelfAbsorptionMode(),
        _moveScalarProperty("DynamicStateOptions", "hasDynamicState", "MonteCarloSimulation", "iterateMediumState"),
        # -- radiation field options
        _addMediumSystemOptions("RadiationFieldOptions", "ExtinctionOnlyOptions"),
        _addMediumSystemOptions("RadiationFieldOptions", "DustEmissionOptions"),
        _moveScalarProperty("ExtinctionOnlyOptions", "storeRadiationField", "RadiationFieldOptions"),
        _copyCompoundProperty("ExtinctionOnlyOptions", "radiationFieldWLG", "RadiationFieldOptions"),
        _removeCompoundProperty("ExtinctionOnlyOptions", "radiationFieldWLG"),
        _copyCompoundProperty("DustEmissionOptions", "radiationFieldWLG", "RadiationFieldOptions"),
        _removeCompoundProperty("DustEmissionOptions", "radiationFieldWLG"),
        # -- sampling options
        _addMediumSystemOptions("SamplingOptions", "MediumSystem", "numDensitySamples"),
        _moveScalarProperty("MediumSystem", "numDensitySamples", "SamplingOptions"),
        # -- iteration options
        _addMediumSystemOptions("IterationOptions", "DynamicStateOptions"),
        _addMediumSystemOptions("IterationOptions", "DustSelfAbsorptionOptions"),
        _moveScalarProperty("DynamicStateOptions", "minIterations", "IterationOptions", "minPrimaryIterations"),
        _moveScalarProperty("DynamicStateOptions", "maxIterations", "IterationOptions", "maxPrimaryIterations"),
        _moveScalarProperty("DynamicStateOptions", "iterationPacketsMultiplier",
                            "IterationOptions", "primaryIterationPacketsMultiplier"),
        _moveScalarProperty("DustSelfAbsorptionOptions", "minIterations", "IterationOptions", "minSecondaryIterations"),
        _moveScalarProperty("DustSelfAbsorptionOptions", "maxIterations", "IterationOptions", "maxSecondaryIterations"),
        _moveScalarProperty("DustSelfAbsorptionOptions", "iterationPacketsMultiplier",
                            "IterationOptions", "secondaryIterationPacketsMultiplier"),
        # -- secondary emission options
        _addMediumSystemOptions("SecondaryEmissionOptions", "DustEmissionOptions"),
        _moveScalarProperty("DustEmissionOptions", "storeEmissionRadiationField", "SecondaryEmissionOptions"),
        _moveScalarProperty("DustEmissionOptions", "secondaryPacketsMultiplier", "SecondaryEmissionOptions"),
        _moveScalarProperty("DustEmissionOptions", "spatialBias", "SecondaryEmissionOptions"),
        # -- dust emission options
        _moveScalarProperty("DustSelfAbsorptionOptions", "maxFractionOfPrimary", "DustEmissionOptions"),
        _moveScalarProperty("DustSelfAbsorptionOptions", "maxFractionOfPrevious", "DustEmissionOptions"),
        # obsolete option sets
        _removeCompoundProperty("MediumSystem", "extinctionOnlyOptions"),
        _removeCompoundProperty("MediumSystem", "dustSelfAbsorptionOptions"),

        # SKIRT update (may 2022): introduce form probes, replacing many existing probes
        _changeTypeName("SpatialGridConvergenceProbe", "ConvergenceInfoProbe"),
        _changeTypeName("DefaultMediaDensityCutsProbe", "ConvergenceCutsProbe"),
        # PerCellForm
        _changeToFormProbe("DustTemperaturePerCellProbe", "TemperatureProbe", "PerCellForm",
                            dict(probeName=None, probeAfter="Run"), dict()),
        _changeToFormProbe("ElectronTemperaturePerCellProbe", "TemperatureProbe", "PerCellForm",
                            dict(probeName=None, probeAfter="Setup"), dict()),
        _changeToFormProbe("GasTemperaturePerCellProbe", "TemperatureProbe", "PerCellForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("MetallicityPerCellProbe", "MetallicityProbe", "PerCellForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("MediumVelocityPerCellProbe", "VelocityProbe", "PerCellForm",
                            dict(probeName=None), dict()),
        _changeToFormProbe("MagneticFieldPerCellProbe", "MagneticFieldProbe", "PerCellForm",
                            dict(probeName=None), dict()),
        _changeToFormProbe("CustomStatePerCellProbe", "CustomStateProbe", "PerCellForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("RadiationFieldPerCellProbe", "RadiationFieldProbe", "PerCellForm",
                            dict(probeName=None, writeWavelengthGrid=None), dict()),
        # DefaultCutsForm
        _changeToFormProbe("DefaultDustTemperatureCutsProbe", "TemperatureProbe", "DefaultCutsForm",
                            dict(probeName=None, probeAfter="Run"), dict()),
        _changeToFormProbe("DefaultElectronTemperatureCutsProbe", "TemperatureProbe", "DefaultCutsForm",
                            dict(probeName=None, probeAfter="Setup"), dict()),
        _changeToFormProbe("DefaultGasTemperatureCutsProbe", "TemperatureProbe", "DefaultCutsForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("DefaultMetallicityCutsProbe", "MetallicityProbe", "DefaultCutsForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("DefaultMediumVelocityCutsProbe", "VelocityProbe", "DefaultCutsForm",
                            dict(probeName=None), dict()),
        _changeToFormProbe("DefaultMagneticFieldCutsProbe", "MagneticFieldProbe", "DefaultCutsForm",
                            dict(probeName=None), dict()),
        _changeToFormProbe("DefaultCustomStateCutsProbe", "CustomStateProbe", "DefaultCutsForm",
                            dict(probeName=None, probeAfter=None), dict()),
        _changeToFormProbe("DefaultRadiationFieldCutsProbe", "RadiationFieldProbe", "DefaultCutsForm",
                            dict(probeName=None, writeWavelengthGrid=None), dict()),
        # AtPositionsForm
        _changeToFormProbe("RadiationFieldAtPositionsProbe", "RadiationFieldProbe", "AtPositionsForm",
                            dict(probeName=None, writeWavelengthGrid=None), propsAtPositionsForm),
        # LinearCutForm
        _changeToFormProbe("LinearDustTemperatureCutProbe", "TemperatureProbe", "LinearCutForm",
                            dict(probeName=None, probeAfter="Run"), propsLinearCutForm),
        # MeridionalCutForm
        _changeToFormProbe("MeridionalDustTemperatureCutProbe", "TemperatureProbe", "MeridionalCutForm",
                            dict(probeName=None, probeAfter="Run"), propsMeridionalCutForm),
        # PlanarCutsForm
        _changeToFormProbe("PlanarMediaDensityCutsProbe", "DensityProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarDustTemperatureCutsProbe", "TemperatureProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter="Run"), propsPlanarCutsForm),
        _changeToFormProbe("PlanarElectronTemperatureCutsProbe", "TemperatureProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter="Setup"), propsPlanarCutsForm),
        _changeToFormProbe("PlanarGasTemperatureCutsProbe", "TemperatureProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarMetallicityCutsProbe", "MetallicityProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarMediumVelocityCutsProbe", "VelocityProbe", "PlanarCutsForm",
                            dict(probeName=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarMagneticFieldCutsProbe", "MagneticFieldProbe", "PlanarCutsForm",
                            dict(probeName=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarCustomStateCutsProbe", "CustomStateProbe", "PlanarCutsForm",
                            dict(probeName=None, probeAfter=None), propsPlanarCutsForm),
        _changeToFormProbe("PlanarRadiationFieldCutsProbe", "RadiationFieldProbe", "PlanarCutsForm",
                            dict(probeName=None, writeWavelengthGrid=None), propsPlanarCutsForm),
        # ParallelProjectionForm
        _changeToFormProbe("ProjectedMediaDensityProbe", "DensityProbe", "ParallelProjectionForm",
                            dict(probeName=None, probeAfter=None), propsParallelProjectionForm),
        # AllSkyProjectionForm
        _changeToFormProbe("OpticalDepthMapProbe", "OpacityProbe", "AllSkyProjectionForm",
                            dict(probeName=None, wavelength=None, probeAfter=None), propsAllSkyProjectionForm),

        # SKIRT update (aug 2022): revise dynamic medium state concepts; rename iterateMediumState property
        _changeScalarPropertyName("MonteCarloSimulation", "iterateMediumState", "iteratePrimaryEmission"),

        # SKIRT update (sep 2022): allow specifying characteristic wavelengths for file/list border wavelength grid
        _changeBooleanToEnumeration("FileBorderWavelengthGrid", "log", "characteristic", "Linear", "Logarithmic"),
        _changeBooleanToEnumeration("ListBorderWavelengthGrid", "log", "characteristic", "Linear", "Logarithmic"),

        # SKIRT update (feb 2024): all Mesh classes are now movable; replace type="MoveableMesh" by type="Mesh"
        _changeScalarPropertyValue("meshX", "type", "MoveableMesh", "Mesh"),
        _changeScalarPropertyValue("meshY", "type", "MoveableMesh", "Mesh"),
        _changeScalarPropertyValue("meshZ", "type", "MoveableMesh", "Mesh"),

        # SKIRT update (sep 2025): OpacityProbe expects compound wavelengthGrid instead of scalar wavelength
        _replaceWavelengthByGrid("OpacityProbe"),

        # SKIRT 10 update (oct 2026): remove FullInstrument; replace it by an SEDInstrument and a FrameInstrument
        _replaceFullInstrument(),
    ]

# --------- handling probe to form-probe updates

# Replace the given old probe by the given new probe with the given nested form;
# probeProps and formProps are dictionaries listing the scalar properties for the new probe and form;
# the value for each property in the dictionaries can be:
#  - a string that does not start with "/": this literal value is used,
#  - a string that starts with "/": the property is copied from the element at the specified path,
#  - None: the value is copied from the old probe;
# any compound properties of the old probe are copied into the form.
def _changeToFormProbe(oldProbeName, probeName, formName, probeProps, formProps):
    # create a list of templates for the probe and form properties
    def constructPropsString(props):
        propsString = ""
        for name, value in props.items():
            if value is None:
                value = "@{0}".format(name)
            elif value.startswith("/"):
                value = value+"/@{0}".format(name)
            else:
                value= "'{0}'".format(value)
            propsString += '''<xsl:attribute name="{0}"> <xsl:value-of select="{1}"/> </xsl:attribute>
                           '''.format(name, value)
        return propsString
    probePropsString = constructPropsString(probeProps)
    formPropsString = constructPropsString(formProps)
    # create the full template
    return ('''//{0}'''.format(oldProbeName),
            '''
            <xsl:template match="//{0}">
                <xsl:element name="{1}">
                    {3}
                    <xsl:element name="form">
                        <xsl:attribute name="type">
                            <xsl:value-of select="'Form'"/>
                        </xsl:attribute>
                        <xsl:element name="{2}">
                            {4}
                            <xsl:apply-templates select="node()"/>
                        </xsl:element>
                    </xsl:element>
                </xsl:element>
            </xsl:template>
            '''.format(oldProbeName, probeName, formName, probePropsString, formPropsString))

# --------- handling specific (non-generic) updates

# Replace simulation mode "LyaWithDustExtinction" with "LyaExtinctionOnly"
def _adjustLyaExtinctionMode():
    return ('''//MonteCarloSimulation[@simulationMode='LyaWithDustExtinction']''',
            '''
            <xsl:template match="//MonteCarloSimulation/@simulationMode">
                <xsl:attribute name="simulationMode">
                    <xsl:value-of select="'LyaExtinctionOnly'"/>
                </xsl:attribute>
            </xsl:template>
            ''')

# Replace simulation mode "DustEmissionWithSelfAbsorption" with "DustEmission"
# and add iterateSecondaryEmission="true"
def _adjustDustSelfAbsorptionMode():
    return ('''//MonteCarloSimulation[@simulationMode='DustEmissionWithSelfAbsorption']''',
            '''
            <xsl:template match="//MonteCarloSimulation">
                <xsl:element name="MonteCarloSimulation">
                    <xsl:apply-templates select="@*"/>
                    <xsl:attribute name="simulationMode">
                        <xsl:value-of select="'DustEmission'"/>
                    </xsl:attribute>
                    <xsl:attribute name="iterateSecondaryEmission">
                        <xsl:value-of select="'true'"/>
                    </xsl:attribute>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            ''')

# Add a new MediumSystem options section of the given type if a given type/property exists
def _addMediumSystemOptions(optionsTypeName, oldTypeName, oldPropName=None):
    optionsPropName = optionsTypeName[0].lower() + optionsTypeName[1:]
    condition = "{0}" if oldPropName is None else "{0}/@{1}"
    condition = "boolean(//{0}) and not(//{1})".format(condition, optionsTypeName)
    return (condition.format(oldTypeName, oldPropName),
            '''
            <xsl:template match="//MediumSystem">
                <xsl:element name="MediumSystem">
                    <xsl:apply-templates select="@*"/>
                        <{0} type="{1}">
                            <{1}/>
                        </{0}>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(optionsPropName, optionsTypeName))

# replace scalar wavelength property by compound wavelengthGrid property
def _replaceWavelengthByGrid(typeName):
    return ('''//{0}[@wavelength]'''.format(typeName),
            '''
            <xsl:template match="//{0}[@wavelength]">
                <xsl:element name="{0}">
                    <xsl:apply-templates select="@*[(local-name() != 'wavelength')]"/>
                    <wavelengthGrid type="WavelengthGrid">
                        <ListWavelengthGrid>
                            <xsl:attribute name="wavelengths">
                                <xsl:value-of select="@wavelength"/>
                            </xsl:attribute>
                        </ListWavelengthGrid>
                    </wavelengthGrid>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(typeName))

# Replace each FullInstrument by an SEDInstrument followed by a FrameInstrument with the same name and line of sight;
# placed consecutively, the two instruments receive the same peel-off photon packets. The SEDInstrument receives all
# properties except for the frame properties and has no aperture, like the SED recorded by the FullInstrument. Both
# instruments receive a copy of the instrument-specific wavelength grid, if any.
def _replaceFullInstrument():
    frameProps = ["fieldOfViewX", "numPixelsX", "centerX", "fieldOfViewY", "numPixelsY", "centerY"]
    notFrameProp = " and ".join("local-name() != '{0}'".format(prop) for prop in frameProps)
    return ('''//FullInstrument''',
            '''
            <xsl:template match="//FullInstrument">
                <xsl:element name="SEDInstrument">
                    <xsl:apply-templates select="@*[{0}]"/>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
                <xsl:element name="FrameInstrument">
                    <xsl:apply-templates select="@*|node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(notFrameProp))

# --------- handling types

## Rename a given type.
def _changeTypeName(oldTypeName, newTypeName):
    return ('''//{0}'''.format(oldTypeName),
            '''
            <xsl:template match="//{0}">
                <xsl:element name="{1}">
                    <xsl:apply-templates select="@*|node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(oldTypeName, newTypeName))

# --------- handling scalar properties

## Add a scalar property to a given type.
def _addScalarProperty(typeName, propName, propStringValue):
    return ('''//{0}[not(@{1})]'''.format(typeName, propName),
            '''
            <xsl:template match="//{0}[not(@{1})]">
                <xsl:element name="{0}">
                    <xsl:apply-templates select="@*"/>
                    <xsl:attribute name="{1}">
                        <xsl:value-of select="{2}"/>
                    </xsl:attribute>
                    <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(typeName, propName, propStringValue))

## Remove a scalar property from a given type.
def _removeScalarProperty(typeName, oldPropName):
    return ('''//{0}/@{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/@{1}">
            </xsl:template>
            '''.format(typeName, oldPropName))

## Remove a scalar property with a value starting with a given string from a given type.
def _removeScalarPropertyWithValue(typeName, oldPropName, oldValue):
    return ('''//{0}[starts-with(@{1},'{2}')]'''.format(typeName, oldPropName, oldValue),
            '''
            <xsl:template match="//{0}[starts-with(@{1},'{2}')]/@{1}">
            </xsl:template>
            '''.format(typeName, oldPropName, oldValue))

## Change the name of a scalar property for a given type.
def _changeScalarPropertyName(typeName, oldPropName, newPropName):
    return ('''//{0}/@{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/@{1}">
                <xsl:attribute name="{2}">
                    <xsl:value-of select="."/>
                </xsl:attribute>
            </xsl:template>
            '''.format(typeName, oldPropName, newPropName))

## Change the value of a given scalar property for a given type.
def _changeScalarPropertyValue(typeName, propName, oldValue, newValue):
    return ('''//{0}[@{1}='{2}']'''.format(typeName, propName, oldValue),
            '''
            <xsl:template match="//{0}[@{1}='{2}']/@{1}">
                <xsl:attribute name="{1}">
                    <xsl:value-of select="'{3}'"/>
                </xsl:attribute>
            </xsl:template>
            '''.format(typeName, propName, oldValue, newValue))

## Change Boolean property to enumeration property.
def _changeBooleanToEnumeration(typeName, oldPropName, newPropName, newValueForFalse, newValueForTrue):
    return ('''//{0}/@{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/@{1}">
                <xsl:attribute name="{2}">
                    <xsl:choose>
                        <xsl:when test=" .='true' or .='True' or .='1' ">
                            <xsl:value-of select="'{4}'"/>
                        </xsl:when>
                        <xsl:otherwise>
                            <xsl:value-of select="'{3}'"/>
                        </xsl:otherwise>
                    </xsl:choose>
                </xsl:attribute>
            </xsl:template>
            '''.format(typeName, oldPropName, newPropName, newValueForFalse, newValueForTrue))

## Move a scalar property and its value from one type to another, optionally using
# a new property name on the target type (the target type is assumed to exist).
def _moveScalarProperty(oldTypeName, oldPropName, newTypeName, newPropName=None):
   return ('''//{0}/@{1}'''.format(oldTypeName, oldPropName),
           '''
           <xsl:template match="//{2}">
               <xsl:element name="{2}">
                   <xsl:apply-templates select="@*"/>
                   <xsl:attribute name="{3}">
                       <xsl:value-of select="//{0}/@{1}"/>
                   </xsl:attribute>
                   <xsl:apply-templates select="node()"/>
               </xsl:element>
           </xsl:template>
           <xsl:template match="//{0}/@{1}">
           </xsl:template>
           '''.format(oldTypeName, oldPropName, newTypeName,
                      oldPropName if newPropName is None else newPropName))

# --------- handling compound properties

## Change the name of a compound property for a given type.
def _changeCompoundPropertyName(typeName, oldPropName, newPropName):
    return ('''//{0}/{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/{1}">
                <xsl:element name="{2}">
                    <xsl:apply-templates select="@*|node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(typeName, oldPropName, newPropName))

## Change the base type of a compound property for a given type.
def _changeCompoundPropertyBaseType(typeName, propName, oldBaseType, newBaseType):
    return ('''//{0}/{1}[@type='{2}']'''.format(typeName, propName, oldBaseType),
            '''
            <xsl:template match="//{0}/{1}/@type">
                <xsl:attribute name="type">
                    <xsl:value-of select="'{2}'"/>
                </xsl:attribute>
            </xsl:template>
            '''.format(typeName, propName, newBaseType))

## Remove a compound property from a given type.
def _removeCompoundProperty(typeName, oldPropName):
    return ('''//{0}/{1}'''.format(typeName, oldPropName),
            '''
            <xsl:template match="//{0}/{1}">
            </xsl:template>
            '''.format(typeName, oldPropName))

## Copy a compound property and its value from one type to another
# (the target type is assumed to exist)
def _copyCompoundProperty(oldTypeName, oldPropName, newTypeName):
    return ('''//{0}/{1}'''.format(oldTypeName, oldPropName),
            '''
            <xsl:template match="//{2}">
                <xsl:element name="{2}">
                     <xsl:apply-templates select="@*"/>
                     <xsl:apply-templates select="//{0}/{1}"/>
                     <xsl:apply-templates select="node()"/>
                </xsl:element>
            </xsl:template>
            '''.format(oldTypeName, oldPropName, newTypeName))

# -----------------------------------------------------------------
