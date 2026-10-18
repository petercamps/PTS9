#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.utils.path Special paths related to PTS
#
# This module allows retrieving some special paths related to PTS, such as for example
# the path to the pts repository. It also offers the BundlePath class, which represents
# the equivalent of a regular file stored as a bundle in an HDF5 file.
#

# -----------------------------------------------------------------

import functools
import inspect
import pathlib
from .error import UserError

# -----------------------------------------------------------------

## This function returns the absolute path to the \c data subdirectory for the PTS package
# in which the given object is defined.
def dataPath(object):
    packagedir = pathlib.Path(inspect.getfile(object)).parent
    if packagedir.name == "do": packagedir = packagedir.parent
    return packagedir / "data"

## This function returns the absolute path to the top-level directory of the pts repository
# as a pathlib.Path object.
def ptsPath():
    return pathlib.Path(inspect.getfile(ptsPath)).parent.parent

## This function returns the absolute path to the PTS resources directory as part of the project structure including
# this PTS source file, or None if no such directory is found.
def ptsResourcesPath():
    path = ptsPath().parent / "resources"
    if path.is_dir(): return path
    return None

## This function returns the absolute path to the directory containing the SKIRT and PTS
# project directories as a pathlib.Path object. In the canonical SKIRT/PTS developer
# directory structure, this is the grandparent of the pts repository.
def projectParentPath():
    return ptsPath().parent.parent

## This function returns the absolute path to the SKIRT executable built as part of the project structure including
# this PTS source file, or None if no executable is found. The function assumes that a SKIRT release version is built
# under a \c ~/SKIRT9 or \c SKIRT directory in the directory returned by the projectParent() function. The \c SKIRT9
# directory is tried first.
def skirtPath():
    project = projectParentPath()
    path = project / "SKIRT9" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    path = project / "SKIRT" / "release" / "SKIRT" / "main" / "skirt"
    if path.is_file(): return path
    return None

## This function returns the absolute path to the SKIRT resources directory as part of the project structure including
# this PTS source file, or None if no such directory is found. The function assumes that SKIRT is installed
# under a \c ~/SKIRT9 or \c SKIRT directory in the directory returned by the projectParent() function. The \c SKIRT9
# directory is tried first.
def skirtResourcesPath():
    project = projectParentPath()
    path = project / "SKIRT9" / "resources"
    if path.is_dir(): return path
    path = project / "SKIRT" / "resources"
    if path.is_dir(): return path
    return None

## This function returns the absolute canonical path corresponding to the given path. The input path may
# be specified as a string or as a pathlib.Path object. If the input path is relative, it is interpreted
# relative to the current working directory. If it starts with a tilde, the tilde is expanded to the home
# directory of the current user. In all cases, the path is streamlined by removing and "." and ".." segments
# and following symbolic links.
def absPath(inpath):
    path = pathlib.Path(inpath).expanduser()
    if not path.is_absolute(): path = pathlib.Path.cwd() / path
    return path.resolve()

## This function returns the absolute canonical path corresponding to the file path that will be used to save
# a result to file, such as a PDF plot or a PNG image. The returned path is derived from the input arguments:
#  - \em defFilePath: required argument; specifies the default save file path and is used in case the other arguments
#    do not fully specify the path. This path is typically constructed by the immediate caller, e.g. a plot function.
#    It can also be a BundlePath instance, which then stands for its regular file path, so that a result derived
#    from a bundle in an HDF5 file is saved next to the HDF5 file.
#  - \em outDirPath: optional argument; overrides the directory path in \em defFilePath.
#  - \em outFileName: optional argument; overrides the filename portion of the path in \em defFilePath. It is valid
#    for both \em outDirPath and \em outFileName to be specified, in which case \em defFilePath is ignored. However,
#    in this case it may be more natural to use the \em outFilePath argument instead.
#  - \em outFilePath: optional argument; overrides the complete path, which means that \em defFilePath, \em outDirPath
#    and \em outFileName are ignored. The three optional arguments are typically (but not necessarily) supplied by
#    the caller of the immediate caller.
#  - \em suffix: required argument; string or sequence of strings listing the allowed suffixes for the returned path
#    (including the leading "." in each suffix).
#    If the path resulting from the above rules does not already has one of the specified suffixes, the suffix
#    of the path will be replaced by the first (or only) specified suffix.
#
# The three input paths may be specified as a string or as a pathlib.Path object. If the input path is relative,
# it is interpreted relative to the current working directory. If it starts with a tilde, the tilde is expanded
# to the home directory of the current user.
def savePath(defFilePath, suffix, *, outDirPath=None, outFileName=None, outFilePath=None):

    # a bundle stands for its regular file path
    if isinstance(defFilePath, BundlePath): defFilePath = defFilePath.regularFilePath()

    # construct the path based on defFilePath, outDirPath, outFileName, and outFilePath
    result = None
    if outFilePath is not None: result = absPath(outFilePath)
    elif outDirPath is not None and outFileName is not None: result = absPath(outDirPath) / outFileName
    elif outFileName is not None: result = absPath(defFilePath).with_name(outFileName)
    elif outDirPath is not None: result = absPath(outDirPath) / absPath(defFilePath).name
    else: result = absPath(defFilePath)

    # adjust the suffix as requested
    if isinstance(suffix, str): suffix = [ suffix ]
    if not (result.suffix in suffix):
        result = result.with_suffix(suffix[0])

    return result

# -----------------------------------------------------------------

## An instance of the BundlePath class represents a bundle in an HDF5 file, i.e. the equivalent of a regular file
# stored in an HDF5 file (see the pts.hdf5 package). For example, the functions in the pts.simulation package that
# list the output files of a SKIRT simulation return a BundlePath instance for each output file that SKIRT wrote as
# a bundle, and the functions in that package that load the contents of an output file accept such an instance.
#
# A bundle is identified by its <em>bundle path</em>, a string of the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt>,
# where \c \<filepath\> is the path of the HDF5 file, including the mandatory \c .hdf5 filename extension;
# \c \<suite\> is the optional path of the suite within the file; and \c \<bundle\> is the name of the bundle, which
# usually is the name of the regular file it replaces. A BundlePath instance also has a <em>regular file path</em>:
# the path of a file with the name of the bundle in the directory that holds the HDF5 file. This is where the
# regular file would reside if it had not been stored in the HDF5 file: SKIRT writes its regular output files next
# to its HDF5 output file, and looks for a regular input file next to its HDF5 input file before reading the bundle.
#
# So that code written for regular files can handle bundles as well, a BundlePath instance offers the properties and
# functions of a pathlib.Path object that do not access the file system, applied to its regular file path: the
# properties \em name, \em stem, \em suffix, and \em parent, and the functions with_name(), with_stem(), and
# with_suffix(), which return a pathlib.Path object. A result derived from a bundle, such as a plot, thus gets a
# regular file path next to the HDF5 file; the savePath() function accepts a BundlePath instance for the same reason.
# On the other hand, a BundlePath instance cannot be used to open a file: a function that reads a file must
# recognize a bundle (see the isBundlePath() function) and read it through the pts.hdf5 package.
#
# Converting a BundlePath instance to a string yields its bundle path, with an absolute file path. BundlePath
# instances are equal if their bundle paths are equal. They are ordered on their regular file path, also with
# respect to pathlib.Path objects, so that a list of output files can be sorted whether or not they are bundles.
#
@functools.total_ordering
class BundlePath:

    ## The constructor accepts a bundle path as described above, as a string or as an object that converts to such
    # a string, such as a BundlePath instance. The file path is interpreted as described for the absPath() function.
    # The constructor raises a user error if the argument does not refer to a bundle in an HDF5 file.
    def __init__(self, bundlePath):
        text = str(bundlePath)
        index = text.find(_hdf5Extension + ":")
        groupPath = ""
        if index >= 0:
            index += len(_hdf5Extension)
            groupPath = "/".join([ name for name in text[index+1:].split("/") if name ])
        if not groupPath:
            raise UserError("Path does not refer to a bundle in an HDF5 file: '{}'".format(text))
        self._hdf5FilePath = absPath(text[:index])
        self._groupPath = groupPath
        self._regularFilePath = self._hdf5FilePath.parent / groupPath.rpartition("/")[2]

    ## This function returns the absolute path of the HDF5 file that holds the bundle, as a pathlib.Path object.
    def hdf5FilePath(self):
        return self._hdf5FilePath

    ## This function returns the path of the bundle within the HDF5 file, i.e. the path of the suite, if any,
    # followed by a forward slash and the name of the bundle.
    def groupPath(self):
        return self._groupPath

    ## This function returns the regular file path of the bundle, as a pathlib.Path object.
    def regularFilePath(self):
        return self._regularFilePath

    ## This property returns the name of the bundle, which is the filename of its regular file path.
    @property
    def name(self):
        return self._regularFilePath.name

    ## This property returns the stem of the filename of the regular file path.
    @property
    def stem(self):
        return self._regularFilePath.stem

    ## This property returns the filename extension of the regular file path.
    @property
    def suffix(self):
        return self._regularFilePath.suffix

    ## This property returns the directory that holds the HDF5 file, which is the parent of the regular file path.
    @property
    def parent(self):
        return self._regularFilePath.parent

    ## This function returns the regular file path with the filename replaced by the specified name.
    def with_name(self, name):
        return self._regularFilePath.with_name(name)

    ## This function returns the regular file path with the stem of the filename replaced by the specified stem.
    def with_stem(self, stem):
        return self._regularFilePath.with_stem(stem)

    ## This function returns the regular file path with the filename extension replaced by the specified suffix.
    def with_suffix(self, suffix):
        return self._regularFilePath.with_suffix(suffix)

    ## This function returns the bundle path, with an absolute file path.
    def __str__(self):
        return "{}:{}".format(self._hdf5FilePath, self._groupPath)

    ## This function returns a representation of the bundle path for use in debugging.
    def __repr__(self):
        return "BundlePath('{}')".format(self)

    ## This function returns True if the other object is a BundlePath instance with the same bundle path.
    def __eq__(self, other):
        return isinstance(other, BundlePath) and str(self) == str(other)

    ## This function returns a hash value consistent with the equality operator.
    def __hash__(self):
        return hash(str(self))

    ## This function returns True if this bundle precedes the other object, which can be a BundlePath instance or a
    # pathlib.Path object, in the order described for this class.
    def __lt__(self, other):
        if not isinstance(other, (BundlePath, pathlib.PurePath)): return NotImplemented
        return _orderingKey(self) < _orderingKey(other)

## This function returns True if the specified path refers to a bundle in an HDF5 file, and False otherwise. The path
# refers to a bundle if it is a BundlePath instance, or a string or pathlib.Path object that includes the \c .hdf5
# filename extension followed by a colon, which is the rule used by SKIRT.
def isBundlePath(path):
    return isinstance(path, BundlePath) or (_hdf5Extension + ":") in str(path)

## This helper function returns the key on which BundlePath instances and pathlib.Path objects are ordered.
def _orderingKey(path):
    if isinstance(path, BundlePath): return (path.regularFilePath(), str(path))
    return (pathlib.Path(path), "")

# the filename extension of an HDF5 file
_hdf5Extension = ".hdf5"

# -----------------------------------------------------------------
