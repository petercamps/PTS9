#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.skirt Executing the SKIRT command line application
#
# An instance of the Skirt class in this module represents a particular SKIRT executable,
# and allows invoking it with given command line arguments.
#

# -----------------------------------------------------------------

import subprocess
import sys
import pts.utils as ut
from .simulation import Simulation

# -----------------------------------------------------------------

## An instance of the Skirt class represents a particular SKIRT executable, and allows invoking it with
# given command line arguments for execution on the local host.
class Skirt:

    ## The constructor accepts an optional argument specifying the file path to the SKIRT executable to be used.
    # If the path is specified, it should include the "skirt" filename of the executable. Otherwise the path
    # is interpreted as described for the pts.utils.absPath() function.
    # If the path is not specified, the constructor looks for a SKIRT executable in the project structure including
    # this PTS source file, assuming that SKIRT is built under a \c ~/SKIRT9 or \c SKIRT directory residing next
    # to the PTS9 or PTS directory.
    def __init__(self, path=None):

        # set the SKIRT path
        if path is None:
            self._path = ut.skirtPath()
            if self._path is None:
                raise ValueError("Cannot locate default SKIRT executable in PTS/SKIRT project directory structure")
        else:
            self._path = ut.absPath(path)
            if not self._path.is_file():
                raise ValueError("Specified SKIRT executable does not exist: {}".format(self._path))

        # initialize execution state
        self._process = None

    ## This function returns the absolute file path of the SKIRT executable used by this instance.
    def path(self):
        return self._path

    ## This function returns a string with version information on the SKIRT executable represented by this
    # object. The function invokes SKIRT with an incorrect command line argument to obtain this information.
    def version(self):
        # execute skirt with incorrect argument list and get its output
        process = subprocess.run([str(self._path), "--version"], capture_output=True, text=True)
        output = process.stdout

        # return the relevant portion of the output
        return "SKIRT" + output.splitlines()[0].partition("SKIRT")[2]

    ## This function invokes the SKIRT executable with the simulation and command line options corresponding to the
    #  values of the function arguments, as described below. The function returns a pts.simulation.simulation.Simulation
    #  instance corresponding to the simulation being performed (or attempted, in case of failure). The function
    #  supports asynchronous execution, allowing the caller to perform other tasks while the simulation is running.
    #
    # - \em skiFilePath: the file path for the configuration file (\em ski file) to be executed, interpreted as
    #   as described for the pts.utils.absPath() function.
    # - \em inDirPath: a string specifying the absolute or relative directory path for simulation input files,
    #   optionally followed by the name of an HDF5 file in that directory and the path of a suite in that file, as
    #   for the \c -i command line option of SKIRT.
    # - \em outDirPath: a string specifying the absolute or relative directory path for simulation output files,
    #   optionally followed by the name of an HDF5 file in that directory and the path of a suite in that file, as
    #   for the \c -o command line option of SKIRT. The returned Simulation instance refers to the directory, in
    #   which SKIRT always writes the log file.
    # - \em skiRelative: if \c True, the simulation input/output directory paths are relative to the directory of the
    #   \em ski file being executed; if \c False or missing, they are relative to the current directory.
    #
    # - \em numThreadsPerProcess: a positive integer specifying the number of parallel threads per process for the
    #   simulation. If zero or missing the number of logical cores on the host computer is used.
    # - \em numProcesses: a positive integer specifying the number of parallel MPI processes to be launched,
    #   or a string indicating that the number of processes and other aspects of the MPI configuration should be
    #   obtained at run time from the corresponding multi-node batch queueing system. Supported string values include
    #   'lsf' for the LSF queueing system and 'srun' for the SLURM queueing system. If the option is missing or
    #   has a value of zero or one, MPI is not used and just a single SKIRT process is launched.
    # - \em verbose: this option has effect only if the number of processes is larger than one. If set to \c True, each
    #   process creates its own complete log file. If missing or set to \em False, only the root process creates a
    #   full log file, and the other processes only create a log file when there are errors or warnings.
    #
    # - \em wait: if \c True or missing, the function waits until SKIRT execution completes; if \c False the function
    #   returns immediately after launching SKIRT without waiting for it.
    # - \em console: a string controlling console logging in synchronous execution mode: 'regular' means all SKIRT
    #   messages are logged to the console; 'brief' means only summary messages are logged to the console, and
    #   'silent' means no SKIRT messages are logged to the console. If \em wait is \c True or missing, the default
    #   logging mode is 'regular'. If \em wait is \c False, the logging mode is always 'silent'.
    #
    def execute(self, skiFilePath, *, inDirPath="", outDirPath="", skiRelative=False,
                numThreadsPerProcess=0, numProcesses=1, verbose=False, wait=True, console='regular'):

        if self.isRunning():
            raise ValueError("Calling execute on SKIRT object that is still executing")

        # --- build the argument list ---

        # mpi support
        if isinstance(numProcesses, str):
            if numProcesses == 'lsf':
                arguments = ["mpirun", "-lsf"]
            elif numProcesses == 'srun':
                arguments = ["mpirun", "-srun"]
            else:
                raise ValueError("Unsupported string value for numProcesses: {}", numProcesses)
        else:
            numProcesses = int(numProcesses)
            if numProcesses > 1:
                arguments = ["mpirun", "-np", str(numProcesses)]
            else:
                arguments = []

        # skirt executable
        arguments += [str(self._path)]

        # ski file
        arguments += [str(ut.absPath(skiFilePath))]

        # i/o path options
        if skiRelative:
            base = ut.absPath(skiFilePath).parent
            inpath = ut.absPath(base / inDirPath)
            outpath = ut.absPath(base / outDirPath)
        else:
            inpath = ut.absPath(inDirPath)
            outpath = ut.absPath(outDirPath)
        arguments += ["-i", str(inpath)]
        arguments += ["-o", str(outpath)]

        # parallelization options
        numThreadsPerProcess = int(numThreadsPerProcess)
        if numThreadsPerProcess > 0:
            arguments += ["-t", str(numThreadsPerProcess)]
        if verbose:
            arguments += ["-v"]

        # logging options
        if console!='regular' or not wait:
            arguments += ["-b"]

        # --- launch SKIRT ---

        if wait:
            self._process = None
            if console=='silent':
                subprocess.run(arguments, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                # we pipe the SKIRT output to the console ourselves because this also works in Jupyter notebooks
                popen = subprocess.Popen(arguments, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         universal_newlines=True)
                for line in iter(popen.stdout.readline, ""): print(line, end="")
                popen.stdout.close()
                popen.wait()
        else:
            self._process = subprocess.Popen(arguments, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        return Simulation(skiFilePath=skiFilePath, inDirPath=_directoryPath(inpath),
                          outDirPath=_directoryPath(outpath), process=self._process)

    ## This function returns True if the simulation started with the most recent call to the execute() function
    # is still running, and False otherwise.
    def isRunning(self):
        if self._process is not None:
            if self._process.poll() is None:
                return True
            else:
                # allow the process object to be freed
                self._process = None
        return False

# -----------------------------------------------------------------

## This helper function returns the directory part of the specified absolute path, given as the value of the \c -i
# or \c -o command line option of SKIRT: the path itself if it is an existing directory or if it does not name an
# HDF5 file, and otherwise the directory that holds the HDF5 file. As in SKIRT, the path names an HDF5 file if it
# includes the \c .hdf5 filename extension followed by a colon, or if it ends with that extension.
def _directoryPath(path):
    if path.is_dir(): return path
    text = str(path)
    index = text.find(".hdf5:")
    if index >= 0: return ut.absPath(text[:index+5]).parent
    if text.endswith(".hdf5"): return path.parent
    return path

# -----------------------------------------------------------------
