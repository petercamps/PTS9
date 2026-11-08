#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.simulation.healpix Handling of HEALPix grids
#
# This module offers functionality to work with HEALPix grids output by the HEALPixSkyInstrument, and by an
# EscapeInstrument with a HEALPixEscapeBinning and the DataCubes output style, which uses the same layout of the
# HEALPix pixels in its data cubes. For the latter, the angles theta and phi are the inclination and azimuth of the
# escape direction, and the center of a projection map without rotation lies at inclination 90° and azimuth 180°.
#

# -----------------------------------------------------------------

import numpy as np
from pts.utils.error import UserError
import logging
from scipy.spatial import cKDTree

# -----------------------------------------------------------------

## Get the shortest distance (in radians) between two points (A, B) on the celestial sphere with the
# given zenith and azimuth.
#
# The distance can easily be derived by realizing that the angle between the unit vectors pointing to
# A and B is given by
# \f[
#    \delta{} = \arccos{\vec{r}_A \dot{} \vec{r}_B}.
# \f]
# In spherical coordinates, the dot product of the two unit vectors with sky coordinates \f$(\theta{}_A,\phi{}_A)\f$
# and \f$(\theta{}_B,\phi{}_B)\f$ yields
# \f[
#    \vec{r}_A \dot{} \vec{r}_B =
#      \sin{\theta{}_A}\sin{\theta{}_B} \left(
#        \cos{\phi{}_A}\cos{\phi{}_B} + \sin{\phi{}_A}\sin{\phi{}_B}
#      \right) + \cos{\theta{}_A}\cos{\theta{}_B} =
#      \sin{\theta{}_A}\sin{\theta{}_B} \cos{(\phi{}_A-\phi{}_B)}
#        + \cos{\theta{}_A}\cos{\theta{}_B} =
#      \frac{1}{2} \left[
#        \cos{(\theta{}_A-\theta{}_B)} - \cos{(\theta{}_A+\theta{}_B)}
#      \right] \cos{(\phi{}_A-\phi{}_B)} + \frac{1}{2} \left[
#        \cos{(\theta{}_A-\theta{}_B)} + \cos{(\theta{}_A+\theta{}_B)}
#      \right].
# \f]
#
def getDistanceOnSky(thetaA, phiA, thetaB, phiB):
    cosThetaDiff = np.cos(thetaA - thetaB)
    cosThetaSum = np.cos(thetaA + thetaB)
    cosPhiDiff = np.cos(phiA - phiB)
    return np.arccos(
        0.5 * (cosThetaDiff - cosThetaSum) * cosPhiDiff
        + 0.5 * (cosThetaDiff + cosThetaSum)
    )


# -----------------------------------------------------------------

## This class represents a HEALPix grid. It contains a reference to the raw
# HEALPixSkyInstrument data cube and some additional geometrical information
# about the grid that is derived from it.
#
# It offer functionality to access, manipulate and project the HEALPix pixels.
#
class HEALPixGrid:

    ## Constructor.
    #
    # Responsible for initialising the geometrical information about the HEALPix grid.
    #
    def __init__(self, HEALPixCube):
        self._HEALPixCube = HEALPixCube
        self._NSide = HEALPixCube.shape[-1] // 4
        self._order = int(np.log2(self._NSide))
        self._j = None
        self._i = None
        self._theta = None
        self._phi = None
        self._positions = None
        self._tree = None

    ## Print some information about the HEALPix grid to the console.
    def printInfo(self):
        logging.info(
            "HEALPix grid has side length {0} and order {1}.".format(
                self._NSide, self._order
            )
        )
        if len(self._HEALPixCube) == 3:
            logging.info(
                "It contains {0} variables per pixel.".format(
                    self._HEALPixCube.shape[0]
                )
            )

    ## This function returns the angle arrays corresponding to the given input index arrays.
    #
    # The index array j contains the ring indices, the index array i contains the pixel-in-ring
    # indices. The return arrays contain the centers of the corresponding pixels.
    #
    # This function is the inverse of getHEALPixIndices(). No checks are done on the input values
    # to ensure that they are actually valid indices.
    #
    def getHEALPixAngles(self, j, i):

        # make sure the input values are NumPy arrays
        if not j is np.ndarray:
            j = np.array(j)
        if not i is np.ndarray:
            i = np.array(i)

        # reserve empty result arrays
        theta = np.zeros(j.shape)
        phi = np.zeros(i.shape)

        # cache the total number of pixels that is used below
        Ntot = 12 * self._NSide ** 2

        # now compute the angles
        # we treat the poles and equatorial region separately

        # north pole
        iNPol = j < self._NSide - 1
        z = 1.0 - 4.0 * (j[iNPol] + 1.0) ** 2 / Ntot
        theta[iNPol] = np.arccos(z)
        phi[iNPol] = 0.5 * np.pi * (i[iNPol] + 0.5) / (j[iNPol] + 1.0)

        # equator
        iEq = (j >= self._NSide - 1) & (j < 3 * self._NSide)
        z = 2.0 * (2 * self._NSide - j[iEq] - 1.0) / (3 * self._NSide)
        theta[iEq] = np.arccos(z)
        phi[iEq] = (
            (i[iEq] + 0.5 * np.bitwise_and(j[iEq], 1)) * 1.5 * np.pi / (3 * self._NSide)
        )

        # south pole
        iSPol = j >= 3 * self._NSide
        z = 4.0 * (4 * self._NSide - j[iSPol] - 1.0) ** 2 / Ntot - 1.0
        theta[iSPol] = np.arccos(z)
        phi[iSPol] = 0.5 * np.pi * (i[iSPol] + 0.5) / (4 * self._NSide - j[iSPol] - 1.0)

        return theta, phi

    ## This function returns the index arrays corresponding to the given input zenith and azimuth angles for the
    # HEALPix grid in SKIRT's modified ring order.
    #
    # For each pair of input angles, the returned pair of output indices points to the pixel that contains
    # those angles. Using those indices to access the raw HEALPix data will return the corresponding pixel
    # value.
    #
    # Note that the code here is (and always should be) identical to the implementation used internally by
    # SKIRT to guarantee consistency. We do however make use of NumPy arrays and NumPy array indexing for
    # increased efficiency.
    #
    def getHEALPixIndices(self, theta, phi):

        # make sure the input angles are NumPy arrays
        if not theta is np.ndarray:
            theta = np.array(theta)
        if not phi is np.ndarray:
            phi = np.array(phi)

        z = np.cos(theta)
        za = np.abs(z)
        tt = np.mod(2.0 * phi / np.pi, 4.0)

        # empty ring (j) and pixel-in-ring (i) index arrays
        j = np.zeros(za.shape, dtype=int)
        i = np.zeros(za.shape, dtype=int)

        # we need different expressions for HEALPix pixels in the equatorial plane
        # (za <= 2./3.) and in the polar caps (za > 2./3.)
        # to exploit NumPy efficiency, we use index arrays

        # equatorial plane
        ieq = za <= 2.0 / 3.0
        t1 = self._NSide * (0.5 + tt[ieq])
        t2 = 0.75 * self._NSide * z[ieq]
        jp = np.floor(t1 - t2)
        jm = np.floor(t1 + t2)
        j[ieq] = self._NSide + 1 + jp - jm
        kshift = 1 - np.bitwise_and(j[ieq], 1)
        t1 = jp + jm + kshift + 1 + 7 * self._NSide
        i[ieq] = np.bitwise_and(np.right_shift(t1.astype(int), 1), 4 * self._NSide - 1)
        j[ieq] += self._NSide - 2

        # polar caps
        ipol = ~ieq
        ntt = np.floor(tt[ipol])
        tp = tt[ipol] - ntt
        tmp = np.zeros(tp.shape)
        # deal with small values of sin(theta)
        tmp[za[ipol] >= 0.99] = (
            self._NSide
            * np.sin(theta[ipol & (za >= 0.99)])
            / np.sqrt((1.0 + za[ipol & (za >= 0.99)]) / 3.0)
        )
        tmp[za[ipol] < 0.99] = self._NSide * np.sqrt(3.0 * (1 - za[ipol & (za < 0.99)]))
        jp = np.floor(tp * tmp)
        jm = np.floor((1.0 - tp) * tmp)
        j[ipol] = jp + jm + 1
        i[ipol] = np.floor(tt[ipol] * j[ipol])
        # distinguish between north and south pole
        j[ipol & (z > 0)] -= 1
        j[ipol & (z < 0)] = 4 * self._NSide - j[ipol & (z < 0)] - 1

        return j, i

    ## This function returns the RING indices for the pixels that contain the given
    # positions on the sky.
    #
    def getRINGIndex(self, theta, phi):
        j, i = self.getHEALPixIndices(theta, phi)
        ringIndex = np.zeros(j.shape, dtype=int)

        # north pole
        iNPol = j < self._NSide - 1
        ringIndex[iNPol] = 2 * (j[iNPol] + 1) * j[iNPol] + i[iNPol]

        # equator
        iEq = (j >= self._NSide - 1) & (j < 3 * self._NSide)
        ringIndex[iEq] = (
            2 * (self._NSide ** 2 - self._NSide)
            + 4 * (j[iEq] - self._NSide + 1) * self._NSide
            + i[iEq]
        )

        # south pole
        iSPol = j >= 3 * self._NSide
        j[iSPol] = 4 * self._NSide - 1 - j[iSPol]
        ringIndex[iSPol] = (
            12 * self._NSide ** 2 - 2 * (j[iSPol] + 1) * j[iSPol] + i[iSPol]
        )

        return ringIndex

    ## This function returns arrays of ring indices j and pixel-in-ring indices i that,
    # when traversed in order, traverse the HEALPix grid in RING order.
    #
    # These values are cached internally for improved efficiency.
    #
    def getPixelIndices(self):
        if self._j is None:
            # we first overdo it a bit
            j, i = np.mgrid[0 : 4 * self._NSide - 1 : 1, 0 : 4 * self._NSide : 1]
            # now we mask out the top and bottom triangles for the missing pixels in the polar regions
            select = np.ones(j.shape, dtype=bool)
            select[i + 1 > 4 * (j + 1)] = False
            select[i + 1 > 4 * (4 * self._NSide - j - 1)] = False
            # and remove them from the arrays
            j = j[select]
            i = i[select]

            self._j = j
            self._i = i

        return self._j, self._i

    ## This function returns arrays of zenith angles theta and azimuth angles phi in the
    # same RING order returned by getPixelIndices().
    #
    # These values are cached internally for improved efficiency.
    #
    def getPixelAngles(self):
        if self._theta is None:
            j, i = self.getPixelIndices()
            theta, phi = self.getHEALPixAngles(j, i)
            self._theta = theta
            self._phi = phi

        return self._theta, self._phi

    ## Access the pixels corresponding to the given ring and pixel-in-ring indices.
    #
    def getPixel(self, j, i):
        if len(self._HEALPixCube.shape) == 2:
            return self._HEALPixCube[j, i]
        elif len(self._HEALPixCube.shape) == 3:
            return self._HEALPixCube[:, j, i]

    ## Access the pixels corresponding to the given RING indices.
    #
    def getPixelFromRINGIndex(self, index):
        j, i = self.getPixelIndices()
        return getPixel(self, j[index], i[index])

    ## This function returns the projection of the HEALPix data map, using the given
    # number of vertical pixels for the resulting projection image, and the given projection
    # transformation. The optional parameters thetaCenter and phiCenter allow to select a
    # different central position than the original crosshair of the HEALPixSkyInstrument.
    #
    # The function first sets up the linear coordinates for the projection image and then converts
    # them to the appropriate zenith and azimuth angles using the appropriate projection. These angles
    # are then rotated according to thetaCenter and phiCenter so that they point to the correct angles
    # on the rotated HEALPix sphere. Finally, the angles are passed on to getHEALPixIndices() to get
    # the corresponding HEALPix pixels.
    #
    # Because of the nature of the Mollweide projection, some pixels in the four corners of the output
    # image are unused. The values in these pixels are set to zero.
    #
    # The angles thetaCenter and phiCenter lead to vertical and horizontal rotations respectively. The
    # horizontal rotation is simply to the right. The vertical rotation goes upward in the left part of
    # the image and downward in the right part. Pixels that rotate through the pole move from the left
    # part to the right part. Due to the projection, thetaCenter can go through a full 360 degrees
    # rotation before the projection image is the same; thetaCenter angles larger than 180 degrees
    # correspond to projection images for smaller thetaCenter values that are seen upside down.
    #
    def getProjectionMap(
        self, nPixelY, projection="Mollweide", thetaCenter=0.0, phiCenter=0.0
    ):

        # convert the azimuth angle from degrees to radians
        phiCenter *= np.pi / 180.0
        thetaCenter *= np.pi / 180.0

        # make sure we oversample the HEALPix pixels initially, we will resample after the projection
        if nPixelY < 2 * self._NSide:
            nPixelYHigh = 2 * self._NSide
            nPixelYHigh += nPixelY - nPixelYHigh % nPixelY
            resFac = nPixelYHigh // nPixelY
        else:
            nPixelYHigh = nPixelY
            resFac = 1

        # initialize the image pixels
        nPixelXHigh = 2 * nPixelYHigh
        yMin = -1.0 + 1.0 / nPixelYHigh
        dY = 2.0 / nPixelYHigh
        xMin = -2.0 + 2.0 / nPixelXHigh
        dX = 4.0 / nPixelXHigh
        y, x = np.mgrid[yMin:1.0:dY, xMin:2.0:dX]

        # create an empty pixel map
        if len(self._HEALPixCube.shape) == 2:
            image = np.zeros(x.shape)
        elif len(self._HEALPixCube.shape) == 3:
            image = np.zeros((self._HEALPixCube.shape[0], x.shape[0], x.shape[1]))
        else:
            raise ValueError(
                "Expected HEALPixCube of dimension 2 or 3, but instead got cube with shape {0}".format(
                    self._HEALPixCube.shape
                )
            )

        logging.info(
            "Creating a projected data cube with shape {0}".format(image.shape)
        )

        if projection == "Mollweide":
            # compute the Mollweide angles for each pixel
            temp = np.arcsin(y)
            theta = np.arcsin((2.0 * temp + np.sin(2.0 * temp)) / np.pi) + 0.5 * np.pi
            phi = np.pi + 0.5 * np.pi * x / np.cos(temp)

            # filter out invalid angles in the corners
            iValid = (0.5 * x) ** 2 + (y) ** 2 < 1.0
            theta = theta[iValid]
            phi = phi[iValid]
        elif projection == "HammerAitoff":
            x *= np.sqrt(2.0)
            y *= np.sqrt(2.0)
            # compute the Hammer-Aitoff angles for each pixel
            temp = np.sqrt(1.0 - (0.25 * x) ** 2 - (0.5 * y) ** 2)
            theta = np.arcsin(temp * y) + 0.5 * np.pi
            phi = 2.0 * np.arctan(0.5 * temp * x / (2.0 * temp ** 2 - 1.0)) + np.pi

            # filter out invalid angles in the corners
            iValid = (0.5 * x) ** 2 + (y) ** 2 < 2.0
            theta = theta[iValid]
            phi = phi[iValid]
        else:
            raise UserError(
                "Unknown projection ("
                + projection
                + ")! Possible values are Mollweide, HammerAitoff."
            )

        # invert the azimuth and zenith angles, so that applying this projection yields exactly
        # the same image as would be returned by an AllSkyInstrument with the same projection
        # at the same position
        # The real reason for the inversion is that geographical coordinates assume a different
        # direction than our astrophyical coordinates:
        #  - Colatitude is 0 at the North Pole and increases to 180 degrees at the South Pole,
        #    while geographical latitude is 0 at the equator and increases towars the North Pole
        #  - While both galactic and geographical longitude increase from West to East, astronomers
        #    put the East left and the West right, since they see the sphere from the inside.
        phi = 2.0 * np.pi - phi
        theta = np.pi - theta

        # we have successfully applied the projection to obtain the angular coordinates of each
        # pixel of the projection image on the unit sphere, centred on (theta, phi) = (0, 0)
        # now we need to transform these angles to the unit sphere centred on (thetaCenter, phiCenter)
        # we do this by applying two rotations:
        #  - a rotation over phiCenter along the z axis. This rotation only affects phi
        phi += phiCenter
        #  - a rotation over thetaCenter along the (rotated) x axis. This rotation needs to be done
        #    in Cartesian space
        sinTheta = np.sin(theta)
        x = sinTheta * np.cos(phi)
        y = sinTheta * np.sin(phi)
        z = np.cos(theta)
        sinThetaCenter = np.sin(thetaCenter)
        cosThetaCenter = np.cos(thetaCenter)
        temp = y * cosThetaCenter + z * sinThetaCenter
        z = -y * sinThetaCenter + z * cosThetaCenter
        y = temp
        phi = np.arctan2(y, x)
        theta = np.arccos(z)

        # make sure the final angles are within range for the HEALPix grid
        phi = np.mod(phi, 2.0 * np.pi)
        theta = np.mod(theta, np.pi)

        j, i = self.getHEALPixIndices(theta, phi)

        # copy the HEALPix pixel values into the corresponding image pixels (only for valid pixels)
        if len(self._HEALPixCube.shape) == 2:
            image[iValid] = self._HEALPixCube[j, i]
            # explicitly set the invalid pixels to NaN to mask them
            image[~iValid] = np.nan
        elif len(self._HEALPixCube.shape) == 3:
            image[:, iValid] = self._HEALPixCube[:, j, i]
            # explicitly set the invalid pixels to NaN to mask them
            image[:, ~iValid] = np.nan

        # resample the image onto the desired resolution
        if resFac > 1:
            if len(self._HEALPixCube.shape) == 2:
                temp = image.reshape(
                    (image.shape[0] // resFac, resFac, image.shape[1] // resFac, resFac)
                )
                image = np.sum(temp, axis=(1, 3))
            elif len(self._HEALPixCube.shape) == 3:
                temp = image.reshape(
                    (
                        image.shape[0],
                        image.shape[1] // resFac,
                        resFac,
                        image.shape[2] // resFac,
                        resFac,
                    )
                )
                image = np.sum(temp, axis=(2, 4))
            logging.info("Resampled projected cube to shape {0}".format(image.shape))

        return image

    ## This function returns a new HEALPixGrid that contains a degraded copy of this one.
    # The degraded copy combines 2*2^degradeFactor pixels into a single pixel using the
    # given operator (default is sum) to combine pixel values. This effectively changes
    # the order of the HEALPixGrid by -degradeFactor.
    #
    def degrade(self, degradeFactor, operator=np.sum):

        # compute the properties of the degraded grid
        newOrder = self._order - degradeFactor
        newNSide = 1 << newOrder

        # allocate the degraded data cube
        if len(self._HEALPixCube.shape) == 2:
            newHEALPixCube = np.zeros((4 * newNSide - 1, 4 * newNSide))
        elif len(self._HEALPixCube.shape) == 3:
            newHEALPixCube = np.zeros(
                (self._HEALPixCube.shape[0], 4 * newNSide - 1, 4 * newNSide)
            )

        # create the degraded grid
        newHEALPixGrid = HEALPixGrid(newHEALPixCube)

        # get the indices and angles for the centres of the old pixels
        j, i = self.getPixelIndices()
        theta, phi = self.getPixelAngles()

        # now get the pixels that contain these angles
        newj, newi = newHEALPixGrid.getHEALPixIndices(theta, phi)
        # what we want to do at this point, is something like
        #  newHEALPixCube[newj, newi] += self._HEALPixCube[j, i]
        # (assuming our operator is a summation)
        # this does not work, and that is very sensible, since this is like
        # trying to add to the same variable with multiple threads in a shared
        # memory context (and for all we know, that could be exactly what happens
        # inside the NumPy library)
        # we can however bypass this issue by using our knowledge of the grid:
        # we know that every new pixel corresponds to exactly 2*2^degradeFactor
        # old pixels, so if we could do the summation 2*2^degradeFactor times,
        # and only sum the corresponding fraction of pixels in one go, then
        # we are done
        # this only works if the new pixels are ordered sensibly, which they
        # will not be in our RING pixel ordering (the first 4 indices for
        # example *always* correspond to different pixels, no matter what
        # the resolution of the HEALPix grid is)
        # we can however find a good ordering by computing the pixel index
        # for our new pixels, and argsorting on this. We don't even care
        # about properly accounting for the corners, as we are only
        # interested in the relative index ordering, not in retrieving the
        # correct RING pixel index
        # So:
        #  - compute the RING pixel indices within the degraded grid
        pixel = newj * 4 * newNSide + newi
        #  - argsort these. If we index the j and i arrays with the pixSort
        #    array, then 2*2^degradeFactor consecutive elements in the
        #    resulting arrays will correspond to the same pixel in the degraded
        #    grid
        pixSort = np.argsort(pixel)
        #  - now do the summations
        #    we could do this using a loop, but that would make it harder to
        #    apply our operator (not all operators can be applied sequentially)
        #    and would possibly be inefficient
        #    instead, we do some clever reshaping and directly put the new pixel
        #    values into their corresponding new pixels
        numOldPerNew = (1 << degradeFactor) ** 2
        # we do need the new pixel indices (in the right order) for this to work
        jsmall, ismall = newHEALPixGrid.getPixelIndices()
        # now get the old pixels sorted per new pixel index and reshaped
        # appropriately so that we can simply apply our operator to the last axis
        if len(self._HEALPixCube.shape) == 2:
            sortedPixels = self._HEALPixCube[j[pixSort], i[pixSort]].reshape(
                (12 * newNSide ** 2, numOldPerNew)
            )
            newHEALPixCube[jsmall, ismall] = operator(sortedPixels, axis=1)
        elif len(self._HEALPixCube.shape) == 3:
            sortedPixels = self._HEALPixCube[:, j[pixSort], i[pixSort]].reshape(
                (self._HEALPixCube.shape[0], 12 * newNSide ** 2, numOldPerNew)
            )
            newHEALPixCube[:, jsmall, ismall] = operator(sortedPixels, axis=2)

        return newHEALPixGrid

    ## This function returns an array containing the position of each pixel in
    # 3D space, assuming a celestial sphere with radius 1.
    #
    # These values are cached internally for improved efficiency.
    #
    def getPositions(self):
        if self._positions is None:
            theta, phi = self.getPixelAngles()
            self._positions = np.zeros((theta.shape[0], 3))
            sinTheta = np.sin(theta)
            self._positions[:, 0] = sinTheta * np.cos(phi)
            self._positions[:, 1] = sinTheta * np.sin(phi)
            self._positions[:, 2] = np.cos(theta)

        return self._positions

    ## This function returns a scipy.spatial.cKDTree containing all pixels in the
    # grid. This tree can be used to perform spatial queries on the pixels (in
    # Cartesian 3D space).
    #
    # The tree is cached fro improved efficiency.
    #
    def getKDTree(self):
        if self._tree is None:
            logging.info(
                "Creating KD tree to speed up spatial queries."
                " This might take a while..."
            )
            positions = self.getPositions()
            self._tree = cKDTree(positions)
            logging.info("Done creating KD tree.")

        return self._tree

    ## Get a mask that selects all pixels within an annulus with the given inner
    # and outer radius surrounding the given central pixel.
    #
    # When indexed with this mask, any grid related array will be reduced to an array
    # that only contains values for pixels within the annulus.
    #
    # When the optional argument useTree is set to True, the function will make use of
    # a KD-tree (scipy.spatial.cKDTree) to speed up the neighbour finding. This is only
    # useful for (many) repeated calls to this function and for large grids. In this case,
    # the first call to this function will take a lot longer, because it will need to
    # construct the tree. The tree is then cached for later reuse.
    #
    def getAnnulusAroundPixel(
        self, centralPixel, innerRadius, outerRadius, useTree=False
    ):
        theta, phi = self.getPixelAngles()
        if useTree:
            positions = self.getPositions()
            tree = self.getKDTree()
            annulus = np.array(
                tree.query_ball_point(positions[centralPixel], outerRadius)
            )
            delta = getDistanceOnSky(
                theta[centralPixel], phi[centralPixel], theta[annulus], phi[annulus]
            )
            selection = (delta >= innerRadius) & (delta <= outerRadius)
            annulus = annulus[selection]
            return annulus
        else:
            delta = getDistanceOnSky(theta[centralPixel], phi[centralPixel], theta, phi)
            selection = (delta >= innerRadius) & (delta <= outerRadius)
            return selection


# -----------------------------------------------------------------

## Construct a new, empty HEALPixGrid of the given order and with the given number
# of values per pixel.
#
# The order \f$k\f$ relates to the HEALPix parameter \f$N_{side} = 2^k\f$, which gives
# the number of pixels along one axis of the 12 HEALPix base pixels. The total number
# of pixels is then given by \f$12N_{side}^2\f$. Each pixel has the exact same size,
# \f$\frac{\pi{}}{3N_{side}^2}\f$. Assuming square pixels (which is approximately true),
# the angular resolution \f$R\f$ can be computed from
# \f[
#    R \approx{} \sqrt{\frac{\pi{}}{3N_{side}^2}}
#      = \frac{1}{N_{side}} \left(\frac{10800}{\sqrt{3\pi{}}}\right) '
#      \approx{} \left(\frac{3520}{N_{side}}\right) '
# \f]
#
# For example, most Planck maps are published with an order \f$k=10\f$, which
# corresponds to a resolution of \f$\approx{}3.4'\f$ (nominally \f$5'\f$). \f$k=6\f$ would
# correspond to a resolution of \f$\approx{}1^\circ{}\f$.
#
def getEmptyHEALPixGrid(order, numValuesPerPixel=1):
    NSide = 1 << order
    if numValuesPerPixel == 1:
        dataCube = np.zeros((4 * NSide - 1, 4 * NSide))
    else:
        dataCube = np.zeros((numValuesPerPixel, 4 * NSide - 1, 4 * NSide))
    return HEALPixGrid(dataCube)


# -----------------------------------------------------------------
