#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visual.rgbimage Managing a color image in various representations
#
# An instance of the RGBImage class in this module manages a color image in one or more representations.

# -----------------------------------------------------------------

import astropy.io.fits as fits
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import PIL.Image
import pts.utils as ut

# -----------------------------------------------------------------
#  RGBImage class
# -----------------------------------------------------------------

## An instance of the RGBImage class represents a single RGB color image, i.e. a 2D array of pixels with three
# color values per pixel. The class allows arbitrary floating point color values, and supports conversion to and
# from more standard image representations using 8-bit or 16-bit integer values. The floating point representation
# offers a wide dynamic range and enables computations on pixel values, such as rescaling by applying a log function.
# It is the natural representation for working with color planes loaded from frames in a FITS file. The integer
# representation is relevant to standard image processing, for example to load or store a color image as a PNG file.
#
# RGBImage properties include the shape of the image (number of horizontal pixels, number of vertical pixels),
# the range of the pixel values (minimum value, maximum value), and some representation of the pixel data itself.
# To avoid unnecessary conversions, and to keep the image representation as general as possible, an RGBImage
# maintains several variants of the pixel data. Conversions between these data types happen only when needed
# and are transparent to the user.
#
# The implementation maintains the following variants of the pixel data:
#  - a PIL.Image instance representing a standard pixel image, used to support various image file formats;
#  - a buffer with RGBA pixels, as can be extracted from a matplotlib Figure instance;
#  - a numpy array with shape (nx, ny, 3) of floating point value type.
#
# An RGBImage can be created from the following sources:
#  - an image file in one of these formats: JPEG, PNG, TIFF;
#  - a FITS file with one or more data frames;
#  - a matplotlib figure generated using an AGG back-end;
#  - a 3D numpy array with shape (nx, ny, 3).
#
# An RGBImage can be repurposed as follows:
#  - save to an image file in one of these formats: JPEG, PNG, TIFF, PDF;
#  - save to a FITS file with three data frames;
#  - plot to the current matplotlib figure;
#  - get pixeldata as a 3D numpy array with shape (nx, ny, 3).
#
# Other functions include:
#  - apply the log function to the pixel values and their range;
#  - change the pixel range without scaling the values, but clipping the values where needed;
#  - scale pixel values to a new range (this happens automatically during some conversions).
#
# There are several private properties:
#  - shape: a 2-tuple of integer numbers specifying the number of horizontal and vertical pixels in the image;
#  - data properties dpil, dbuf, darr: representation of the pixel data, each in its own specific way;
#  - rangearr: a 2-tuple of floating point numbers specifying the minimum and maximum allowed pixel value.
#
# A data property either holds a relevant value or it set to None. At all times, at least one data property holds
# a relevant value; if multiple data properties hold a relevant value, they represent the exact same image.
#
# The allowed range for pixel values in the dpil and dbuf representations is fixed to 0..255
# and thus does not need to be explicitly stored. For pixel data in the form of a numpy array however,
# the allowed range is flexible, and may be larger than the range actually used by pixel values.
# Thus the allowed range cannot be derived from darr.min()/darr.max() and must be explicitly stored.
# In other words, if darr is not None, then rangearr specifies the minimum and maximum allowed pixel values,
# and if darr is None, rangearr is also None.
#
class RGBImage:
    # ---------- Constructing from source -----------------------------

    ## The constructor creates a new RGBImage from the source specified by the sole argument, as follows:
    #  - A pathlib.Path or string with the filepath of an image file in one of the supported formats (JPEG, PNG, TIFF);
    #    the filename extension \em must be one of .jpg, .jpeg, .png, .tif, or .tiff (case insensitive).
    #  - A pathlib.Path or string with the filepath of a FITS file; the filename extension \em must be .fits
    #    (case insensitive). By default, matching the SKIRT convention of storing data cube frames in order of
    #    increasing wavelength, the R channel receives the last frame in the fits file, the G channel the middle frame,
    #    and the B channel the first frame. If there is only one frame, the resulting image is effectively greyscale.
    #    If there are more than three frames, the RGB image is constructed from widely separated frames.
    #    You can override this default behavior by explicitly specifying the frames to be loaded in the RGB
    #    channels with the \em frameIndices argument. Provide a 3-tuple (R,G,B) containing the zero-based frame indices
    #    for respectively the red, green and blue channels; for example: <tt>frameIndices=(34,12,0)</tt>.
    #  - A matplotlib Figure instance generated using an AGG back-end (use plt.get_backend() to verify).
    #  - A three-dimensional numpy array with shape (nx, ny, 3), using any numeric value type.
    #    The first two indices of the array respectively correspond to the x and y axes, so that
    #    arr[0,0,:] represents the lower left pixel, and arr[nx-1,ny-1,:] the upper right pixel. Furthermore
    #    arr[:,:,0] is the red plane, arr[:,:,1] the green plane, and arr[:,:,2] the blue plane.
    #
    def __init__(self, source, frameIndices=None):
        # initialize private data members here to avoid warnings in the development environment
        self._shape = None
        self._dpil = None
        self._dbuf = None
        self._darr = None
        self._rangearr = None

        # numpy array
        if isinstance(source, np.ndarray) and source.ndim == 3 and source.shape[2] == 3:
            self._shape = source.shape[0:2]
            self._setArr(source)

        # matplotlib Figure
        elif isinstance(source, matplotlib.figure.Figure):
            source.canvas.draw()    # flush the figure to the buffer
            buffer, self._shape = source.canvas.print_to_buffer()
            self._setBuf(buffer)

        # file path
        elif isinstance(source, (str,pathlib.Path)):
            path = ut.absPath(source)

            # standard image file
            if path.suffix.lower() in (".jpg",".jpeg",".png",".tif",".tiff"):
                self._setPil(PIL.Image.open(path))
                self._shape = self._dpil.size

            # FITS file
            elif path.suffix.lower() == ".fits":
                data = fits.getdata(path).T  # returns an array with shape (nx, ny) or (nx, ny, nlambda)
                if data.ndim == 2:
                    self._setArr(np.dstack((data, data, data)))
                elif data.ndim == 3:
                    if frameIndices is None:
                        n = data.shape[2]
                        frameIndices = (n-1, n//2, 0)
                    r,g,b = frameIndices
                    self._setArr(np.dstack((data[:, :, r], data[:, :, g], data[:, :, b])))
                else:
                    raise ValueError("Data in FITS file '{}' has unsupported shape".format(path))
                self._shape = self._darr.shape[0:2]

        # unsupported type
        else:
            raise ValueError("Image source '{}' has unsupported type".format(source))

    # ---------- Basic attributes -------------------------------------

    ## This function returns the image's shape as a tuple (nx, ny).
    def shape(self):
        return self._shape

    ## This function returns the range of the image's pixel values as a tuple (min, max).
    def pixelRange(self):
        if self._darr is None: return 0, 255.999
        else: return self._rangearr

    ## This function determines the range of the image's pixel values, ignoring zero and negative values, and omitting
    # outliers on both sides of the scale as specified by the \em from and \em to percentiles. The resulting range is
    # returned as a tuple (min, max).
    def percentilePixelRange(self, from_percentile=0, to_percentile=100):
        self._ensureArr(invalidate=False)
        # determine the range, ignoring zero and negative values, and removing the outliers
        nonzero = self._darr[self._darr > 0]
        if len(nonzero)==0: return 0, 0
        smallest = np.percentile(nonzero, from_percentile)
        largest = np.percentile(nonzero, to_percentile)
        return smallest, largest

    # ---------- Repurposing to target --------------------------------

    ## This function saves the image to a file with the specified /em outFilePath in one of the supported formats.
    # The file format is chosen based on the filename extension in the /em outFilePath, which /em must be one of
    # .jpg, .jpeg, .png, .tif, .tiff, .pdf or .fits (case insensitive). The image is saved in RGB mode.
    # If the \em tiff16bit argument is \c True and the filename extension is .tif or .tiff, the image
    # is saved with 16 bits per sample rather than the standard 8 bits per sample; however in this case
    # the image data is not compressed.
    # JPEGs are compressed with quality 80 on a scale of 1 to 100. FITS files contain a data cube with
    # three data frames in the order BGR, matching the SKIRT convention of storing data cube frames in
    # order of increasing wavelength.
    def saveTo(self, outFilePath, tiff16bit=False):
        path = ut.absPath(outFilePath)
        suffix = path.suffix.lower()

        # 16-bit TIFF file
        if tiff16bit and suffix in (".tif",".tiff"):
            # tiff header as uint16 words
            lsNX,msNX = _split16(self._shape[0])
            lsNY,msNY = _split16(self._shape[1])
            lsBYTES,msBYTES = _split16(self._shape[0] * self._shape[1] * 6)
            header = (
                0x4949, 42, 8, 0,                  #   0: TIFF header (little endian)
                12,                                #   8: number of directory entries
                                                   #  (directory entry: tag,type,count,0,value/offset x 2)
                256, 4, 1, 0, lsNX, msNX,          #  10: ImageWidth, 1 LONG
                257, 4, 1, 0, lsNY, msNY,          #  22: ImageLength, 1 LONG
                258, 3, 3, 0, 158, 0,              #  34: BitsPerSample, 3 SHORT (-> offset!)
                259, 3, 1, 0, 1, 0,                #  46: Compression, 1 SHORT
                262, 3, 1, 0, 2, 0,                #  58: PhotometricInterpretation, 1 SHORT
                273, 4, 1, 0, 180, 0,              #  70: StripOffsets, 1 LONG
                277, 3, 1, 0, 3, 0,                #  82: SamplesPerPixel, 1 SHORT
                278, 4, 1, 0, lsNY, msNY,          #  94: RowsPerStrip, 1 LONG
                279, 4, 1, 0, lsBYTES, msBYTES,    # 106: StripByteCounts, 1 LONG
                282, 5, 1, 0, 164, 0,              # 118: XResolution, 1 RATIONAL (-> offset!)
                283, 5, 1, 0, 172, 0,              # 130: YResolution, 1 RATIONAL (-> offset!)
                296, 3, 1, 0, 2, 0,                # 142: ResolutionUnit, 1 SHORT
                0, 0,                              # 154: IFD list terminator
                16, 16, 16,                        # 158: BitsPerSample value
                72, 0, 1, 0,                       # 164: XResolution value
                72, 0, 1, 0 )                      # 172: YResolution value
                                                   # 180: Image data
            out = open(path, 'wb')
            out.write(np.array(header,dtype=np.uint16).tobytes())
            data = self.scaledPixelArray(0,65535.99)
            out.write(np.flipud(np.rollaxis(data,1)).astype(np.uint16).tobytes())
            out.close()

        # standard 8-bit image file
        elif suffix in (".bmp",".gif",".jpg",".jpeg",".png",".tif",".tiff",".pdf"):
            self._ensurePil(invalidate=False)
            if suffix in (".jpg",".jpeg"): self._dpil.save(path, format="JPEG", quality=80)
            elif suffix in (".tif",".tiff"): self._dpil.save(path, format="TIFF")
            elif suffix == ".png": self._dpil.save(path, format="PNG")
            elif suffix == ".pdf": self._dpil.save(path, format="PDF")

        # FITS file
        elif suffix == ".fits":
            self._ensureArr(invalidate=False)
            data = np.dstack((self._darr[:, :, 2], self._darr[:, :, 1], self._darr[:, :, 0]))
            fits.writeto(path, data.T, overwrite=True)

        # unsupported type
        else:
            raise ValueError("OutFilePath has unsupported filename extension")

    ## This function plots the image to the current axes in the current matplotlib figure. The image is
    # rescaled to fit the axes. If \em fill is True, the image is stretched to fit the axes in both directions
    # changing the image aspect ratio if needed. If \em fill is False (the default), the axes aspect ratio
    # is adjusted so that the image aspect ratio is preserved.
    def plot(self, fill=False):
        # if we have no array but we have a PIL image, use the PIL image; otherwise use the array
        if self._darr is None and self._dpil is not None:
            plt.imshow(self._dpil, aspect='auto' if fill else 'equal', interpolation='bicubic', origin='upper')
        else:
            data = np.rollaxis(self.scaledPixelArray(0,1), 1)
            plt.imshow(data, aspect='auto' if fill else 'equal', interpolation='bicubic', origin='lower')

    ## This function returns the image's pixeldata as a 3D numpy array with shape (nx, ny, 3).
    def pixelArray(self):
        self._ensureArr(invalidate=False)
        return self._darr

    ## This function returns the image's pixel data as a 3D numpy array with shape (nx, ny, 3),
    # after scaling the data to the specified range, but without affecting the data in the image itself.
    def scaledPixelArray(self, newmin, newmax):
        self._ensureArr(invalidate=False)
        if newmin != self._rangearr[0] or newmax != self._rangearr[1]:
            scale = float(newmax-newmin) / float(self._rangearr[1] - self._rangearr[0])
            offset = newmin - scale*self._rangearr[0]
            result = offset + scale*self._darr
            np.clip(result, newmin, newmax, out=result)  # clip away out-of-range values due to inaccuracies
            return result
        else:
            return self._darr

    ## This function returns the image's pixel data as a bytes object in the rgba format that can also be
    # obtained from mtplotlib (3 RGB channels plus a dummy alpha channel, 8 bits each).
    def bytesArray(self):
        self._ensureBuf(invalidate=False)
        return self._dbuf

    # ---------- Adjusting the image ----------------------------------

    ## This function applies the natural logarithm function to the pixel values of the image, and adjusts
    # the pixel range accordingly. This function raises an error if the pixel range includes zero or negative
    # values when this function is called. This can be avoided, for example, by calling the setrange()
    # function with a nonzero minimum value before calling the applylog() function.
    def applyLog(self):
        self._ensureArr(invalidate=True)
        # verify that all values are positive
        if self._rangearr[0] <= 0: raise ValueError("Log can't be applied to negative or zero pixel values")
        # apply logarithm to array values and range
        self._darr = np.log(self._darr)
        self._rangearr = (np.log(self._rangearr[0]), np.log(self._rangearr[1]))

    ## This function applies the square root function to the pixel values of the image, and adjusts
    # the pixel range accordingly. This function raises an error if the pixel range includes zero or negative
    # values when this function is called. This can be avoided, for example, by calling the setrange()
    # function with a nonzero minimum value before calling the applylog() function.
    def applySqrt(self):
        self._ensureArr(invalidate=True)
        # verify that all values are positive
        if self._rangearr[0] <= 0: raise ValueError("Sqrt can't be applied to negative or zero pixel values")
        # apply logarithm to array values and range
        self._darr = np.sqrt(self._darr)
        self._rangearr = (np.sqrt(self._rangearr[0]), np.sqrt(self._rangearr[1]))

    ## This function applies a transformation defined by a cubic spline curve to the (scaled) pixel values,
    # and replaces the image data with the result. The pixel values are scaled to the interval [0,1].
    # The curve is composed of three cubic spline segments connecting the four points
    # \f$(0,0),\,(x_1,y_1),\,(x_2,y_2),\,(1,1)\f$ with \f$0<x_1<x_2<1\f$,
    # where the outer points are fixed and the inner points are provided as arguments to this function.
    def applyCurve(self, point1=(0.25,0.16), point2=(0.80,0.86)):
        self.scaleValues(0., 1.)
        curve = _CubicSpline(point1, point2)
        self._darr = curve.ay(self._darr)

    ## This function applies the specified color map to the (scaled) pixel values of the red channel,
    # and replaces the image data with the result. The \em cmap argument can be one of the following:
    #  - the name of a standard matplotlib color map (see figure below), possibly
    #    followed by "_r" to get the reversed map;
    #  - a matplotlib ColorMap object;
    #  - any callable object taking a single array argument with floating point values in range [0,1]
    #    and returning a tuple of four RGBA values in range [0,1] for each element in the array -- in
    #    other words newshape = oldshape+(4,) -- where the A value is ignored.
    #
    # \image html cmaps.png
    #
    def applyColorMap(self, cmap):
        if not hasattr(cmap, "__call__"): cmap = plt.get_cmap(cmap)
        self.scaleValues(0., 1.)
        self._darr = np.delete(cmap(self._darr[:, :, 0]), 3, 2)

    ## This function changes the pixel value range without affecting the actual pixel values except that,
    #  if the new range is smaller than the previous range, the pixel values are clipped to the new range.
    #  One can specify a new minimum, a new maximum, or both.
    def setRange(self, newmin=None, newmax=None):
        self._ensureArr(invalidate=True)
        if newmin is None: newmin = self._rangearr[0]
        if newmax is None: newmax = self._rangearr[1]
        if newmin > self._rangearr[0] or newmax < self._rangearr[1]:
            np.clip(self._darr, newmin, newmax, out=self._darr)
        self._rangearr = (newmin, newmax)

    ## This function applies a linear transformation to the pixel values and to the pixel value range,
    #  such that the previous range is exactly mapped onto the new specified range.
    def scaleValues(self, newmin, newmax):
        self._ensureArr(invalidate=True)
        self._darr = self.scaledPixelArray(newmin, newmax)
        self._rangearr = (newmin, newmax)

    # ---------- Combining images -------------------------------------

    ## This function places the specified RGBImage to the right of the receiving image.
    # It is the caller's responsibility to ensure that both images have the same height (in pixels).
    def addRight(self, image):
        data = np.vstack(( self.scaledPixelArray(0.,1.), image.scaledPixelArray(0.,1.) ))
        self._setArr(data)
        self._shape = self._darr.shape[0:2]

    ## This function places the specified RGBImage below the receiving image.
    # It is the caller's responsibility to ensure that both images have the same width (in pixels).
    def addBelow(self, image):
        data = np.hstack(( image.scaledPixelArray(0.,1.), self.scaledPixelArray(0.,1.) ))
        self._setArr(data)
        self._shape = self._darr.shape[0:2]

    ## This function adds pixel rows and/or columns to the receiving image if needed to match a given shape.
    # Each direction (horizontal and vertical) is handled separately. Extra pixels are added equally on both sides
    # of the image (left/right or top/bottom), and are given a value of zero. Pixels are never removed.
    # The target shape can be specified as a 2-tuple (nx,ny), or by providing another RGBImage object.
    def enlargeCanvas(self, shape):
        if isinstance(shape, RGBImage): shape = shape.shape()
        xdif = shape[0]-self._shape[0]
        ydif = shape[1]-self._shape[1]
        if xdif > 0 or ydif > 0:
            myArr = self.pixelArray()
            if xdif > 0:
                arrLeft  = np.zeros(( int(xdif/2), myArr.shape[1], 3 ))
                arrRight = np.zeros(( xdif - int(xdif/2), myArr.shape[1], 3))
                myArr = np.vstack(( arrLeft, myArr, arrRight ))
            if ydif > 0:
                arrUp   = np.zeros(( myArr.shape[0], int(ydif/2), 3))
                arrDown = np.zeros(( myArr.shape[0], ydif - int(ydif/2), 3))
                myArr = np.hstack(( arrDown, myArr, arrUp))
            self._setArr(myArr)
            self._shape = self._darr.shape[0:2]

    # ---------- Private utility functions ----------------------------

    ## This private function sets the specified PIL image as internal representation, invalidating
    #  the other representations. The image is converted to RGB mode if needed.
    def _setPil(self, pil):
        self._dpil = pil.convert("RGB")
        self._dbuf = None
        self._darr = None
        self._rangearr = None

    ## This private function sets the specified RGBA buffer or string as internal representation,
    #  invalidating the other representations.
    def _setBuf(self, buf):
        self._dpil = None
        self._dbuf = buf
        self._darr = None
        self._rangearr = None

    ## This private function sets the specified numpy array as internal representation, invalidating
    #  the other representations. The array values are converted to 8-byte floating point if needed.
    def _setArr(self, arr):
        self._dpil = None
        self._dbuf = None
        self._darr = arr.astype(np.float64)
        self._rangearr = (self._darr.min(), self._darr.max())

    ## This private function ensures that there is a valid PIL image representation, converting from
    #  one of the other representations if necessary, and invalidating the other representations if requested.
    def _ensurePil(self, invalidate=True):
        if self._dpil is None:
            if self._dbuf is not None:
                self._dpil = PIL.Image.frombytes("RGB", self._shape, self._dbuf, "raw", "RGBX", 0, 1)
            elif self._darr is not None:
                data = self.scaledPixelArray(0,255.999)
                buf = np.rollaxis(data,1).astype(np.uint8).tobytes()
                self._dpil = PIL.Image.frombytes("RGB", self._shape, buf, "raw", "RGB", 0, -1)
            else:
                raise ValueError("No source data for conversion to PIL image")
        if invalidate:
            self._dbuf = None
            self._darr = None
            self._rangearr = None

    ## This private function ensures that there is a valid buffer representation, converting from
    #  one of the other representations if necessary, and invalidating the other representations if requested.
    def _ensureBuf(self, invalidate=True):
        if self._dbuf is None:
            if self._dpil is not None:
                self._dbuf = self._dpil.tobytes("raw", "RGBX", 0, 1)
            elif self._darr is not None:
                data = self.scaledPixelArray(0,255.999)
                self._dbuf = np.dstack((np.flipud(np.rollaxis(data, 1)).astype(np.uint8),
                                        np.zeros(self._shape[::-1], np.uint8))).tobytes()
            else:
                raise ValueError("No source data for conversion to buffer")
        if invalidate:
            self._dpil = None
            self._darr = None
            self._rangearr = None

    ## This private function ensures that there is a valid numpy array representation, converting from
    #  one of the other representations if necessary, and invalidating the other representations if requested.
    def _ensureArr(self, invalidate=True):
        if self._darr is None:
            if self._dpil is not None:
                self._darr = np.frombuffer(self._dpil.tobytes("raw", "RGB", 0, -1), np.uint8).astype(np.float64)
                self._darr = np.rollaxis(np.reshape(self._darr, (self._shape[1], self._shape[0], 3)), 1)
            elif self._dbuf is not None:
                self._darr = np.frombuffer(self._dbuf, np.uint8).astype(np.float64)
                self._darr = np.delete(np.reshape(self._darr, (self._shape[1], self._shape[0], 4)), 3, 2)
                self._darr = np.rollaxis(np.flipud(self._darr), 1)
            else:
                raise ValueError("No source data for conversion to array")
            self._rangearr = (0, 255.999)
        if invalidate:
            self._dpil = None
            self._dbuf = None

# -----------------------------------------------------------------

## This helper function returns a 2-tuple containing the least and most significant 16-bit portion
# of the specified unsigned 32-bit integer value.
def _split16(value32):
    return value32 & 65535, (value32 >> 16) & 65535

# -----------------------------------------------------------------
#  CubicSpline class
# -----------------------------------------------------------------

## An instance of the CubicSpline class implements a function \f$y=f(x)\f$ composed of three cubic spline segments
# of the form \f$y=ax^3+bx^2+cx+d\f$.
# The function is constrained as follows:
#  - the function is defined only for \f$x\f$ values in the interval [0,1]
#  - the function goes through the four points \f$(0,0),\,(x_1,y_1),\,(x_2,y_2),\,(1,1)\f$ with \f$0<x_1<x_2<1\f$,
#    where the outer points are fixed and the inner points are provided as arguments to the constructor
#  - the function's first and second derivatives are continuous at the inner points
#  - the function's second derivative is zero at the outer points ("natural bounding conditions")
#  - regardless of the form dictated by the spline segments, the \f$y\f$ value is always clipped to the interval [0,1]
#
class _CubicSpline:
    # This constructor calculates the coefficients for the three spline segments based on the specified points
    # and on the constraints described in the documentation for this class. The points must be given as (x,y) tuples.
    def __init__(self, point1, point2):
        # copy the specified values for convenience and verify the constraints
        x1,y1 = point1
        x2,y2 = point2
        if not (0 < x1 < x2 < 1): raise ValueError("invalid x-value(s)")
        # remember the segment breaks
        self.x1 = x1
        self.x2 = x2
        # determine the coefficients for the three segments
        # (12 linear equations for 12 unknowns, 2 of which are trivially zero and thus excluded from the matrix)
        x1p2 = x1*x1
        x1p3 = x1p2*x1
        x2p2 = x2*x2
        x2p3 = x2p2*x2
        mat = [ [   x1p3,     x1,      0,      0,      0,      0,      0,      0,      0,      0 ],
                [      0,      0,   x1p3,   x1p2,     x1,      1,      0,      0,      0,      0 ],
                [      0,      0,   x2p3,   x2p2,     x2,      1,      0,      0,      0,      0 ],
                [      0,      0,      0,      0,      0,      0,   x2p3,   x2p2,     x2,      1 ],
                [      0,      0,      0,      0,      0,      0,      1,      1,      1,      1 ],
                [ 3*x1p2,      1,-3*x1p2,  -2*x1,     -1,      0,      0,      0,      0,      0 ],
                [      0,      0, 3*x2p2,   2*x2,      1,      0,-3*x2p2,  -2*x2,     -1,      0 ],
                [   6*x1,      0,  -6*x1,     -2,      0,      0,      0,      0,      0,      0 ],
                [      0,      0,   6*x2,      2,      0,      0,  -6*x2,     -2,      0,      0 ],
                [      0,      0,      0,      0,      0,      0,      6,      2,      0,      0 ] ]
        vec = [ y1, y1, y2, y2, 1, 0, 0, 0, 0, 0 ]
        sol = np.linalg.solve(np.array(mat,dtype=np.float64), np.array(vec,dtype=np.float64))
        self.a1 = sol[0]
        self.b1 = 0.
        self.c1 = sol[1]
        self.d1 = 0.
        self.a2 = sol[2]
        self.b2 = sol[3]
        self.c2 = sol[4]
        self.d2 = sol[5]
        self.a3 = sol[6]
        self.b3 = sol[7]
        self.c3 = sol[8]
        self.d3 = sol[9]

    # This function returns the value of the represented function \f$y=f(x)\f$ where x is a scalar value.
    def y(self, x):
        if not (0 <= x <= 1): raise ValueError("invalid x-value")
        xp2 = x*x
        xp3 = xp2*x
        if x < self.x1:   result = self.a1*xp3 + self.b1*xp2 + self.c1*x + self.d1
        elif x < self.x2: result = self.a2*xp3 + self.b2*xp2 + self.c2*x + self.d2
        else:             result = self.a3*xp3 + self.b3*xp2 + self.c3*x + self.d3
        return max(0,min(1,result))

    # This function returns the value of the represented function \f$y=f(x)\f$ where x is a numpy array.
    def ay(self, x):
        if np.any(x<0) or np.any(x>1): raise ValueError("invalid x-value")
        xp2 = x*x
        xp3 = xp2*x
        mask1 = x < self.x1
        mask3 = x > self.x2
        mask2 = np.logical_and(~mask1,~mask3)
        result = np.empty_like(x)
        result[mask1] = self.a1*xp3[mask1] + self.b1*xp2[mask1] + self.c1*x[mask1] + self.d1
        result[mask2] = self.a2*xp3[mask2] + self.b2*xp2[mask2] + self.c2*x[mask2] + self.d2
        result[mask3] = self.a3*xp3[mask3] + self.b3*xp2[mask3] + self.c3*x[mask3] + self.d3
        return np.clip(result,0.,1.,out=result)

# -----------------------------------------------------------------
