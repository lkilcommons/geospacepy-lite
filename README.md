# Geospacepy-lite #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)

A small library of python functions for doing space science data analysis

### Rules of the Road ###

Geospacepy-lite is MIT open source licensed software. The usual legal no-warranties, no-guarantees provisions apply. Many caluclations have been verified in the unit tests using literature reference values, but not all algorithms have or can be tested in this way.

## Modules in Geospacepy-lite ###

### OMNIREADER ###

Moved to it's own seperate package (https://github.com/lkilcommons/nasaomnireader)

### SPECIALDATETIME ###

Utilities for handling numpy arrays of Python datetimes, and converting from and to other time formats, such as MATLAB datenumber, day-of-year, and Julian date.

### SUN ###

Calculate various solar-position-related parameters

* Solar position (imlpements the Astronomical Almanac Low-Accuracy Algorithm)
* Greenwich hour angle
* Local hour angle
* Local mean solar time
* Solar zenith angle

### ROTATIONS ###

Rotate coordinate axes of any number of 3 component vectors

### TERRESTRIAL_SPHERICAL ###

* Coordinate transformations for earth centered inertial (ECI) and earth centered earth fixed (ECEF) coordinate frames
* Coordinate transformations for geographic latitudes/longitudes assuming a spherical earth

### SPHERICAL_GEOMETERY ###

Calculations on the surface of a sphere:

* Great circle distance between two lat/lon points
* Midpoint between two lat/lon points on a great circle
* Area of a spherical 'rectangle'
* Surface integral of data on a regular lat/lon grid

### SATPLOTTOOLS ###

Plotting library building on the functionality provided by matplotlib. Methods for polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and plot coorindates.

### Installation Instructions ###

* Clone the repository
* `python setup.py install`

### Who do I talk to? ###

* This repository was created and is managed by Liam M. Kilcommons at CU Boulder
