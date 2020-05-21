# Geospacepy #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)

### What is this repository for? ###

Geospacepy is a small library of python functions for doing space science data analysis. It is also a dependency of several other applications.

### Rules of the Road ###

Geospacepy is MIT open source licensed software. The usual legal no-warranties, no-guarantees provisions apply. Though this software is indended for use in academic research it has not been peer-reviewed, there may be errors in calculations.  

### OMNIREADER ###

OmniReader used to be part of this package, but it has been moved to it's own seperate package (https://github.com/lkilcommons/nasaomnireader)

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

Calculations for positions and measurements on the earth and in near-earth space. 

* Coordinate transformations for earth centered inertial (ECI) and earth centered earth fixed (ECEF) coordinate frames
* Coordinate transformations for geographic latitudes/longitudes assuming a spherical earth

### SPHERICAL_GEOMETERY ###

Calculations (distances,areas,integrals) over the surface of a sphere

### SATPLOTTOOLS ###

Plotting library building on the functionality provided by matplotlib. Methods for polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and plot coorindates.

### Installation Instructions ###

* Clone the repository
* `python setup.py install`

### Who do I talk to? ###

* This repository was created and is managed by Liam M. Kilcommons at CU Boulder
