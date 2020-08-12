# Geospacepy-lite #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)
[![Documentation Status](https://readthedocs.org/projects/geospacepy-lite/badge/?version=stable)](https://geospacepy-lite.readthedocs.io/en/stable/?badge=stable)
[![Coverage Status](https://coveralls.io/repos/github/lkilcommons/geospacepy-lite/badge.svg?branch=master&kill_cache=1)](https://coveralls.io/github/lkilcommons/geospacepy-lite?branch=master)

A small library of python functions useful for geospace data analysis.

Click the 'docs' badge to see the documentation, for interactive examples go to the [notebooks repo](https://github.com/lkilcommons/geospacepy-notebooks).

## Modules in Geospacepy-lite ###

### OMNIREADER ###

Moved to it's own seperate package (https://github.com/lkilcommons/nasaomnireader)

### SPECIALDATETIME ###

Utilities for handling numpy arrays of Python datetimes, and converting from and to other time formats, such as MATLAB datenumber, day-of-year, and Julian date.

### SUN ###

Calculate various solar-position-related parameters

* Solar position (implements the Astronomical Almanac Low-Accuracy Algorithm)
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

Plotting library building on the functionality provided by matplotlib. Methods for polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and plot coordinates.

### Installation Instructions ###

* Clone the repository
* `python setup.py install`

### Who do I talk to? ###

* This repository was created and is managed by Liam M. Kilcommons at CU Boulder
