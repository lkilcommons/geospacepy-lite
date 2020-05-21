# Geospacepy #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)

### What is this repository for? ###

Geospacepy is a small library of python functions for doing space science data analysis. It is also a dependency of several other applications.

### Rules of the Road ###

Geospacepy is MIT open source licensed software. The usual legal no-warranties, no-guarantees provisions apply. Though this software is indended for use in academic research it has not been peer-reviewed, there may be errors in calculations.  

### OMNIREADER ###

OmniReader has been moved to it's own seperate package (https://github.com/lkilcommons/nasaomnireader)

### SPECIALDATETIME ###

This is a very basic set of utilities for handling numpy arrays of Python datetimes, and converting from and to other time formats, such as MATLAB datenumber, day-of-year, and Julian date. There are other tools which do this, but specialdatetime.py has the advantage of being very lightweight.

### SATPLOTTOOLS ###

This is my plotting library. It builds on the functionality provided by matplotlib. It has methods for doing things such as polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and cartesian.

### Installation Instructions ###

* Clone the repository
* `python setup.py install`

### Who do I talk to? ###

* This repository was created and is managed by Liam M. Kilcommons at CU Boulder
