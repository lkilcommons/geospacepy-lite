# Geospacepy #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)

### What is this repository for? ###

Geospacepy is a small library of python functions for doing space science data analysis. It is also a dependency of several other applications.

### Rules of the Road ###

Geospacepy is MIT open source licensed software. The usual legal no-warranties, no-guarantees provisions apply. Though this software is indended for use in academic research it has not been peer-reviewed, there may be errors in calculations.  

### Python Version WARNING ###

Geospacepy was developed using Python 2.7, Python 3 support was added programmatically using 2to3 script. I haven't tested extensively using Python 3. Unit tests pass but test coverage is pretty spotty still.

### SPECIALDATETIME ###

This is a very basic set of utilities for handling numpy arrays of Python datetimes, and converting from and to other time formats, such as MATLAB datenumber, day-of-year, and Julian date. There are other tools which do this, but specialdatetime.py has the advantage of being very lightweight.

### SATPLOTTOOLS ###

This is my plotting library. It builds on the functionality provided by matplotlib. It has methods for doing things such as polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and cartesian. 

### ASTRODYNAMICS2 ###

This is a library I made during a graduate-level astrodynamics course. It implements a few of the basics algorithms such as a very basic ECI-to-ECEF, as described in Vallado, Fundamentals for Astrodynamics and Applications. 

### Installation Instructions ###

* Clone the repository
* `python setup.py install`

Example code using omnireader:
```{python}
from geospacepy import omnireader
#Create a time window
sTimeIMF = datetime.datetime(2010,1,1)
eTimeIMF = datetime.datetime(2010,1,3)

#omni_interval is a dictionary-like object 
#that you can use to get the omni data for
#any variable as a numpy array 
#for any span of time
omniInt = omnireader.omni_interval(sTimeIMF,eTimeIMF,'5min')
t = omniInt['Epoch'] #datetime timestamps
By,Bz = omniInt['BY_GSM'],omniInt['BZ_GSM']

```
### Who do I talk to? ###

* This repository was created and is managed by Liam M. Kilcommons at CU Boulder
