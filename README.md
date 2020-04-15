# Geospacepy #

[![Build Status](https://travis-ci.org/lkilcommons/geospacepy-lite.svg?branch=master)](https://travis-ci.org/lkilcommons/geospacepy-lite)

### What is this repository for? ###

Geospacepy is a small library of python functions for doing space science data analysis. It is also a dependency of several other applications.

### Rules of the Road ###

Geospacepy is MIT open source licensed software. The usual legal no-warranties, no-guarantees provisions apply. Though this software is indended for use in academic research it has not been peer-reviewed, there may be errors in calculations.  

### OMNIREADER ###

I envision that the most useful part of this library to others will be the omnireader.py code. It provides an functionality for on-demand downloading and reading of solar wind data from [NASA OmniWeb](http://omniweb.gsfc.nasa.gov/).

If you need fast OMNI data reading I recommend installing:

1. [The NASA CDF Library](http://cdf.gsfc.nasa.gov/)

2. [Spacepy](https://pypi.python.org/pypi/SpacePy), for it's pyCDF python interface to the NASA CDF Library

If geospacepy-lite does not detect SpacePy on your computer, it will use the NASA OMNIWeb ASCII text files.

### SPECIALDATETIME ###

This is a very basic set of utilities for handling numpy arrays of Python datetimes, and converting from and to other time formats, such as MATLAB datenumber, day-of-year, and Julian date. There are other tools which do this, but specialdatetime.py has the advantage of being very lightweight.

### SATPLOTTOOLS ###

This is my plotting library. It builds on the functionality provided by matplotlib. It has methods for doing things such as polar plots with latitude as the radial direction and local time as the azimuthal (so called 'dialplots') and convenience functions for conversion between latitude and longitude/local time and cartesian.

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
