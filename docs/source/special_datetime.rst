special_datetime
================

This module converts between various time formats. It has two forms of each type of time conversion, a function which takes and returns only Python primatives (float,int) or datetime.datetime objects. These are the 'scalar' versions of each time conversion. There are also 'vector' implementations of each unit conversion, which take and return numpy arrays of the same types as their 'scalar' equivalents.

Time types handled (abbreviation: description) :

* ymdhms: year,month,day,hour,minute,second
* sod: second of day (0-86399) 
* doy: decimal day of year
* jd: julian date (Vallado, Algorithm 22) (jd)
* j2000: j2000 (julian date - julian date of 12:00 1-1-2000)
* datenum: matlab datenumber 

.. note::

	The code is only guaranteed to be able to convert between datetime and each of the above systems (it doesn't implement, for 
	example, julian date to MATLAB datenumber)

The naming convention for special_datetime functions is:

.. code-block:: python

	src_abbrev + ('arr' if vector else '') + '2' + dest_abbrev

For example:
	
* datetime2jd - converts from datetime.datetime to julian date
* sodarr2datetime - converts an array of seconds of the day to an array of datetimes

API
---

.. automodule:: geospacepy.special_datetime
    :members: