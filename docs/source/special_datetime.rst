special_datetime
================

This module converts between Python datetime and other time represenations.

It has two forms of each time conversion:

* scalar : function takes and returns only Python primatives (float,int) or datetime.datetime objects
* vector : functions generated programatically from scalar functions which take numpy arrays or lists and return numpy arrays

Array shape convention for vector functions
-------------------------------------------

Vector functions output arrays which match the shape of the input.
If the inputs are lists or other iterables that don't have a shape, a flat
(1D) array is returned.

Time types handled (abbreviation: description)
----------------------------------------------

* ymdhms: year,month,day,hour,minute,second
* sod: decimal (float) second of day
* soy: decimal (float) second of year  
* doy: decimal (float) day of year
* jd: julian date (days since 12:00 PM on January 1, 4713 B.C.)
* j2000: j2000 (julian date, but relative to January 1, 2000 at 12:00 PM)
* datenum: matlab datenumber

Naming convention
-----------------

.. code-block:: python

	src_abbrev + ('arr' if vector else '') + '2' + dest_abbrev

For example:
	
* datetime2jd - converts from datetime.datetime to julian date
* sodarr2datetime - converts an array of seconds of the day to an array of datetimes

API
---

.. automodule:: geospacepy.special_datetime
    :members: