.. geospacepy-lite documentation master file, created by
   sphinx-quickstart on Fri Dec 13 12:00:49 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geospacepy-lite's documentation!
===========================================

Geospacepy-lite is a toolbox of loosely-related modules which were originally
created for analyzing in-situ sensed electrodynamics and particle precipitation data from spacecraft, but are general-purpose enough to useful for various common geospace data analysis tasks.

.. note::
	
	Where multiple versions of an algorithm exist, this code generally implements the most approximate simple-minded version. It's a place to start and get an approximate solution to a particular problem.

The dependacies of the package are limited as much as possible to the core scientific python packages like the standard library, numpy, and matplotlib.

Algorithms implemented in Geospacepy-lite come from many sources, and I attempt to cite each algorithm within it's docstring, but a majority of them are from:

Vallado, D. A., & McClain, W. D. (2007). Fundamentals of Astrodynamics and Applications (3rd edition). Hawthorne, Calif.: Microcosm Press/Springer.

Modules in Geospacepy-lite
++++++++++++++++++++++++++

special_datetime
----------------

This module converts between various time formats. It has two forms of each type of time conversion, a function which takes and returns only Python primatives (float,int) or datetime.datetime objects. These are the 'scalar' versions of each time conversion. There are also 'vector' implementations of each unit conversion, which take a return numpy arrays of the same types as their 'scalar' equivalents.

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



.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
