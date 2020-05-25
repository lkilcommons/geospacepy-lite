.. geospacepy-lite documentation master file, created by
   sphinx-quickstart on Fri Dec 13 12:00:49 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geospacepy-lite's documentation!
===========================================

Geospacepy-lite is a toolbox of loosely-related modules which were originally
created for analyzing in-situ sensed electrodynamics and particle precipitation data from spacecraft, but are general-purpose enough to useful for various common geospace data analysis tasks.

.. note::
	
	Where multiple versions of an algorithm exist, this code generally implements the most approximate simple-minded version.

The dependacies of the package are limited as much as possible to the core scientific python packages like the standard library, numpy, and matplotlib.

Algorithms implemented in Geospacepy-lite come from many sources. I attempt to cite each algorithm within it's docstring. A majority are from:

Vallado, D. A., & McClain, W. D. (2007). Fundamentals of Astrodynamics and Applications (3rd edition). Hawthorne, Calif.: Microcosm Press/Springer.

Several conventions are used throughout geospacepy. These conventions are enforced in preference to adding, for instance, additional optional arguments.

- Angles are assumed radians *except*:

  * Latitudes (assumed degrees)
  * Longitudes (assumed degrees)
  * Localtimes (when used as azimuth, i.e. replacing Longitude)(assumed hours)

- Sign/Nomenclature conventions reflect ISO standards (when possible):

  * Longitudes are assumed to follow ISO 6709 (negative is west) convention

- Times are assumed to be approximate Universal Time (*leap seconds NOT handled*)
- Constants are specified in SI units (meters for lengths, kilograms for mass)

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   special_datetime
   sun
   rotations
   terrestrial_spherical
   spherical_geometry





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
