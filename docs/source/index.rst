.. geospacepy-lite documentation master file, created by
   sphinx-quickstart on Fri Dec 13 12:00:49 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geospacepy-lite's documentation!
===========================================

Geospacepy-lite is a toolbox of loosely-related modules which were originally
created for analyzing in-situ sensed electrodynamics and particle precipitation data from spacecraft, but are general-purpose enough to useful for various common geospace data analysis tasks.

The dependacies of the package are limited to the core scientific python packages like the standard library, numpy, and matplotlib.

Algorithms
----------

Most algorithms come from:

Vallado, D. A., & McClain, W. D. (2007). Fundamentals of Astrodynamics and Applications (3rd edition). Hawthorne, Calif.: Microcosm Press/Springer.

See the API documentation for additional algorithm sources.

Conventions
-----------

To keep geospacepy-lite maintainable, these conventions are enforced:

- Input or output arguments which are clearly angles will always be in radians *except* :

  * Latitudes (assumed degrees)
  * Longitudes (assumed degrees)
  * Localtimes (when used as azimuth, i.e. replacing Longitude)(assumed hours)

- ISO standards:

  * ISO 6709 (negative longitude means west longitude)

- Times are assumed to be Universal Time

- SI units for physical quantities (e.g. earth radius is in meters)

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
