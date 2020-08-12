.. geospacepy-lite documentation master file, created by
   sphinx-quickstart on Fri Dec 13 12:00:49 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geospacepy-lite's documentation!
===========================================

Geospacepy-lite is a toolbox of loosely-related modules which were originally
created for analyzing in-situ sensed electrodynamics and particle precipitation data from spacecraft, but are general-purpose enough to useful for various common geospace data analysis tasks.

The dependacies of the package are limited to the core scientific python packages like the standard library, numpy, and matplotlib.

Conventions
-----------

- Angles are generally assumed to be in **radians** *except* :

  * Latitudes (assumed degrees)
  * Longitudes (assumed degrees)
  * Localtimes (when used as azimuth for a plot, i.e. instead of longitude)(assumed hours)

- ISO standards:

  * ISO 6709 (negative longitude means west longitude)
  * ISO 31-11 (spherical coordinate nomenclature: r, theta, phi, phi is azimuth)

- Times are assumed to be Universal Time (no timezones)

- Preferred time representation is julian date (use special_datetime to convert)

- SI units for physical quantities (e.g. earth radius is in meters)

Contribution
------------
    
Contributions are encouraged and welcomed. 

The Geospacepy-lite project maintainer pledges to abide by and expects contributors to also abide by accepted open source community kindness, professionalism and courtesy standards (e.g. as outlined in the  `contributor covenant <https://www.contributor-covenant.org/version/2/0/code_of_conduct//>`_ )

**Please make your contribution pull requests to the develop branch** 

**Please ensure contributions include a unit test**

Algorithms
----------

Most algorithms come from:

Vallado, D. A., & McClain, W. D. (2007). Fundamentals of Astrodynamics and Applications (3rd edition). Hawthorne, Calif.: Microcosm Press/Springer.

See the API documentation for additional algorithm sources.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   special_datetime
   sun
   rotations
   terrestrial_spherical
   terrestrial_ellipsoidal
   spherical_geometry
   satplottools

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
