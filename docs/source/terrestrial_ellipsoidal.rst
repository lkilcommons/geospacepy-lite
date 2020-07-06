terrestrial_ellipsoidal
=======================

This module contains routines which more accurately represent positions
on the earth than in terrestrial_spherical. Here the earth is modeled
as an ellipsoid.

Primarily this module deals with expression of positions in geodetic coordinates.

.. warning::
    
    The algorithms for dealing with positions on an ellipsoidal earth are
    less forgiving than for a spherical earth. Pay careful attention to
    details about inputs and outputs (e.g. units)

API
---

.. automodule:: geospacepy.terrestrial_ellipsoidal
    :members: