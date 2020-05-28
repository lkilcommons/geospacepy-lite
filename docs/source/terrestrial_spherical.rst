terrestrial_spherical
=====================

This module contains various routines for manipulating vectors related to
positions and observations on and around the earth *under the approximation that
the earth is a sphere with radius 6371.2 kilometers*

The spherical assumption implies the following:
    - All latitudes are geo*centric* not geo*detic*
    - All altitudes are measured relative to the center of the earth, not perpendicular to the earth's surface

Converting Multiple Vectors
---------------------------

If an input is described as a vector, it is expected to be a numpy array with the vector components as columns.

For example:

* One vector would be an array of shape ( 1 , 3 )
* *n* vectors would be an array of shape ( *n* , 3 )

Other inputs (not described as vectors) are broadcast to match the shape
of the vector input(s), if possible.

API
---

.. automodule:: geospacepy.terrestrial_spherical
    :members:
