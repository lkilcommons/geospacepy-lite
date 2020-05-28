rotations
=========

This module provides a mechanism for simple coordinate transformations.
"Simple" here means coordinate tranformations which can be described 
as any number of serial rotations of a different coordinate frame 
about the x,y, or z axes.

There are many subtleties involved in describing coordinate rotation
mathematically (https://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities)

This module aims to emulate the convention used in the Vallado textbook:
    * Rotations are of the coordinate system "beneath" the vector
    * Rotation matrices ( M ) are to be "pre-multiplied" ( Mv ) by column vectors (v)

API
---

.. automodule:: geospacepy.rotations
    :members:

