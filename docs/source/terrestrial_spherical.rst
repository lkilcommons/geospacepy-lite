terrestrial_spherical
---------------------

This module contains various routines for manipulating vectors related to
positions and observations on and around the earth *under the approximation that
the earth is a sphere with radius 6371.2 meters*

The module assumes the following:
    - All latitudes are geo*centric* not geo*detic*
    - All altitudes are measured relative to the center of the earth, 
        not perpendicular to the earth's surface

