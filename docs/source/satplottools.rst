satplottools
============

Plotting functions originally created for plotting in-situ sensed data from
low-earth-orbit spacecraft (Defense Meterology Satellite Program).

.. note::

    This is one of the oldest parts of geospacepy-lite. There is a lot of 
    technical debt surrounding this code. It was written for matplotlib
    versions which did not yet have well-behaved polar plotting, so
    it bypasses the matplotlib projection=polar option, though with
    modern matplotlib, that would be the right solution.

.. automodule:: geospacepy.satplottools
    :members: