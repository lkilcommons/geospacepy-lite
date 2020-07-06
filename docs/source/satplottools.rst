satplottools
============

Plotting functions originally created for plotting in-situ sensed data from
low-earth-orbit spacecraft (Defense Meterology Satellite Program).

.. warning::

    This is one of the oldest parts of geospacepy-lite. There is a lot of 
    technical debt surrounding this code. It does not use matplotlib's 
    projection=polar to make polar plots. It manually makes polar plots
    on a standard matplotlib axes.

.. autofunction:: geospacepy.satplottools.draw_dialplot

.. autofunction:: geospacepy.satplottools.latlt2cart

.. autofunction:: geospacepy.satplottools.latlon2cart

.. autofunction:: geospacepy.satplottools.multiline_timelabels

.. autofunction:: geospacepy.satplottools.vectorplot

.. autofunction:: geospacepy.satplottools.polarbinplot

