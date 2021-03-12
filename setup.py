# Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# May, 2016
# (C) 2016 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group

import os
import glob
import textwrap

os.environ['DISTUTILS_DEBUG'] = "1"

from setuptools import setup, Extension
from setuptools.command import install as _install

setup(name='geospacepy',
      version = "0.2.2",
      description = "A small library for geospace data analysis",
      author = "Liam Kilcommons",
      author_email = 'liam.kilcommons@colorado.edu',
      url = "https://github.com/lkilcommons/geospacepy-lite",
      download_url = "https://github.com/lkilcommons/geospacepy-lite",
      long_description = textwrap.dedent("""Geospacepy-lite aims to be
                        a toolbox library for working with observations
                        from low-earth-orbiting spacecraft focusing on
                        complete vectorization of algorithms and 
                        consistant numpy array shape handling.
                        Submodules focus on transformations between various time
                        types, earth centered coordinate systems,
                        distance and area calculations using spherical geometry,
                        and solar-position-dependant calcuations."""),
      install_requires=['numpy','matplotlib'],
      packages=['geospacepy'],
      package_dir={'geospacepy' : 'geospacepy'},
      license='LICENSE.txt',
      zip_safe = False,
      classifiers = [
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Programming Language :: Python"
            ]
      )
