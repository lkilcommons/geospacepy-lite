# Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# May, 2016
# (C) 2016 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group

import os
import glob

os.environ['DISTUTILS_DEBUG'] = "1"

from setuptools import setup, Extension
from setuptools.command import install as _install

setup(name='geospacepy',
      version = "0.1.0",
      description = "A collection of standalone modules for geospace data analysis",
      author = "Liam Kilcommons",
      author_email = 'liam.kilcommons@colorado.edu',
      url = "https://bitbucket.org/amienext/geospacepy-lite",
      download_url = "https://bitbucket.org/amienext/geospacepy-lite",
      long_description = "Just a toolbox of useful functions for geospace data analysis. Check individual modules for more info.",
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['spacepy','ephem','numpy','matplotlib','scipy'],
      packages=['geospacepy'],
      package_dir={'geospacepy' : 'geospacepy'},
      license='LICENSE.txt',
      zip_safe = False,
      classifiers = [
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "License :: MIT",
            "Natural Language :: English",
            "Programming Language :: Python"
            ]
      )