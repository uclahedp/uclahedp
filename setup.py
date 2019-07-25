#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: peter
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
          return f.read()

setup(name='uclahedp',
      version='0.1',
      description='Analysis programs written by the UCLA HEDP group.',
      url = 'https://github.com/phyzicist/uclahedp', 
      author = 'Peter Heuer',
      author_email = 'pheuer@physics.ucla.edu', 
      license = 'MIT',
      packages=['uclahedp'],
      install_requires=[
      	'h5py',
      	'bapsflib',
      	'pyqt5',
      	'matplotlib',
      	'numpy'
      ],
      zip_safe = False )
