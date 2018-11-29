#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ndfClassAgain.py: Another take on the ndf class

[Work in progress]

Created by Scott Feister on Wed Nov 28 20:26:47 2018
"""

# More general learning info on Python classes and special methods:
#  https://dbader.org/blog/python-dunder-methods
#  https://micropyramid.com/blog/python-special-class-methods-or-magic-methods/
#  http://getpython3.com/diveintopython3/special-method-names.html

# Destruction operators:
#  https://eli.thegreenplace.net/2009/06/12/safely-using-destructors-in-python/

import h5py

class ndf:
    """A simple UCLA Lab HDF5 file class. Works for either Raw or Full files."""
    
    def __init__(self, filename): # Initialize object
        self.filename = filename
        self.f = h5py.File(filename, 'r')
        
    def __del__(self): # Close object
        self.f.close()
    
    def load(self):
        """ Load the HDF5 file into RAM """
        # TODO
    
    def save(self, filename=None):
        """ Save the RAM object into an HDF5 file """
        if filename is None:
            filename = self.filename # Overwrite the current object
            self.f.close()
        
        with h5py.File(filename, 'w') as f:
            # TODO: Write the contents of the object into the HDF5 file
        
if __name__ == "__main__":
    pass
