#sdf class -> "standard data format"
#sdfattrs -> list of attributes

import h5py
from astropy import units as u
import copy

class sdfattrs:
    ##########################################################################
    # Basic Class Methods
    ##########################################################################
    def __init__(self, attrs={} ):
        """Initialize an sdfattrs object given keyword parameters"""
        self._setup(attrs)
        
    def close(self):
        pass
    
    def copy(self):
        """Return a deep copy of this sdfattrs object"""
        return copy.deepcopy(self)
        
    def _setup(self, attrs):
        """Setup the sdfarr, making some assumptions based on the keyword
        parameters given """
        self.attrs = attrs
        
    def _unpack(self, g):
        """Unpack the object from an HDF file group g"""
        self.attrs = {}
        keys = list(g.attrs)
        for key in keys:
            self.attrs[key] = g.attrs[key]
            
    def _pack(self, g):
        """Pack the object into an HDF file group g"""
        #Delete all the current attributes in the file
        for key in g.attrs.keys():
            del(g.attrs[key])

        for key in self.attrs.keys():
            g.attrs[key] =self.attrs[key]
            
    def readHDF(self, filepath, hdfpath = ''):
        with h5py.File(filepath, 'r') as f: 
            g = f.require_group(hdfpath)
            self._unpack(g)

    def saveHDF(self, filepath, hdfpath = ''):
        with h5py.File(filepath) as f: 
            g = f.require_group(hdfpath) 
            self._pack(g)
            
            
    def __getattr__(self, key):
        key = key.lower().strip()
        # If the axis is called, give the full astropy quantity object
        if key in self.attrs.keys():
            return self.attrs[key]

        #If you get this far, the key must be invalid
        print("Invalid Key: " + str(key))
        raise AttributeError
        
    def __add__(self, other):
        self.attrs.update(other.attrs)
        return self
        
        
    def addAttr (self, key, value):
        key = key.lower().strip()
        self.attrs[key] = value
        
    def delAttr(self, key):
        key = key.lower().strip()
        if key in self.attrs.keys():
            del self.attrs[key]
        else:
            raise Exception ("Key does not exist, so it can't be deleted!: " + str(key) )
            
    def containsKeys(self, *args):
        for key in args:
            if key.lower().strip() not in self.attrs.keys():
                return False
        return True


    ##########################################################################
    # Attribute list methods
    ##########################################################################
    
    
            
            
if __name__ == "__main__":
    sname = r"C:\Users\Peter\Desktop\TempData\testsave.h5"
    
    fname = r"C:\Users\Peter\Desktop\TempData\run56_LAPD1_pos_raw.h5"
    

    attrs = {'one': 1, 'two' : [1,2,3],
             'red':"string"}
    attrs2 = {'set':'sail', 'one':100}
    
    obj = sdfattrs(attrs=attrs)
    
    print(obj.containsKeys('one', 'two', 'orange'))
  

    obj.saveHDF(sname, '/test/')
    obj.close()