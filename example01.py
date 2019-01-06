#Example analyisis 

import numpy as np
import astropy.units as u
from uclahedp.sdfClass import sdfarr, sdfattrs

sname = r"C:\Users\Peter\Desktop\testsave.h5"
groupname = '/testrun/testprobe/'

groupname = '/'
    
z = np.ones([10,100,200,3])*u.V

a = [{'label':'x', 'axis':  np.arange(10)*u.cm},
      {'label':'y', 'axis': np.arange(100)*u.cm},
      {'label':'z', 'axis': np.arange(200)*u.cm},
      {'label':'rep', 'axis': np.arange(3)}]
   

#Create an array object
obj = sdfarr(data=z, axes=a)
   
#Convert some units 
print(obj.z.unit)
obj.convertAxisUnit('z', u.mm)
print(obj.z.unit)
obj.convertDataUnit('mV')


#Print out some information
print(obj.z.shape)
print(obj.dz)

#Perform some shape manipulations
print(obj.data.shape)
obj.avgDim('rep')
print(obj.data.shape)
obj.collapseDim('x', 3*u.cm)
print(obj.data.shape)
obj.thinDim('z', bin=5)
print(obj.data.shape)

#Save the data in an HDF file
obj.saveHDF(sname, groupname)
obj.close()



#Actually, let's re-open that arr object plot it
obj = sdfarr() #Initialize an empty object
obj.readHDF(sname, hdfpath=groupname)
print(obj.data.shape)

#Plot the object (It's now a 2D array)
obj.plot()

obj.close()


#Now lets add some attributes
a1 = {'one': 1, 'two' : [1,2,3] }
a2 = {'red': 'foo', 'blue' : 'bar' }

#Say we have two lists: maybe run and probe attributes
obj = sdfattrs(attrs=a1)
obj2 = sdfattrs(attrs=a2)

#Concatenate them into one
obj = obj + obj2

#Add a new attribute
obj.addAttr('newattr', 'Added!')

#Test some parameters
print(obj.containsKeys('one', 'red', 'blue'))

#Now save the data onto the same location as the arr data in the HDF
obj.saveHDF(sname, groupname)
obj.close()










    




