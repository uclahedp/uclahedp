import numpy as np
import h5py

import os

from uclahedp import hdftools

def shotClosestTo(pos, point=(0,0,0) ):
    dist = np.sqrt( (pos[:,0] - point[0])**2 + 
                  (pos[:,1] - point[1])**2 + 
                  (pos[:,2] - point[2])**2 )
    
    return np.argmin(dist)



def makeAxes(pos, precision=.1):
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    xaxes, yaxes, zaxes = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    return xaxes, yaxes, zaxes


def gridShotIndList(pos, precision=.1):
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    nshots = len(xpos)
    
    xaxes, yaxes, zaxes = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    nx, ny, nz  = len(xaxes), len(yaxes), len(zaxes)
    
    nreps =np.floor(nshots/(nx*ny*nz))
    
    #This array will be used to count the number of reps we've found in each place
    #that will make sure each one gets a unique rep number
    countgrid = np.zeros( (nx,ny,nz) )
    
    #shotinds[i,:] = (xi, yi, zi, irep)
    shotind = np.zeros([nshots, 4],  dtype=np.int32)
    
    for i in range(nshots):
        xi = np.argmin( np.abs(xaxes - xpos[i]) )
        yi = np.argmin( np.abs(yaxes - ypos[i]) )
        zi = np.argmin( np.abs(zaxes - zpos[i]) )
        irep = countgrid[xi,yi,zi]
        
        # Only add up to nreps at each spot.
        # After that, shots will just be thrown away
        if irep < nreps:
            shotind[i,:] = np.array( (xi,yi,zi, irep )   )
            countgrid[xi,yi,zi] = countgrid[xi,yi,zi] + 1

    return shotind, xaxes, yaxes, zaxes



if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_PL11B.hdf5") )
    
    with h5py.File(src.file, 'a') as sf:
        pos = sf[src.group]['pos'][:]
        
        
        
    shotind, xaxes, yaxes, zaxes = gridShotIndList(pos)
        
    point = (4.5,0,-9.5)
    print('Requested position: ' + str(point) )
    shot = shotClosestTo(pos, point=point ) 
    print('Closest shot: ' +  str(shot)  )
    print('Actual position of shot ' + str(np.round(pos[shot,:],1))   )
        
        
    
    gridind, xaxes, yaxes, zaxes = gridShotIndGrid(pos)
    print(gridind.shape)
    xi, yi, zi = 9, 0, 7
    predicted_pos = (xaxes[xi], yaxes[yi], zaxes[zi]  )
    print('Predicted position: ' + str(predicted_pos))
    shots = gridind[xi, yi, zi, :]
    print('Shot numbers: ' + str(shots))
    print('Actual positions at those shot numbers: \n' + str(pos[shots,:]) )