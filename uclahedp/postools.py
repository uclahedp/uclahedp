import numpy as np
import h5py
import hdftools



def shotClosestTo(pos, point=(0,0,0) ):
    dist = np.sqrt( (pos[:,0] - point[0])**2 + 
                  (pos[:,1] - point[1])**2 + 
                  (pos[:,2] - point[2])**2 )
    
    return np.argmin(dist)



def grid(pos, precision=.1):
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    nshots = len(xpos)
    
    xaxes, yaxes, zaxes = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    nx, ny, nz  = len(xaxes), len(yaxes), len(zaxes)
    nreps = int(np.floor(nshots/(nx*ny*nz)))
    gridind = np.zeros([nx, ny, nz, nreps],  dtype=np.int32)
    
    #This nested for loop isn't as bad as it looks, since the numbers are small
    #Still, there's probably a better way
    #...but what is it?
    for xi in range(nx):
        for yi in range(ny):
            for zi in range(nz):
                dist = np.sqrt( (xaxes[xi] - xpos)**2 + 
                      (yaxes[yi] - ypos)**2  +
                      (zaxes[zi] - zpos)**2 )
                si = np.array( np.where(dist <= precision), dtype=np.int32)
                gridind[xi,yi,zi,:] = si[0:nreps-1]
                
                
    return gridind, xaxes, yaxes, zaxes

if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_PL11B.hdf5") )
    
    with h5py.File(src.file, 'a') as sf:
        pos = sf[src.group]['pos'][:]
        
    point = (4.5,0,-9.5)
    print('Requested position: ' + str(point) )
    shot = shotClosestTo(pos, point=point ) 
    print('Closest shot: ' +  str(shot)  )
    print('Actual position of shot ' + str(np.round(pos[shot,:],1))   )
        
        
    
    gridind, xaxes, yaxes, zaxes = grid(pos)
    xi, yi, zi = 9, 0, 7
    predicted_pos = (xaxes[xi], yaxes[yi], zaxes[zi]  )
    print('Predicted position: ' + str(predicted_pos))
    shots = gridind[xi, yi, zi, :]
    print('Shot numbers: ' + str(shots))
    print('Actual positions at those shot numbers: \n' + str(pos[shots,:]) )