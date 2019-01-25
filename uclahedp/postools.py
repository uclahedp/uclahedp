import numpy as np
import h5py
import hdftools



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


def gridShotIndGrid(pos, precision=.1):
    
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    xaxes, yaxes, zaxes = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    nshots = len(xpos)
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


def gridShotIndList(pos, precision=.1):
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    nshots = len(xpos)
    
    xaxes, yaxes, zaxes = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    nx, ny, nz  = len(xaxes), len(yaxes), len(zaxes)
    nreps = int(np.floor(nshots/(nx*ny*nz)))
    
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