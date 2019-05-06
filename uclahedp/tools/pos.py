import numpy as np
import h5py

import os

from uclahedp.tools import hdf as hdftools

def shotClosestTo(pos, point=(0,0,0) ):
    dist = np.sqrt( (pos[:,0] - point[0])**2 + 
                  (pos[:,1] - point[1])**2 + 
                  (pos[:,2] - point[2])**2 )
    return np.argmin(dist)


def calcNpoints(pos, xaxis, yaxis, zaxis):
    nx,ny,nz = len(xaxis), len(yaxis), len(zaxis)
    nshots = len(pos)
    npos = nx*ny*nz
    nreps = int(np.floor(nshots/npos))
    
    return nx, ny, nz, nreps


def makeAxes(nx, ny, nz, dx, dy, dz, x0, y0, z0):
    """
    This function maxes axes a priori from grid specifications
    """
    #TODO: Do these formulas work if npoints is even??
    #We don't usually do this, but it's possible...
    xstep = x0 - dx*(nx - 1)/2.0
    ystep = y0 - dy*(ny - 1)/2.0
    zstep = z0 - dz*(nz - 1)/2.0
    
    xaxis = np.arange(nx)*dx + xstep
    yaxis = np.arange(ny)*dy + ystep
    zaxis = np.arange(nz)*dz + zstep
    
    #Sort the axes: even if dxyz was negative, we want the axis to be in
    #ascending order.
    xaxis.sort()
    yaxis.sort()
    zaxis.sort()
    
    return xaxis, yaxis, zaxis


def guessAxes(pos, precision=0.1):
    """
    This function guesses the axes for a given grid based on the position array
    """
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    xaxis, yaxis, zaxis = np.unique(xpos), np.unique(ypos), np.unique(zpos)
    
    return xaxis, yaxis, zaxis


def fuzzyGrid(pos, xaxis, yaxis, zaxis, precision=0.1):
    """
    Grid's positional data by trying to match each point to the nearest grid
    point. Nreps is calculated and enforced by throwing away any data that
    overfills a gridpoint.
    """
    #Round to within the precision given
    #+0.0 prevents weird -0's being in the resulting array
    xpos = np.round(pos[:,0]/precision, 0)*precision + 0.0
    ypos = np.round(pos[:,1]/precision, 0)*precision + 0.0
    zpos = np.round(pos[:,2]/precision, 0)*precision + 0.0
    
    
    nx, ny, nz, nreps = calcNpoints(pos, xaxis, yaxis, zaxis)
    
    nshots = int(nx*ny*nz*nreps)
    
    
    #This array will be used to count the number of reps we've found in each place
    #that will make sure each one gets a unique rep number
    countgrid = np.zeros((nx,ny,nz))
    

    #This array will hold the grid indices for each shot
    shotindarr = np.zeros([nshots, 4],  dtype=np.int32)
    
    for i in range(nshots):
        xi = np.argmin( np.abs(xaxis - xpos[i]) )
        yi = np.argmin( np.abs(yaxis - ypos[i]) )
        zi = np.argmin( np.abs(zaxis - zpos[i]) )
        irep = countgrid[xi,yi,zi]
        
        # Only add up to nreps at each spot.
        # After that, shots will just be thrown away
        if irep < nreps:
            shotindarr[i,:] = np.array( (xi,yi,zi, irep )   )
            countgrid[xi,yi,zi] = countgrid[xi,yi,zi] + 1

    return shotindarr




def strictGrid(nx,ny,nz,nreps):
    """
    Grid's positional data by assuming a number of reps and that the probe
    was moved in the standard X->Y->Z order. 
    """
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    nreps = int(nreps)
    nshots = int(nx*ny*nz*nreps)
    
    shotindarr = np.zeros([nshots, 4],  dtype=np.int32)
    
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for r in range(nreps):
                    shot = r + x*nreps + y*nreps*nx + z*nreps*nx*ny
                    shotindarr[shot, :] = [x,y,z,r]
    return shotindarr
    
    


if __name__ == "__main__":
 
    print(makeAxes(21,1,1,-1,0,0,0,0,0))
        
        
