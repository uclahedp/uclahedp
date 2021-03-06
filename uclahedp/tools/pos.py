import numpy as np
import h5py
import os

from uclahedp.tools import csv as csvtools


def grid(pos, attrs, strict_axes=False, strict_grid=False, grid_precision=0.1,
         invert=False):
    """
    This program is a high-level encapsulation of all of the gridding typically
    used by one of the probe programs. 
    
    The default output of this function is: 
         shotgridind, xaxis, yaxis, zaxis, nx, ny, nz, nreps, nshots
         
    Where shotgridind is an array [nshots, 4]
    where each shot entry contains [xind, yind, zind, irep]. These indicies can
    then be used to place shot data in the correct grid locations.
    
    pos (array): Position array from the source file [nshots, 3] where the
    third dimension is xyz

    attrs (dict): Attribute dictionary for the dataset. Used to get probe
    origin, etc.
    
    strict_axes: If true, strict axis creation will be used.
    
    strict_grid: If true, strict gridding will be used. 
    
    grid_precision (float): If fuzzy gridding is used, values will be rounded
    to this precision before sorting into the grid.
    
    invert (bool): If true, returns an 'inverted' version of the shotgridind
    array which has size [nx, ny, nz, nreps], and each value in the array is
    the shot number that belongs there. This is useful when you need to load
    all reps at once, eg. if you need to average Langmuir traces prior to 
    the full analysis phase.
    
    
    
    """

    req_keys = []
        
    #Determine the motion format: default is cartesian
    if 'motion_format' in attrs.keys():
        motion_format = attrs['motion_format'][0]
    else:
        attrs['motion_format'] = 'cartesian'
      
    #Require any necessary keywords
    if motion_format == 'cartesian':
        pass
    elif motion_format == 'fixed_pivot':
        req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
        
        
        
    csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        
        
    #Set default values for keywords if they aren't present
    
    #origin is the location of the probe's origin in the main coordinates system
    origin = np.zeros(3)
    origin_keys = ['probe_origin_x','probe_origin_y','probe_origin_z']
    for i, key in enumerate(origin_keys):
         if key in attrs.keys():
              origin[i] = attrs[key][0]
         else:
              origin[i] = 0.0

    #probe_ax pol is the polarization of the probe's axes relative to the main
    #coordinates. 1 means they are aligned, -1 means anti-aligned.
    ax_pol = np.ones(3)
    ax_pol_keys = ['ax_pol_x','ax_pol_y','ax_pol_z']
    for i, key in enumerate(ax_pol_keys):
         if key in attrs.keys():
              ax_pol[i] = attrs[key][0]
         else:
              ax_pol[i] = 1
    

    #Adjust the probe for the origin and possible direction reversal
    pos[:,0] = pos[:,0]*ax_pol[0] + origin[0]
    pos[:,1] = pos[:,1]*ax_pol[1] + origin[1]
    pos[:,2] = pos[:,2]*ax_pol[2] + origin[2]
        

    #Generate the grid axes
    #TODO support non-cartesian axis gridding here? With a keyword or something?
    if strict_axes is True:
        try:
            print("Applying stric axis creation")
            nx = attrs['nx'][0]
            ny = attrs['ny'][0]
            nz = attrs['nz'][0]
            dx = attrs['dx'][0]
            dy = attrs['dy'][0]
            dz = attrs['dz'][0]
            x0 = attrs['x0'][0]
            y0 = attrs['y0'][0]
            z0 = attrs['z0'][0]
            xaxis, yaxis, zaxis = makeAxes(nx,ny,nz,
                                                    dx,dy,dz,
                                                    x0,y0,z0)
            
            #Apply these transformations to the made axes
            xaxis = xaxis*ax_pol[0] + origin[0]
            yaxis = yaxis*ax_pol[1] + origin[1]
            zaxis = zaxis*ax_pol[2] + origin[2]
            
        except KeyError:
            print("Missing axis parameters: attempting fuzzy axis creation")
            #Continue through to the strict_axes =False case if this happens
            strict_axes = False
                    
    if strict_axes is False:
        print("Applying fuzzy axis creation")
        xaxis,yaxis,zaxis = guessAxes(pos, precision=grid_precision)
        
                
    #Calculate length of axes
    nx, ny, nz, nreps = calcNpoints(pos, xaxis, yaxis, zaxis)
            

            
    #This line SHOULD be redundent, but sometimes it's necessary.
    #It is possible, when combining two motion lists, to get
    #some extra shots at the end of the datarun that don't have
    #positions associated. Recalculating this here ensures we only
    #take the ones relevant to this grid
    nshots = nx*ny*nz*nreps
            
    #If grid precision is zero, apply strict gridding
    #Otherwise, apply fuzzy gridding
    if strict_grid:
        print("Applying strict gridding")
        shotgridind = strictGrid(nx,ny,nz,nreps)  
    else:
        print("Applying fuzzy gridding")
        shotgridind = fuzzyGrid(pos, xaxis, yaxis, zaxis, 
                                         precision=grid_precision)
        
        
    if invert:
         shotgridind = invertShotGrid(shotgridind, xaxis, yaxis, zaxis)

    return shotgridind, xaxis, yaxis, zaxis, nx, ny, nz, nreps, nshots
        
        
        
        
            



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
    



def invertShotGrid(shotgridind, xaxis, yaxis, zaxis):
     """
     Takes a grid produces by one of the other functions and "inverts" it
     to give an array of the shot numbers at each position. This is useful
     if you are going to average over reps before saving full data (e.g.
     vsweep langmuir probe raw to full processing)
     """
     
     nx, ny, nz = len(xaxis), len(yaxis), len(zaxis)
     nshots = shotgridind.shape[0]
     npos = nx*ny*nz
     nreps = int(np.floor(nshots/npos))
     shots = np.zeros([nx, ny, nz, nreps], dtype=np.int32)
     
     for i in range(nshots):
          xi = shotgridind[i,0]
          yi = shotgridind[i,1]
          zi = shotgridind[i,2]
          ri = shotgridind[i,3]
          shots[xi,yi,zi,ri] = i
     return shots









     
    


if __name__ == "__main__":
 
    f=  os.path.join("F:", "LAPD_Mar2018", "RAW","run104_JanusBaO_raw.hdf5")
    with h5py.File(f, 'a') as sf:
          pos = sf['pos'][:]
        
        
    xaxis, yaxis, zaxis = guessAxes(pos)
    shotgridind = fuzzyGrid(pos, xaxis, yaxis, zaxis)
    shots = invertShotGrid(shotgridind, xaxis, yaxis, zaxis)
    
    
    #print(shotgridind[0:5, :])
    print(shots.shape)
    
    print(shots[0,0,0,:])
