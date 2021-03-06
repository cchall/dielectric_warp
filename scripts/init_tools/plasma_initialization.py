import numpy as np
from scipy.constants import c
from warp import PicklableFunction

class PlasmaInjector( object ):

    def __init__(self, elec, ions, w3d, top, dim, p_nx, p_ny, p_nz,
                 p_zmin, p_zmax, p_xmax, p_ymax, dens_func=None,
                 ux_m=0., uy_m=0., uz_m=0., ux_th=0., uy_th=0., uz_th=0.,
                 ncells_from_edge=2 ):
        """
        Initialize an injector for the plasma.
                 
        Parameters
        ----------
        elec: a Species object, as defined in Warp
           The particles that represent the electrons

        ions: a dictionary containing a list of ion species
            (one list per element in the dictionary; one ion species per
            charge state in each list) or None (no ions)
            Only the lowest charge state of each element is filled

        w3d, top: Forthon objects

        dim : str
           The dimension of the simulation (either "2d", "circ" or "3d")
        
        p_nx, p_ny, p_nz: int
           The number of particle per cell along each direction

        p_zmin, p_zmax, p_ymax, p_zmax: floats
           The positions between which the electrons should be initialized
           (The electrons stop being injected after p_zmax, even for
           continuous injection. p_zmax is given at t=0.)

        dens_func: callable
           A function of the form dens_func( x, y, z )
           where x, y, z are 1darrays of macroparticles positions
           (arrays of the same length) and which returns a 1darray of the
           same length containing the relative density

        ux_m, uy_m, uz_m: floats (dimensionless)
           Normalized mean momenta in each direction

        ux_th, uy_th, uz_th: floats (dimensionless)
           Normalized thermal momenta in each direction
           Only the electrons have these momenta added to their mean momentum.

        ncells_from_edge: int
           The number of cells between the initial injection position and
           the right end of the simulation box.
           For periodic simulations, it is advised to use ncells_from_edge=0
           so that the plasma fills the whole.
           For simulations with moving window and continuous injection, it is
           advised to keep the default (ncells_from_edge=2) so that the injector
           never leaves the plasma
        """
        # Register the species to be injected        
        self.elec = elec
        self.ions = ions
        self.dim = dim
        if dens_func is not None:
            # --- The PicklableFunction allows the dens_func reference to be pickled
            self.dens_func = PicklableFunction(dens_func)
        else:
            self.dens_func = None
        self.w3d = w3d
        self.top = top
        self.p_nz = p_nz

        # Momenta parameters
        self.ux_m = ux_m
        self.uy_m = uy_m
        self.uz_m = uz_m
        self.gamma_inv_m = 1./np.sqrt( 1 + ux_m**2 + uy_m**2 + uz_m**2 )
        self.ux_th = ux_th
        self.uy_th = uy_th
        self.uz_th = uz_th
        
        # Get the 1d arrays of evenly-spaced positions for the particles
        # - Positions along x, in one given slice
        dx = w3d.dx / p_nx  # Spacing between particles
        nx_local = (w3d.xmmaxlocal - w3d.xmminlocal) / dx
        x_reg = w3d.xmminlocal + dx*( np.arange( nx_local ) + 0.5 )
        self.x_reg = x_reg[ (x_reg>=-p_xmax) & (x_reg<=p_xmax) ]
        # - Positions along y, in the 3d case
        if dim == "3d":
            dy = w3d.dy / p_ny
            ny_local = (w3d.ymmaxlocal - w3d.ymminlocal) / dy
            y_reg = w3d.ymminlocal + dy*( np.arange( ny_local ) + 0.5 )
            self.y_reg = y_reg[ (y_reg>=-p_ymax) & (y_reg<=p_ymax) ]
        # - Angular positions, in the circ case
        elif dim == "circ":
            dtheta = 2*np.pi / p_ny
            self.theta_reg = dtheta * np.arange( p_ny )

        # For the continuous injection:
        zmmax =  w3d.zmmax + top.zgrid
        # Continuously varying injection position (moves with moving window)
        self.z_inject = zmmax - ncells_from_edge*w3d.dz
        # Maximum position after which the plasma stops being injected
        # (moves with the plasma)
        self.p_zmax = p_zmax
        # Position of the right end of the plasma (moves with the plasma)
        self.z_end_plasma = min( p_zmax, zmmax - ncells_from_edge*w3d.dz )
        # Mean speed of the end of the plasma
        self.v_plasma = c * uz_m * self.gamma_inv_m
            
        # Inject particles between p_zmin and self.z_end_plasma
        # (within the limits of the domain)
        self.load_plasma( self.z_end_plasma,
                          max(p_zmin, w3d.zmminlocal + top.zgrid),
                          min(p_zmax, w3d.zmmaxlocal + top.zgrid) )

    def load_plasma( self, z_end_plasma, zmin, zmax ):
        """
        Load plasma between zmin and zmax.
        
        The positions of the particles along the z axis are of the form
        z = z_end_plasma - i*dz - 0.5*dz
        and satisfy zmin <= z < zmax
        where i is an integer, and dz is the spacing
        between particles

        Parameters
        ----------
        z_end_plasma : float
           Position of the global end of the plasma
        
        zmin, zmax : floats (meters)
           Positions between which the plasma is to be loaded, in the
           local domain
        """
        # Get 1d array of evenly-spaced positions for the particles along z
        dz = self.w3d.dz / self.p_nz
        # Get the min and max indices i for z = z_end_plasma - i*dz - 0.5*dz
        i_min = int( (z_end_plasma-zmin)/dz - 0.5 ) 
        i_max = int( (z_end_plasma-zmax)/dz + 0.5 )
        i_max = max( i_max, 0 )
        z_reg = z_end_plasma - dz*( np.arange( i_max, i_min+1 ) + 0.5 )

        # Get the corresponding particle positions at injection
        if self.dim == "3d":
            zp, xp, yp = np.meshgrid( z_reg, self.x_reg, self.y_reg )
            z0 = zp.flatten()
            x0 = xp.flatten()
            y0 = yp.flatten()
            r0 = np.sqrt( x0**2 + y0**2 )
        elif self.dim == "circ":
            # (copy=True is important here, since it allows to
            # change the theta angles individually)
            zp, rp, thetap = np.meshgrid( z_reg, self.x_reg,
                                          self.theta_reg, copy=True )
            # Prevent the particles from being aligned along any direction
            theta = unalign_angles( thetap )
            r0 = rp.flatten()
            z0 = zp.flatten()
            x0 = r0 * np.cos( theta )
            y0 = r0 * np.sin( theta )
        elif self.dim == "2d":
            zp, xp = np.meshgrid( z_reg, self.x_reg )
            z0 = zp.flatten()
            x0 = xp.flatten()
            y0 = np.zeros_like(x0)
            r0 = abs(x0)
        elif self.dim == "1d":
            z0 = z_reg
            x0 = y0 = np.zeros_like(z0)
            r0 = abs(x0)

        # Modulate the weights according to the density
        # (Take into account the motion of the plasma to retrieve the
        # the position where this slice of plasma was at t=0, this
        # is because the dens_func is given at t=0)
        if self.dens_func is not None:
            w = self.dens_func( x0, y0, z0 - self.v_plasma*self.top.time )
        else:
            w = np.ones_like( z0 - self.v_plasma*self.top.time )
        # In circ, the particles have larger weight at higher radius
        if self.dim == "circ":
            w = w * r0 / self.w3d.xmmax

        # Add random momenta for the particles
        ux = self.ux_m
        uy = self.uy_m
        uz = self.uz_m
        gamma_inv = self.gamma_inv_m
        n_part = len(z0)
        if self.ux_th != 0:
            ux += self.ux_th * np.random.normal( size=n_part )
        if self.uy_th != 0:
            uy += self.uy_th * np.random.normal( size=n_part )
        if self.uz_th != 0:
            uz += self.uz_th * np.random.normal( size=n_part )
        if (self.ux_th !=0) or (self.uy_th !=0) or (self.uy_th !=0):
            gamma_inv = 1./np.sqrt( 1 + ux**2 + uy**2 + uz**2 )

        # Load the electrons and ions (on top of each other)
        if self.elec is not None:
            # Use the random momenta
            self.elec.addpart( x=x0, y=y0, z=z0,
                    vx=c*ux*gamma_inv, vy=c*uy*gamma_inv,
                    vz=c*uz*gamma_inv, gi=gamma_inv, w=w,
                    lallindomain=False )
        if self.ions is not None:
            # For each element, only add particles to the lowest charge state
            for element in self.ions.keys():
                # Use only the mean momenta
                lowest_state_species = self.ions[ element ][0]
                lowest_state_species.addpart( x=x0, y=y0, z=z0,
                    vx=c*self.ux_m*self.gamma_inv_m,
                    vy=c*self.uy_m*self.gamma_inv_m,
                    vz=c*self.uz_m*self.gamma_inv_m,
                    gi=self.gamma_inv_m, w=w,
                    lallindomain=False )
                
    def continuous_injection(self):
        """
        Routine which is called by warp at each timestep
        """
        # Move the injection position with the moving window
        self.z_inject += self.top.vbeamfrm * self.top.dt
        # Move the position of the end of the plasma by its mean velocity
        self.z_end_plasma += self.v_plasma * self.top.dt
        # Move the position of the limit beyond which no plasma is injected
        # (It moves along with the plasma because the user gives p_zmax at t=0)
        self.p_zmax += self.v_plasma * self.top.dt

        # Add slices filled with plasma
        while (self.z_end_plasma < self.z_inject) and \
            (self.z_end_plasma < self.p_zmax) :

            # Add one slice
            self.load_plasma( self.z_end_plasma + self.w3d.dz,
                              self.z_end_plasma,
                              self.z_end_plasma + self.w3d.dz )

            # One slice has been added ; increment the position of plasma end
            self.z_end_plasma += self.w3d.dz
                


def unalign_angles( thetap ) :
    """
    Shifts the angles by a random amount

    Parameter
    ---------
    thetap: 3darray
       An array of shape Nr, Nz, Ntheta, where Nr, Nz, Ntheta are the
       number of regularly-spaced macroparticles along each direction.
       The angles evenly sample [0,2 \pi], when changing the last index

    Returns
    -------
    A 1darray of length Nr*Nz*Ntheta, where a random angle has been added to
    each set of Ntheta macroparticles that are at a given position in r and z
    """
    # Generate random angles between 0 and 2 pi
    # (but only for the two first axis, to the last)
    theta_shift = 2*np.pi*np.random.rand( thetap.shape[0], thetap.shape[1] )

    # For each r and z position, add the same random shift for all Ntheta values
    theta = thetap + theta_shift[:,:,np.newaxis]
    
    # Flatten
    return( theta.flatten() )







    


    
