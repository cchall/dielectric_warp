"""
This file defines the class BoostedParticleDiagnostic

Major features:
- The class reuses the existing methods of FieldDiagnostic
  as much as possible, through class inheritance
- The class implements memory buffering of the slices, so as
  not to write to disk at every timestep

Remaining questions:
- Should one use the IO collectives when only a few proc modify a given file?
- Should we just have the proc writing directly to the file ?
  Should we gather on the first proc ?
- Is it better to write all the attributes of the openPMD file
  with only one proc ?
"""
import os
import numpy as np
import time
import h5py
from scipy.constants import c
from particle_diag import ParticleDiagnostic
from data_dict import z_offset_dict

class BoostedParticleDiagnostic(ParticleDiagnostic):
	def __init__(self, zmin_lab, zmax_lab, v_lab, dt_snapshots_lab,
                 Ntot_snapshots_lab, gamma_boost, period, em, top, w3d,
                 comm_world=None,  species = {"electrons": None},
                 particle_var=["position", "momentum", "weighting"],
                 select=None, write_dir=None, lparallel_output=False ):
	"""
	Initialize diagnostics that retrieve the data in the lab frame,
        as a series of snapshot (one file per snapshot),
        within a virtual moving window defined by zmin_lab, zmax_lab, v_lab.
                 
        Parameters
        ----------
        zmin_lab, zmax_lab: floats (meters)
            Positions of the minimum and maximum of the virtual moving window,
            *in the lab frame*, at t=0

        v_lab: float (m.s^-1)
            Speed of the moving window *in the lab frame*

        dt_snapshots_lab: float (seconds)
            Time interval *in the lab frame* between two successive snapshots

        Ntot_snapshots_lab: int
            Total number of snapshots that this diagnostic will produce

        period: int
            Number of iterations for which the data is accumulated in memory,
            before finally writing it to the disk. 
            
        See the documentation of ParticleDiagnostic for the other parameters

	"""

	if write_dir is None:
    		write_dir='lab_diags'

    	#initialize Particle diagnostic normal attributes
    	ParticleDiagnostic.__init__(period, top, w3d, comm_world,
                 species, particle_data, select, write_dir, lparallel_output)

    	# Register the boost quantities
        self.gamma_boost    = gamma_boost
        self.inv_gamma_boost= 1./gamma_boost
        self.beta_boost     = np.sqrt( 1. - self.inv_gamma_boost**2 )
        self.inv_beta_boost = 1./self.beta_boost

        # Find the z resolution and size of the diagnostic *in the lab frame*
        # (Needed to initialize metadata in the openPMD file)
        dz_lab              = c*self.top.dt * self.inv_beta_boost*self.inv_gamma_boost
        Nz                  = int( (zmax_lab - zmin_lab)/dz_lab )
        self.inv_dz_lab     = 1./dz_lab
        
        # Create the list of LabSnapshot objects
        self.snapshots      = []
        self.species        =species
        # Record the time it takes
        measured_start      = time.clock()
        print('\nInitializing the lab-frame diagnostics: %d files...' %(
            Ntot_snapshots_lab) )
        # Loop through the lab snapshots and create the corresponding files
        for i in range( Ntot_snapshots_lab ):
            t_lab   = i * dt_snapshots_lab
            snapshot= LabSnapshot( t_lab,
                                    zmin_lab + v_lab*t_lab,
                                    zmax_lab + v_lab*t_lab,
                                    self.write_dir, i )
            self.snapshots.append( snapshot )
            # Initialize a corresponding empty file
            self.create_file_empty_meshes( snapshot.filename, i,
                snapshot.t_lab, Nz, snapshot.zmin_lab, dz_lab, self.top.dt )
        # Print a message that records the time for initialization
        measured_end= time.clock()
        print('Time taken for initialization of the files: %.5f s' %(
            measured_end-measured_start) )

    def write(self):
        	"""
        	Redefines the method write of the parent class ParticleDiagnostic
        	Should be registered with installafterstep in Warp
        	"""
        	# At each timestep, store a slice of the particles in memory buffers
        	self.store_snapshot_slices()
        	# Every self.period, write the buffered slices to disk 
        	if self.top.it % self.period == 0:
        		self.flush_to_disk()

    def store_snapshot_slices( self ):
        """
        Store slices of the particles in the memory buffers of the
        corresponding lab snapshots
        """

        # Loop through the labsnapshots
        for snapshot in self.snapshots:

            # Update the positions of the output slice of this snapshot
            # in the lab and boosted frame (current_z_lab and current_z_boost)
            snapshot.update_current_output_positions( self.top.time,
                            self.inv_gamma_boost, self.inv_beta_boost )

            # For this snapshot:
            # - check if the output position *in the boosted frame*
            #   is in the current local domain
            # - check if the output position *in the lab frame*
            #   is within the lab-frame boundaries of the current snapshot
            self.ParticleCatcher= ParticleCatcher(
            self.gamma_boost, self.beta_boost, snapshot.current_z_boost,
            snapshot.prev_z_boost,snapshot.current_z_lab,
            snapshot.prev_z_lab)
            slice_array         = self.ParticleCatcher.extract_slice(self.species)
            snapshot.register_slice( slice_array, self.inv_dz_lab )
