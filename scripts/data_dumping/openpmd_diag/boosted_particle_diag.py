"""
This file defines the class BoostedParticleDiagnostic

Major features:
- The class reuses the existing methods of ParticleDiagnostic
  as much as possible, through class inheritance
- The class implements memory buffering of the slices, so as
  not to write to disk at every timestep
"""
import os
import numpy as np
import time
from scipy.constants import c
from particle_diag import ParticleDiagnostic
from parallel import gatherarray

class BoostedParticleDiagnostic(ParticleDiagnostic):
	def __init__(self, zmin_lab, zmax_lab, v_lab, dt_snapshots_lab,
				 Ntot_snapshots_lab, gamma_boost, period, 
				 em, top, w3d, comm_world=None, 
				 particle_data=["position", "momentum", "weighting"],
				 select=None, write_dir=None, lparallel_output=False,
				 species = {"electrons": None}):
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

		# Do not leave write_dir as None, as this may conflict with
        # the default directory ('./diags') in which diagnostics in the
        # boosted frame are written
		if write_dir is None:
			write_dir='lab_diags'
		
		#initialize Particle diagnostic normal attributes
		ParticleDiagnostic.__init__(self, period, top, w3d, comm_world,
			species, particle_data, select, write_dir, lparallel_output)

		# Register the boost quantities
		self.gamma_boost = gamma_boost
		self.inv_gamma_boost = 1./gamma_boost
		self.beta_boost = np.sqrt(1. - self.inv_gamma_boost**2)
		self.inv_beta_boost = 1./self.beta_boost
		
		# Find the z resolution and size of the diagnostic *in the lab frame*
		# (Needed to initialize metadata in the openPMD file)
		dz_lab = c*self.top.dt*self.inv_beta_boost*self.inv_gamma_boost
		self.inv_dz_lab = 1./dz_lab
		
		# Create the list of LabSnapshot objects
		self.snapshots = []
		self.species = species

		# Record the time it takes
		measured_start = time.clock()
		print('\nInitializing the lab-frame diagnostics: %d files...' %(
			Ntot_snapshots_lab) )

		# Loop through the lab snapshots and create the corresponding files
		self.ParticleCatcher = ParticleCatcher(
			self.gamma_boost, self.beta_boost,top)
		self.ParticleCatcher.initialize_previous_instant()

		for i in range( Ntot_snapshots_lab ):
			t_lab = i*dt_snapshots_lab
			snapshot = LabSnapshot(	t_lab,
									zmin_lab + v_lab*t_lab,
									top.dt,
									zmax_lab + v_lab*t_lab,
									self.write_dir, i ,self.species)
			self.snapshots.append(snapshot)
			# Initialize a corresponding empty file
			self.create_file_empty_slice(
				snapshot.filename, i, snapshot.t_lab, self.top.dt)
			
		# Print a message that records the time for initialization
		measured_end = time.clock()
		print('Time taken for initialization of the files: %.5f s' %(
			measured_end-measured_start) )

	def write(self ):
		
		"""
		Redefines the method write of the parent class ParticleDiagnostic

		Should be registered with installafterstep in Warp
		"""
		# At each timestep, store a slice of the particles in memory buffers 
		self.store_snapshot_slices()
		# Every self.period, write the buffered slices to disk 
		#print snapshot.buffered_slices
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
			snapshot.update_current_output_positions(self.top.time,
							self.inv_gamma_boost, self.inv_beta_boost)
			
			# Setting up PartcleCatcher attributes with the updated snapshot 
			# attributes
			self.ParticleCatcher.zboost	= snapshot.current_z_boost
			self.ParticleCatcher.zboost_prev = snapshot.prev_z_boost
			self.ParticleCatcher.zlab = snapshot.current_z_lab
			self.ParticleCatcher.zlab_prev = snapshot.prev_z_lab
			self.ParticleCatcher.zmin_lab = snapshot.zmin_lab
			self.ParticleCatcher.zmax_lab = snapshot.zmax_lab

			for species_name in self.species_dict:
				species = self.species_dict[species_name]
				self.ParticleCatcher.species = species
				slice_array = self.ParticleCatcher.extract_slice(
					species, self.select)
				snapshot.register_slice(slice_array, species_name)

	def flush_to_disk(self):
		"""
		Writes the buffered slices of particles to the disk

		Erase the buffered slices of the LabSnapshot objects
		"""
		# Loop through the labsnapshots and flush the data
		
		for snapshot in self.snapshots:
			
			# Compact the successive slices that have been buffered
			# over time into a single array
			for species_name in self.species_dict:

				particle_array = snapshot.compact_slices(species_name)
				# Write this array to disk (if this snapshot has new slices)
				
				if particle_array.size:
					self.write_slices(particle_array, species_name, 
						snapshot, self.ParticleCatcher.particle_to_index)

			# Erase the memory buffers
			snapshot.buffer_initialization(self.species_dict)


	def write_boosted_dataset(self, species_grp, path, data, quantity):
		"""
		Writes each quantity of the buffered dataset to the disk, the 
		final step of the writing
		"""

		dset = species_grp[path]
		index = dset.shape[0]

		# Resize the h5py dataset 
		dset.resize(index + len(data), axis=0)

		# Write the data to the dataset at correct indices
		dset[index:] = data

	def write_slices( self, particle_array, species_name, snapshot, p2i ): 
		"""
		For one given snapshot, write the slices of the
		different fields to an openPMD file

		Parameters
		----------
		particle_array: array of reals
			Array of shape
			- (10, num_part, nslices) 

		iz_min, iz_max: integers
			The indices between which the slices will be written
			iz_min is inclusice and iz_max is exclusive

		snapshot: a LabSnaphot object

		p2i: dict
			Dictionary of correspondance between the particle quantities
			and the integer index in the particle_array
		"""
		# Open the file without parallel I/O in this implementation
		
		f = self.open_file(snapshot.filename)
		particle_path = "/data/%d/particles/%s" %(snapshot.iteration, 
			species_name)
		species_grp = f[particle_path]

		# Loop over the different quantities that should be written
		for particle_var in self.particle_data:
			# Scalar field
			if particle_var == "position":
				for coord in ["x","y","z"]:
					quantity= coord
					path = "%s/%s" %(particle_var, quantity)
					data = particle_array[ p2i[ quantity ] ]
					self.write_boosted_dataset(species_grp, path, data, 
						quantity)
				self.setup_openpmd_species_record(species_grp[particle_var], 
					particle_var)
	 
			elif particle_var == "momentum":
				for coord in ["x","y","z"]:
					quantity= "u%s" %coord
					path = "%s/%s" %(particle_var,coord)
					data = particle_array[ p2i[ quantity ] ]
					self.write_boosted_dataset( species_grp, path, data,
						quantity)
				self.setup_openpmd_species_record(species_grp[particle_var], 
					particle_var)
				
			elif particle_var == "weighting":
			   quantity= "w"
			   path	= 'weighting'
			   data	= particle_array[ p2i[ quantity ] ]
			   self.write_boosted_dataset(species_grp, path, data,
			   		quantity)
			   self.setup_openpmd_species_record(species_grp[particle_var], 
			   		particle_var)
			
		# Close the file
		f.close()


class LabSnapshot:
	"""
	Class that stores data relative to one given snapshot
	in the lab frame (i.e. one given *time* in the lab frame)
	"""
	def __init__(self, t_lab, zmin_lab, dt, zmax_lab, write_dir, i, 
		species_dict):
		"""
		Initialize a LabSnapshot 

		Parameters
		----------
		t_lab: float (seconds)
			Time of this snapshot *in the lab frame*
			
		zmin_lab, zmax_lab: floats
			Longitudinal limits of this snapshot

		write_dir: string
			Absolute path to the directory where the data for
			this snapshot is to be written

		i: int
		   	Number of the file where this snapshot is to be written

		species_dict: dict
			Contains all the species name of the species object 
			(inherited from Warp)
		"""

		# Deduce the name of the filename where this snapshot writes
		self.filename = os.path.join( write_dir, 'hdf5/data%08d.h5' %i)
		self.iteration = i
		self.dt = dt

		# Time and boundaries in the lab frame (constants quantities)
		self.zmin_lab = zmin_lab
		self.zmax_lab = zmax_lab
		self.t_lab = t_lab

		# Positions where the fields are to be registered
		# (Change at every iteration)
		self.current_z_lab = 0
		self.current_z_boost = 0

		self.buffer_initialization(species_dict)
	
	def buffer_initialization(self, species_dict):
		"""
		Initialize the buffer after each flush to disk

		Parameters
		----------
		species_dict: dict
			Contains all the species name of the species object 
			(inherited from Warp)
		"""

		self.buffered_slices = {}
		
		for species in species_dict:
			self.buffered_slices[species]=[]

	def update_current_output_positions( self, t_boost, inv_gamma, inv_beta ):
		"""
		Update the current and previous positions of output for this snapshot,
		so that if corresponds to the time t_boost in the boosted frame

		Parameters
		----------
		t_boost: float (seconds)
			Time of the current iteration, in the boosted frame

		inv_gamma, inv_beta: floats
			Inverse of the Lorentz factor of the boost, and inverse
			of the corresponding beta
		"""

		# Some shorcuts for further calculation's purposes
		t_lab = self.t_lab	
		t_boost_diff = t_boost + self.dt

		# This implements the Lorentz transformation formulas,
		# for a snapshot having a fixed t_lab
		self.current_z_boost = (t_lab*inv_gamma - t_boost)*c*inv_beta
		self.prev_z_boost = (t_lab*inv_gamma - t_boost_diff)*c*inv_beta		
		self.current_z_lab = (t_lab - t_boost*inv_gamma)*c*inv_beta
		self.prev_z_lab = (t_lab - t_boost_diff*inv_gamma)*c*inv_beta 

	def register_slice(self, slice_array, species):
		"""
		Store the slice of fields represented by slice_array
		and also store the z index at which this slice should be
		written in the final lab frame array

		Parameters
		----------
		slice_array: array of reals
			An array of packed fields that corresponds to one slice,
			as given by the SliceHandler object

		inv_dz_lab: float
			Inverse of the grid spacing in z, *in the lab frame*

		species: String, key of the species_dict
			Act as the key for the buffered_slices dictionary
		"""

		# Store the values and the index
		self.buffered_slices[species].append(slice_array)

	def compact_slices(self, species):
		"""
		Compact the successive slices that have been buffered
		over time into a single array, and return the indices
		at which this array should be written.

		Returns
		-------
		paticle_array: an array of reals of shape
		- (8, numPart) regardless of the dimension

		species: String, key of the species_dict
			Act as the key for the buffered_slices dictionary

		Returns None if the slices are empty
		"""
		
		particle_array = np.concatenate(self.buffered_slices[species], axis=1)
		
		return particle_array

class ParticleCatcher:
	"""
	Class that extracts, Lorentz-transforms and writes particles
	"""
	def __init__(self, gamma_boost, beta_boost, top):
		"""
		Initialize the ParticleHandler object

		Parameters
		----------
		gamma_boost, beta_boost: float
			The Lorentz factor of the boost and the corresponding beta

		top: WARP object
		"""

		# Some attributes neccessary for particle selections
		self.gamma_boost = gamma_boost
		self.beta_boost = beta_boost
		self.zboost = 0.0
		self.zboost_prev = 0.0
		self.zlab = 0.0
		self.zlab_prev = 0.0
		self.zlab_min = 0.0
		self.zlab_max = 0.0
		self.top = top
		self.species = None

		#self.top.allspecl = True
		
		# Create a dictionary that contains the correspondance
		# between the field names and array index
		self.particle_to_index = {'x':0, 'y':1, 'z':2, 'ux':3,
				'uy':4, 'uz':5, 'w':6, 'gamma':7}
			   
	def particle_getter(self):
		"""
		Returns the quantities of the particle at each time step

		Returns
		-------
		species : a warp Species object
		an array of quantities that are ready to be written in the buffer
		"""

		# the current quantities
		current_x =	self.get_quantity("x")
		current_y =	self.get_quantity("y")
		current_z =	self.get_quantity("z")
		current_ux = self.get_quantity("ux")
		current_uy = self.get_quantity("uy")
		current_uz = self.get_quantity("uz")
		current_weights = self.get_quantity( "w")

		# the previous quantities
		previous_x = self.get_quantity("x", l_prev = 1)
		previous_y = self.get_quantity("y", l_prev = 1)
		previous_z = self.get_quantity("z", l_prev = 1)
		previous_ux	= self.get_quantity("ux", l_prev = 1)
		previous_uy	= self.get_quantity("uy", l_prev = 1)
		previous_uz	= self.get_quantity("uz", l_prev = 1)
		
		# an array for mapping purpose
		z_array = np.arange(len(current_z))
			
		# For this snapshot:
		# - check if the output position *in the boosted frame*
		#   crosses the zboost in a forward motion
		# - check if the output position *in the lab frame*
		#   rosses the zboost_prev in a backward motion
		ii=np.compress((((current_z >= self.zboost) & 
			(previous_z <= self.zboost_prev)) |
		 	((current_z <= self.zboost) & 
		 	(previous_z >= self.zboost_prev))),z_array)
		
		## particle quantities that satisfy the aforementioned condition
		self.x_captured = np.take(current_x,ii)
		self.y_captured = np.take(current_y,ii)
		self.z_captured = np.take(current_z,ii)
		self.ux_captured = np.take(current_ux,ii)
		self.uy_captured = np.take(current_uy,ii)
		self.uz_captured = np.take(current_uz,ii)
		self.w_captured = np.take(current_weights,ii)
		self.gamma_captured = np.sqrt(1.+(self.ux_captured**2+\
			self.uy_captured**2+self.ux_captured**2)/c**2)

		self.x_prev_captured = np.take(previous_x,ii)
		self.y_prev_captured = np.take(previous_y,ii)
		self.z_prev_captured = np.take(previous_z,ii)
		self.ux_prev_captured = np.take(previous_ux,ii)
		self.uy_prev_captured = np.take(previous_uy,ii)
		self.uz_prev_captured = np.take(previous_uz,ii)
		self.gamma_prev_captured = np.sqrt(1.+(self.ux_prev_captured**2+\
			self.uy_prev_captured**2+self.ux_prev_captured**2)/c**2)

		num_part = np.shape(ii)[0]

		return num_part

	def transform_particles_to_lab_frame(self):
		uzfrm =- self.beta_boost*self.gamma_boost*c
		len_z = np.shape(self.z_captured)[0]
		
		self.top.setu_in_uzboosted_frame3d(len_z, 
				self.ux_captured,
				self.uy_captured, 
				self.uz_captured, 
				1./self.gamma_captured,
				uzfrm, self.gamma_boost)

		len_prev_z=np.shape(self.z_prev_captured)[0]

		self.top.setu_in_uzboosted_frame3d(len_prev_z, 
				self.ux_prev_captured,
				self.uy_prev_captured, 
				self.uz_prev_captured, 
				1./self.gamma_prev_captured,
				uzfrm, self.gamma_boost)
	
		self.z_captured	= (self.z_captured-self.zboost) + self.zlab 
		self.z_prev_captured = (self.z_prev_captured-self.zboost_prev) \
		+ self.zlab_prev
		self.t = self.gamma_boost*self.top.time*np.ones(len_z) \
		-uzfrm*self.z_captured/c**2
		self.t_prev	= self.gamma_boost*(self.top.time-self.top.dt)\
		*np.ones(len_prev_z) -uzfrm*self.z_prev_captured/c**2

	def collapse_to_mid_point(self):
		"""
		collapse the particle quantities to the mid point between 
		t_prev and t_current

		"""
		# Putting particles' current and previous time in an array for 
		# convenience in mean calculation
		time = [self.t, self.t_prev]
		t_mid = np.mean(time, axis=0)

		self.x_captured	 = self.x_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.x_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.y_captured	 = self.y_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.y_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.z_captured	 = self.z_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.z_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.ux_captured = self.ux_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.ux_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.uy_captured = self.uy_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.uy_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.uz_captured = self.uz_prev_captured*(self.t - t_mid)/\
		(self.t - self.t_prev) + self.uz_captured*(t_mid - self.t_prev)/\
		(self.t - self.t_prev)
		self.gamma_captured	 = self.gamma_prev_captured*(self.t - t_mid)/\
		(self.t-self.t_prev) + self.gamma_captured*(t_mid - self.t_prev)/\
		(self.t-self.t_prev)

	def gather_array(self, quantity):
		ar=np.zeros(np.shape(self.x_captured)[0])

		if quantity == "x" :
			ar=np.array(gatherarray(self.x_captured))
		elif quantity == "y":
			ar=np.array(gatherarray(self.y_captured))
		elif quantity == "z":
			ar=np.array(gatherarray(self.z_captured))
		elif quantity == "ux":
			ar=np.array(gatherarray(self.ux_captured*self.species.mass))
		elif quantity == "uy":
			ar=np.array(gatherarray(self.uy_captured*self.species.mass))
		elif quantity == "uz":
			ar=np.array(gatherarray(self.uz_captured*self.species.mass))
		elif quantity == "w":
			ar=np.array(gatherarray(self.w_captured))
		elif quantity == "gamma":
			ar=np.array(gatherarray(self.gamma_captured))
		return ar

	def extract_slice(self, species, select = None):
		"""
		Extract a slice of the particles at z_boost, using interpolation in z

		See the docstring of extract_slice for the parameters.

		Returns
		-------
		An array that packs together the slices of the different particles.
			The shape of this arrays is:
			 (8, num_part,)  

		"""
		# declare an attribute for convenience
		p2i = self.particle_to_index

		# Get the particles
		num_part = self.particle_getter()
		
		# Transform the particles from boosted frame back to lab frame
		self.transform_particles_to_lab_frame()
		
		# Collapse the particle quantities using interpolation to
		# the the midpoint of t and t_prev
		self.collapse_to_mid_point()

		slice_array = np.empty((8, num_part,))
                                                                                            
		for quantity in self.particle_to_index.keys():
			# Here typical values for 'quantity' are e.g. 'z', 'ux', 'gamma'
			slice_array[ p2i[quantity], ... ] = self.gather_array(quantity)

		# Choose the particles based on the select criteria defined by the 
		# users. Notice: this implementation still comes with a cost, 
		# one way to optimize it would be to do the selection before Lorentz
		# transformation back to the lab frame
		if (select is not None) and slice_array.size:
			select_array = self.apply_selection(select, slice_array)
			row, column =  np.where(select_array == True)
			temp_slice_array = slice_array[row,column]

			# Temp_slice_array is a 1D numpy array, we reshape it so that it 
			# has the same size as slice_array
			slice_array = np.reshape(temp_slice_array,(8,-1))

		return slice_array

	def get_quantity(self, quantity, l_prev = 0) :
		"""
		Get a given quantity

		Parameters
		----------
		species : a Species object
			Contains the species object to get the particle data from

		quantity : string
			Describes which quantity is queried
			Either "x", "y", "z", "ux", "uy", "uz" or "w"

		l_prev : boolean
			If 1, then return the quantities of the previous timestep;
			else return quantities of the current timestep
		"""
		# Extract the chosen quantities
		if l_prev:
			if quantity == "x" :
				quantity_array = self.species.getpid(id = self.top.xoldpid-1, 
					gather = 0, bcast = 0)
			elif quantity == "y" :
				quantity_array = self.species.getpid(id = self.top.yoldpid-1, 
					gather = 0, bcast = 0)
			elif quantity == "z" :
				quantity_array = self.species.getpid(id = self.top.zoldpid-1, 
					gather = 0, bcast = 0)
			elif quantity == "ux" :
				quantity_array = self.species.getpid(id = self.top.uxoldpid-1, 
					gather = 0, bcast = 0)
			elif quantity == "uy" :
				quantity_array = self.species.getpid(id = self.top.uyoldpid-1, 
					gather = 0, bcast = 0)
			elif quantity == "uz" :
				quantity_array = self.species.getpid(id = self.top.uzoldpid-1, 
					gather = 0, bcast = 0)
		else:
			if quantity == "x" :
				quantity_array = self.species.getx(gather = False)
			elif quantity == "y" :
				quantity_array = self.species.gety(gather = False)
			elif quantity == "z" :
				quantity_array = self.species.getz(gather = False)
			elif quantity == "ux" :
				quantity_array = self.species.getux(gather = False)
			elif quantity == "uy" :
				quantity_array = self.species.getuy(gather = False)
			elif quantity == "uz" :
				quantity_array = self.species.getuz(gather = False)
			elif quantity == "w" :
				quantity_array = self.species.getweights(gather = False)

		return quantity_array

	def initialize_previous_instant(self):
		"""
		Initialize the top.'quantity'oldpid array. This is used to store
		the previous values of the quantities.

		"""

		if not self.top.xoldpid:self.top.xoldpid = self.top.nextpid()
		if not self.top.yoldpid:self.top.yoldpid = self.top.nextpid()
		if not self.top.zoldpid:self.top.zoldpid = self.top.nextpid()
		if not self.top.uxoldpid:self.top.uxoldpid = self.top.nextpid()
		if not self.top.uyoldpid:self.top.uyoldpid = self.top.nextpid()
		if not self.top.uzoldpid:self.top.uzoldpid = self.top.nextpid()

	def apply_selection(self, select, slice_array) :
		"""
		Apply the rules of self.select to determine which
		particles should be written

		Parameters
		----------
		select : a dictionary that defines all selection rules based
		on the quantities

		Returns
		-------
		A 1d array of the same shape as that particle array
		containing True for the particles that satify all
		the rules of self.select
		"""

		p2i = self.particle_to_index

		# Initialize an array filled with True
		select_array = np.ones( np.shape(slice_array), dtype = 'bool' )

		# Apply the rules successively
		# Go through the quantities on which a rule applies
		for quantity in select.keys() :
			# Lower bound
			if select[quantity][0] is not None :
				select_array = np.logical_and(
					slice_array[p2i[quantity]] >\
					 select[quantity][0], select_array )
			# Upper bound
			if select[quantity][1] is not None :
				select_array = np.logical_and(
					slice_array[p2i[quantity]] <\
					select[quantity][1], select_array )

		return select_array 
