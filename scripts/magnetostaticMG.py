"""Class for doing complete magnetostatic multigrid field solve"""
# ToDo:
#  - modify setj to check if particles are within grid
from warp import *
from fieldsolver import SubcycledPoissonSolver
from generateconductors import installconductors
from find_mgparam import find_mgparam

try:
  import psyco
except ImportError:
  pass

##############################################################################
class MagnetostaticMG(SubcycledPoissonSolver):
  
  __bfieldinputs__ = ['mgparam','downpasses','uppasses',
                      'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                      'lcndbndy','icndbndy','laddconductor',
                      'lcylindrical','lanalyticbtheta'] 
  __f3dinputs__ = ['gridmode','mgparam','downpasses','uppasses',
                   'mgmaxiters','mgtol','mgmaxlevels','mgform','mgverbose',
                   'lcndbndy','icndbndy','laddconductor','lprecalccoeffs'] 

  def __init__(self,**kw):
    self.grid_overlap = 2

    # --- Save input parameters
    self.processdefaultsfrompackage(MagnetostaticMG.__f3dinputs__,f3d,kw)
    self.processdefaultsfrompackage(MagnetostaticMG.__bfieldinputs__,
                                    f3d.bfield,kw)

    # --- Check for cylindrical geometry
    if self.lcylindrical: self.solvergeom = w3d.RZgeom

    SubcycledPoissonSolver.__init__(self,kwdict=kw)
    self.ncomponents = 3
    self.nxguard = 1
    self.nyguard = 1
    self.nzguard = 1
    self.lusevectorpotential = true

    # --- Kludge - make sure that the multigrid3df routines never sets up
    # --- any conductors.
    f3d.gridmode = 1

    # --- If there are any remaning keyword arguments, raise an error.
    assert len(kw.keys()) == 0,"Bad keyword arguemnts %s"%kw.keys()

    # --- Create a conductor object, which by default is empty.
    self.conductors = ConductorType()
    self.conductorlist = []
    self.newconductorlist = []

    # --- Give these variables dummy initial values.
    self.mgiters = zeros(3,'l')
    self.mgerror = zeros(3,'d')

    # --- Make sure that these are arrays
    self.mgmaxiters = ones(3)*self.mgmaxiters
    self.mgmaxlevels = ones(3)*self.mgmaxlevels
    self.mgparam = ones(3)*self.mgparam
    self.mgform = ones(3)*self.mgform
    self.mgtol = ones(3)*self.mgtol
    self.mgverbose = ones(3)*self.mgverbose
    self.downpasses = ones(3)*self.downpasses
    self.uppasses = ones(3)*self.uppasses

    # --- At the start, assume that there are no bends. This is corrected
    # --- in the solve method when there are bends.
    self.linbend = false

  def __getstate__(self):
    dict = SubcycledPoissonSolver.__getstate__(self)
    if self.lreducedpickle:
      del dict['conductors']
      dict['newconductorlist'] += self.conductorlist
      dict['conductorlist'] = []
    return dict

  def __setstate__(self,dict):
    SubcycledPoissonSolver.__setstate__(self,dict)
    if 'newconductorlist' not in self.__dict__:
      self.newconductorlist = self.conductorlist
      self.conductorlist = []
    if self.lreducedpickle and not self.lnorestoreonpickle:
      # --- Regenerate the conductor data
      self.conductors = ConductorType()

  def getconductorobject(self):
    for conductor in self.newconductorlist:
      self.installconductor(conductor)
    self.newconductorlist = []
    return self.conductors

  def getpdims(self):
    # --- Returns the dimensions of the jp, bp, and ap arrays
    return ((3,1+self.nxp,1+self.nyp,1+self.nzp),
            (3,1+self.nxp,1+self.nyp,1+self.nzp),
            (3,3+self.nxp,3+self.nyp,3+self.nzp))

  def getdims(self):
    # --- Returns the dimensions of the j, b, and a arrays
    return ((3,1+self.nxlocal,1+self.nylocal,1+self.nzlocal),
            (3,1+self.nxlocal,1+self.nylocal,1+self.nzlocal),
            (3,3+self.nxlocal,3+self.nylocal,3+self.nzlocal))

  def getj(self):
    'Returns the current density array'
    return self.source

  def getb(self):
    'Returns the B field array'
    return self.field
  
  def geta(self):
    'Returns the a array without the guard cells'
    ix1 = self.nxguard
    if ix1 == 0: ix1 = None
    ix2 = -self.nxguard
    if ix2 == 0: ix2 = None
    ix = slice(ix1,ix2)
    iy1 = self.nyguard
    if iy1 == 0: iy1 = None
    iy2 = -self.nyguard
    if iy2 == 0: iy2 = None
    iy = slice(iy1,iy2)
    iz1 = self.nzguard
    if iz1 == 0: iz1 = None
    iz2 = -self.nzguard
    if iz2 == 0: iz2 = None
    iz = slice(iz1,iz2)
    return self.potential[ix,iy,iz]

  def loadj(self,lzero=None,**kw):
    SubcycledPoissonSolver.loadsource(self,lzero,**kw)

  def fetchb(self,*args):
    SubcycledPoissonSolver.fetchfield(self,*args)

  def setsourcep(self,js,pgroup,zgrid):
    n  = pgroup.nps[js]
    if n == 0: return
    i  = pgroup.ins[js] - 1
    x  = pgroup.xp[i:i+n]
    y  = pgroup.yp[i:i+n]
    z  = pgroup.zp[i:i+n]
    ux = pgroup.uxp[i:i+n]
    uy = pgroup.uyp[i:i+n]
    uz = pgroup.uzp[i:i+n]
    gaminv = pgroup.gaminv[i:i+n]
    q  = pgroup.sq[js]
    w  = pgroup.sw[js]*top.pgroup.dtscale[js]
    if top.wpid > 0: wght = top.pgroup.pid[i:i+n,top.wpid-1]
    else:            wght = zeros((0,),'d')
    self.setsourcepatposition(x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,w)

  def setsourcepatposition(self,x,y,z,ux,uy,uz,gaminv,wght,zgrid,q,w):
    n = len(x)
    if n == 0: return
    if len(wght) > 0:
      nw = len(wght)
    else:
      nw = 0.
      wght = zeros(1,'d')
    setj3d(self.sourcep,self.sourcep,n,x,y,z,zgrid,ux,uy,uz,gaminv,
           q,w,nw,wght,top.depos,
           self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.lcylindrical)

  def fetchfieldfrompositions(self,x,y,z,ex,ey,ez,bx,by,bz,js=0,pgroup=None):
    n = len(x)
    if n == 0: return
    setb3d(self.fieldp,n,x,y,z,self.getzgridprv(),bx,by,bz,
           self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
           self.xmminp,self.ymminp,self.zmminp,
           self.l2symtry,self.l4symtry,self.lcylindrical)

  def fetchpotentialfrompositions(self,x,y,z,a):
    n = len(x)
    if n == 0: return
    fetchafrompositions3d(self.potentialp,n,x,y,z,a,self.getzgrid(),
                          self.nxp,self.nyp,self.nzp,self.dx,self.dy,self.dz,
                          self.xmminp,self.ymminp,self.zmminp,
                          self.l2symtry,self.l4symtry,self.lcylindrical)

  def setsourceforfieldsolve(self,*args):
    SubcycledPoissonSolver.setsourceforfieldsolve(self,*args)
    if self.lparallel:
      SubcycledPoissonSolver.setsourcepforparticles(self,*args)
      setjforfieldsolve3d(self.nxlocal,self.nylocal,self.nzlocal,self.source,
                          self.nxp,self.nyp,self.nzp,self.sourcep,
                          self.nxpguard,self.nypguard,self.nzpguard,
                          self.fsdecomp,self.ppdecomp)

  def getpotentialpforparticles(self,*args):
    """Despite the name, this actually gets the field instead, since that is
       always used in the magnetostatic solver"""
    if not self.lparallel:
      SubcycledPoissonSolver.getfieldpforparticles(self,*args)
    else:
      self.setfieldpforparticles(*args)
      getphipforparticles3d(3,self.nxlocal,self.nylocal,self.nzlocal,self.field,
                            self.nxp,self.nyp,self.nzp,self.fieldp,0,0,0,
                            self.fsdecomp,self.ppdecomp)


  def applysourceboundaryconditions(self):
    if ((self.pbounds[0] == 1 or self.l4symtry) and self.nx > 0 and
        self.solvergeom != w3d.RZgeom and
        self.fsdecomp.ix[self.fsdecomp.ixproc] == 0):
      self.source[:,0,:,:] = 2.*self.source[:,0,:,:]

    if (self.pbounds[1] == 1 and self.nx > 0 and
        self.fsdecomp.ix[self.fsdecomp.ixproc]+self.nxlocal == self.nx):
      self.source[:,-1,:,:] = 2.*self.source[:,-1,:,:]

    if ((self.pbounds[2] == 1 or self.l2symtry or self.l4symtry) and
        self.ny > 0 and self.fsdecomp.iy[self.fsdecomp.iyproc] == 0):
      self.source[:,:,0,:] = 2.*self.source[:,:,0,:]

    if (self.pbounds[3] == 1 and self.ny > 0 and
        self.fsdecomp.iy[self.fsdecomp.iyproc]+self.nylocal == self.ny):
      self.source[:,:,-1,:] = 2.*self.source[:,:,-1,:]

    if (self.pbounds[4] == 1 and self.nz > 0 and
        self.fsdecomp.iz[self.fsdecomp.izproc] == 0):
      self.source[:,:,:,0] = 2.*self.source[:,:,:,0]

    if (self.pbounds[5] == 1 and self.nz > 0 and
        self.fsdecomp.iz[self.fsdecomp.izproc]+self.nzlocal == self.nz):
      self.source[:,:,:,-1] = 2.*self.source[:,:,:,-1]

    if self.pbounds[0] == 2 or self.pbounds[1] == 2:
      if self.nxprocs == 1:
        self.source[:,0,:,:] = self.source[:,0,:,:] + self.source[:,-1,:,:]
        self.source[:,-1,:,:] = self.source[:,0,:,:]
      else:
        tag = 70
        if self.fsdecomp.ixproc == self.fsdecomp.nxprocs-1:
          ip = self.convertindextoproc(ix=self.fsdecomp.ixproc+1,
                                       bounds=self.pbounds)
          mpi.send(self.source[:,self.nxlocal,:,:],ip,tag)
          self.source[:,self.nxlocal,:,:],status = mpi.recv(ip,tag)
        elif self.fsdecomp.ixproc == 0:
          ip = self.convertindextoproc(ix=self.fsdecomp.ixproc-1,
                                       bounds=self.pbounds)
          sourcetemp,status = mpi.recv(ip,tag)
          self.source[:,0,:,:] = self.source[:,0,:,:] + sourcetemp
          mpi.send(self.source[:,0,:,:],ip,tag)

    if self.pbounds[2] == 2 or self.pbounds[3] == 2:
      if self.nyprocs == 1:
        self.source[:,:,0,:] = self.source[:,:,0,:] + self.source[:,:,-1,:]
        self.source[:,:,-1,:] = self.source[:,:,0,:]
      else:
        tag = 71
        if self.fsdecomp.iyproc == self.fsdecomp.nyprocs-1:
          ip = self.convertindextoproc(iy=self.fsdecomp.iyproc+1,
                                       bounds=self.pbounds)
          mpi.send(self.source[:,:,self.nylocal,:],ip,tag)
          self.source[:,:,self.nylocal,:],status = mpi.recv(ip,tag)
        elif self.fsdecomp.iyproc == 0:
          ip = self.convertindextoproc(iy=self.fsdecomp.iyproc-1,
                                       bounds=self.pbounds)
          sourcetemp,status = mpi.recv(ip,tag)
          self.source[:,:,0,:] = self.source[:,:,0,:] + sourcetemp
          mpi.send(self.source[:,:,0,:],ip,tag)

    if self.pbounds[4] == 2 or self.pbounds[5] == 2:
      if self.nzprocs == 1:
        self.source[:,:,:,0] = self.source[:,:,:,0] + self.source[:,:,:,-1]
        self.source[:,:,:,-1] = self.source[:,:,:,0]
      else:
        tag = 72
        if self.fsdecomp.izproc == self.fsdecomp.nzprocs-1:
          ip = self.convertindextoproc(iz=self.fsdecomp.izproc+1,
                                       bounds=self.pbounds)
          mpi.send(self.source[:,:,:,self.nzlocal],ip,tag)
          self.source[:,:,:,self.nzlocal],status = mpi.recv(ip,tag)
        elif self.fsdecomp.izproc == 0:
          ip = self.convertindextoproc(iz=self.fsdecomp.izproc-1,
                                       bounds=self.pbounds)
          sourcetemp,status = mpi.recv(ip,tag)
          self.source[:,:,:,0] = self.source[:,:,:,0] + sourcetemp
          mpi.send(self.source[:,:,:,0],ip,tag)

  def installconductor(self,conductor,
                            xmin=None,xmax=None,
                            ymin=None,ymax=None,
                            zmin=None,zmax=None,
                            dfill=top.largepos):
    if conductor in self.conductorlist: return
    self.conductorlist.append(conductor)
    installconductors(conductor,xmin,xmax,ymin,ymax,zmin,zmax,dfill,
                      self.getzgrid(),
                      self.nx,self.ny,self.nz,
                      self.nxlocal,self.nylocal,self.nzlocal,
                      self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                      self.zmmin,self.zmmax,1.,self.l2symtry,self.l4symtry,
                      solvergeom=self.solvergeom,
                      conductors=self.conductors,decomp=self.fsdecomp)

  def hasconductors(self):
    conductorobject = self.getconductorobject()
    return (conductorobject.interior.n > 0 or
            conductorobject.evensubgrid.n > 0 or
            conductorobject.oddsubgrid.n > 0)

  def clearconductors(self):
    self.conductors.interior.n = 0
    self.conductors.evensubgrid.n = 0
    self.conductors.oddsubgrid.n = 0

  def find_mgparam(self,lsavephi=false,resetpasses=1):
    find_mgparam(lsavephi=lsavephi,resetpasses=resetpasses,
                 solver=self,pkg3d=self)

  def dosolve(self,iwhich=0,*args):
    # --- Setup data for bends.
    rstar = fzeros(3+self.nzlocal,'d')
    if top.bends:
      setrstar(rstar,self.nzlocal,self.dz,self.zmminlocal,self.getzgrid())
      self.linbend = min(rstar) < largepos

    self.source[...] = self.source*mu0*eps0
    conductorobject = self.getconductorobject()

    if self.lcylindrical:
      init_bworkgrid(self.nxlocal,self.nzlocal,self.dx,self.dz,
                     self.xmminlocal,self.zmminlocal,self.bounds,
                     self.lparallel)

    # --- Note that the arrays being passed in are not contiguous, which means
    # --- that copies are being done.
    # --- If only initialization is being done (iwhich==1) then the bvp3d_work
    # --- routine only needs to be called once. Proper arrays are still passed
    # --- though they should never be needed during initialization.
    idmax = 2
    if iwhich == 1: idmax = 0
    for id in range(idmax+1):
      if (self.lanalyticbtheta and
          ((self.lusevectorpotential and (id == 0 or id == 2)) or
          (not self.lusevectorpotential and id == 1))): continue

      if self.lcylindrical:
        multigridrzb(iwhich,id,self.potential[id,:,1,:],
                     self.source[id,:,0,:],
                     self.nxlocal,self.nzlocal,self.mgtol[id])
      else:
        multigrid3dsolve(iwhich,self.nx,self.ny,self.nz,
                         self.nxlocal,self.nylocal,self.nzlocal,
                         self.dx,self.dy,self.dz,
                         self.potential[id,:,:,:],
                         self.source[id,:,:,:],
                         rstar,self.linbend,self.bounds,
                         self.xmmin,self.ymmin,self.zmmin,
                         self.mgparam[id],self.mgform[id],
                         self.mgiters[id],self.mgmaxiters[id],
                         self.mgmaxlevels[id],self.mgerror[id],
                         self.mgtol[id],self.mgverbose[id],
                         self.downpasses[id],self.uppasses[id],
                         self.lcndbndy,self.laddconductor,self.icndbndy,
                         self.gridmode,conductorobject,self.lprecalccoeffs,
                         self.fsdecomp)

  # # --- This is slightly inefficient in some cases, since for example, the
  # # --- MG solver already takes care of the longitudinal BC's.
  # setaboundaries3d(self.potential,self.nx,self.ny,self.nzlocal,
  #                  self.zmminlocal,self.zmmaxlocal,self.zmmin,self.zmmax,
  #                  self.bounds,self.lcylindrical,false)

    # --- Now take the curl of A to get B.
    getbfroma3d(self.potential,self.field,
                self.nxlocal,self.nylocal,self.nzlocal,
                self.dx,self.dy,self.dz,self.xmminlocal,
                self.lcylindrical,self.lusevectorpotential)

    # --- If using the analytic form of Btheta, calculate it here.
    if self.lanalyticbtheta:
      getanalyticbtheta(self.field,self.source,
                        self.nxlocal,self.nylocal,self.nzlocal,
                        self.dx,self.xmminlocal)

    # --- Unscale the current density
    self.source[...] = self.source/(mu0*eps0)


  ##########################################################################
  # Define the basic plot commands
  def genericpf(self,kw,pffunc):
    #kw['conductors'] = self.getconductorobject()
    kw['solver'] = self
    # --- This is a temporary kludge until the plot routines are updated to
    # --- use source and potential instead of rho and phi.
    self.j = self.source
    self.b = self.field
    self.a = self.potential
    pffunc(**kw)
    del self.j
    del self.b
    del self.a
  def pcjzy(self,**kw): self.genericpf(kw,pcjzy)
  def pcjzx(self,**kw): self.genericpf(kw,pcjzx)
  def pcjxy(self,**kw): self.genericpf(kw,pcjxy)
  def pcbzy(self,**kw): self.genericpf(kw,pcbzy)
  def pcbzx(self,**kw): self.genericpf(kw,pcbzx)
  def pcbxy(self,**kw): self.genericpf(kw,pcbxy)
  def pcazy(self,**kw): self.genericpf(kw,pcazy)
  def pcazx(self,**kw): self.genericpf(kw,pcazx)
  def pcaxy(self,**kw): self.genericpf(kw,pcaxy)

# --- This can only be done after MagnetostaticMG is defined.
try:
  psyco.bind(MagnetostaticMG)
except NameError:
  pass







