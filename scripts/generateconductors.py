"""
This module contains classes for generating the conductor data from a
combination of simple geometrical elements.
The following elements are defined:

Box(xsize,ysize,zsize,...)
Cylinder(radius,length,theta=0.,phi=0.,...)
ZCylinder(radius,length,...)
ZCylinderOut(radius,length,...)
ZRoundedCylinderOut(radius,length,radius2,...)
YCylinder(radius,length,...)
XCylinder(radius,length,...)
Sphere(radius,...)
Cone(r_zmin,r_zmax,length,theta=0.,phi=0.,...)
ZCone(r_zmin,r_zmax,length,...)
ZConeOut(r_zmin,r_zmax,length,...)
ConeSlope(slope,intercept,length,theta=0.,phi=0.,...)
ZConeSlope(slope,intercept,length,...)
ZConeOutSlope(slope,intercept,length,...)
ZTorus(r1,r2,...)
Beamletplate(za,zb,z0,thickness,...)

Note that all take the following additional arguments:
voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1

installconductors(a): generates the data needed for the fieldsolve
"""

# The following classes are also defined but should not be directly used.
# Grid
# Assembly
# AssemblyNot
# AssemblyAnd
# AssemblyPlus
# AssemblyMinus
# Delta

from warp import *
if not lparallel: import VPythonobjects

generateconductorsversion = "$Id: generateconductors.py,v 1.20 2003/04/17 18:22:42 dave Exp $"
def generateconductors_doc():
  import generateconductors
  print generateconductors.__doc__

##############################################################################
def installconductors(a,xmin=None,xmax=None,ymin=None,ymax=None,
                        zmin=None,zmax=None,dfill=top.largepos):
  """
Installs the given conductors.
  - a: the assembly of conductors
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These can be set for optimization, to avoid looking
    for conductors where there are none. Also, they can be used crop a
    conductor
  - dfill=largepos: points at a depth in the conductor greater than dfill
                    are skipped.
  """
  # First, create a grid object
  g = Grid(xmin,xmax,ymin,ymax,zmin,zmax)
  # Generate the conductor data
  g.getdata(a,dfill)
  # Then install it
  g.installdata()

##############################################################################
##############################################################################
##############################################################################
class Assembly:
  """
Class to hold assemblies of conductors.  Base class of all conductors.
Should never be directly created by the user.
 - v=0.: voltage on conductor
 - x,y,z=0.,0.,0: center of conductor
 - condid=1: conductor identification number
 - kwlist=[]: list of string names of variable describing conductor
 - generatorf=None: function which generates the distances between the points
                    and the conductors along the axis.
 - generatord=None: function which generates the smallest distance between the
                    points and the conductor surface.
  """

  voltage = 0.
  xcent = 0.
  ycent = 0.
  zcent = 0.

  def __init__(self,v=0.,x=0.,y=0.,z=0.,condid=1,kwlist=[],
                    generatorf=None,generatord=None):
    self.voltage = v
    self.xcent = x
    self.ycent = y
    self.zcent = z
    self.condid = condid
    self.kwlist = kwlist
    self.generatorf = generatorf
    self.generatord = generatord

  def getkwlist(self):
    kwlist = []
    for k in self.kwlist:
      kwlist.append(self.__dict__[k])
    kwlist.append(self.__dict__['xcent'])
    kwlist.append(self.__dict__['ycent'])
    kwlist.append(self.__dict__['zcent'])
    return kwlist

  def griddistance(self,ix,iy,iz,xx,yy,zz):
    result = Delta(ix,iy,iz,xx,yy,zz,voltage=self.voltage,condid=self.condid,
                   generator=self.generatorf,kwlist=self.getkwlist())
    return result

  def distance(self,xx,yy,zz):
    result = Distance(xx,yy,zz,generator=self.generatord,
                      kwlist=self.getkwlist())
    return result

  def isinside(self,xx,yy,zz):
    result = IsInside(xx,yy,zz,generator=self.generatord,
                      condid=self.condid,kwlist=self.getkwlist())
    return result

  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    pass  # virtual function

  # Operations which return an Assembly expression.
  def __mul__(self,right):
    return AssemblyAnd(self,right)
  def __add__(self,right):
    return AssemblyPlus(self,right)
  def __sub__(self,right):
    return AssemblyMinus(self,right)
  def __neg__(self):
    return AssemblyNot(self)


class AssemblyNot(Assembly):
  """
AssemblyNot class.  Represents 'not' of assemblies.
  """
  def __init__(self,l):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (-(self.left.griddistance(ix,iy,iz,xx,yy,zz)))
  def distance(self,xx,yy,zz):
    return (-(self.left.distance(xx,yy,zz)))
  def isinside(self,xx,yy,zz):
    return (-(self.left.isinside(xx,yy,zz)))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)


class AssemblyAnd(Assembly):
  """
AssemblyAnd class.  Represents 'and' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) *
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) *
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) *
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)


class AssemblyPlus(Assembly):
  """
AssemblyPlus class.  Represents 'or' of assemblies.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) +
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) +
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) +
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)


class AssemblyMinus(Assembly):
  """
AssemblyMinus class.
  """
  def __init__(self,l,r):
    Assembly.__init__(self,0.,l.xcent,l.ycent,l.zcent,l.condid)
    self.left = l
    self.right = r
  def griddistance(self,ix,iy,iz,xx,yy,zz):
    return (self.left.griddistance(ix,iy,iz,xx,yy,zz) -
            self.right.griddistance(ix,iy,iz,xx,yy,zz))
  def distance(self,xx,yy,zz):
    return (self.left.distance(xx,yy,zz) -
            self.right.distance(xx,yy,zz))
  def isinside(self,xx,yy,zz):
    return (self.left.isinside(xx,yy,zz) -
            self.right.isinside(xx,yy,zz))
  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    self.left.visualize(xmin,xmax,ymin,ymax)
    self.right.visualize(xmin,xmax,ymin,ymax)

##############################################################################

class Delta:
  """
Class to hold the set of distances in each of the six directions.
Distances have the sign of the outward normal surface vector, i.e.
distances to outside the surface are positive, inside negative.
  """

  def __init__(self,ix=None,iy=None,iz=None,xx=None,yy=None,zz=None,
                    dels=None,vs=None,ns=None,
                    parity=None,voltage=0.,condid=1,generator=None,kwlist=[]):
    if ix is None:
      self.ndata = 0
      nn = 10000
      self.ix = zeros(nn)
      self.iy = zeros(nn)
      self.iz = zeros(nn)
      self.xx = zeros(nn,'d')
      self.yy = zeros(nn,'d')
      self.zz = zeros(nn,'d')
      self.dels = zeros((6,nn),'d')
      self.vs = zeros((6,nn),'d')
      self.ns = zeros((6,nn))
      self.parity = zeros(nn)
      self.mglevel = zeros(nn)
    elif generator is not None:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = zeros((6,self.ndata),'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.dels[0,:],self.dels[1,:],
                                self.dels[2,:],self.dels[3,:],
                                self.dels[4,:],self.dels[5,:]] + [fuzz])
      self.setvoltages(voltage)
      self.setcondids(condid)
      self.setlevels(0)
    else:
      self.ndata = len(ix)
      self.ix = ix
      self.iy = iy
      self.iz = iz
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.dels = dels
      self.vs = vs
      self.ns = int(ns)
      self.parity = parity
      self.setlevels(0)
   
  def setvoltages(self,voltage):
    "Routine to set appropriate voltages."
    self.vs = voltage + zeros((6,self.ndata),'d')
   
  def setcondids(self,condid):
    "Routine to setcondid condids."
    self.ns = int(condid) + zeros((6,self.ndata))
   
  def setlevels(self,level):
    self.mglevel = level + zeros(self.ndata)

  def normalize(self,dx,dy,dz):
    """
Normalizes the data with respect to the grid cell sizes.
dx,dy,dz: the grid cell sizes
    """
    self.dels[:,:] = self.dels/array([dx,dx,dy,dy,dz,dz])[:,NewAxis]

  def setparity(self,dfill):
    """
Set parity. For points inside, this is set to -1. For points near the surface,
this is set to the parity of ix+iy+iz. Otherwise defaults to large integer.
This assumes that the data has already been normalized with respect to the
grid cell sizes.
    """
    self.parity = zeros(self.ndata) + 999
    fuzz = 1.e-9
    # --- A compiled routine is called for optimization
    setconductorparity(self.ndata,self.ix,self.iy,self.iz,
                       self.dels,self.parity,fuzz,f3d.fuzzsign,dfill)

  def clean(self):
    """
Removes the data which is far from any conductors. Assumes that setparity
has already been called.
    """
    ii = compress(self.parity < 2,arange(self.ndata))
    self.ix    = take(self.ix,ii)
    self.iy    = take(self.iy,ii)
    self.iz    = take(self.iz,ii)
    self.xx    = take(self.xx,ii)
    self.yy    = take(self.yy,ii)
    self.zz    = take(self.zz,ii)
    self.dels  = take(self.dels,ii,1)
    self.vs    = take(self.vs,ii,1)
    self.ns    = take(self.ns,ii,1)
    self.parity= take(self.parity,ii)
    self.ndata = len(self.ix)

  def append(self,d):
    n1 = self.ndata
    n2 = d.ndata
    if n1 + n2 > len(self.ix):
      ix = self.ix[:n1]
      iy = self.iy[:n1]
      iz = self.iz[:n1]
      xx = self.xx[:n1]
      yy = self.yy[:n1]
      zz = self.zz[:n1]
      dels = self.dels[:,:n1]
      vs = self.vs[:,:n1]
      ns = self.ns[:,:n1]
      parity = self.parity[:n1]
      mglevel = self.mglevel[:n1]

      newn = max(int(2*len(self.ix)),n1+n2)
      self.ix = zeros(newn)
      self.iy = zeros(newn)
      self.iz = zeros(newn)
      self.xx = zeros(newn,'d')
      self.yy = zeros(newn,'d')
      self.zz = zeros(newn,'d')
      self.dels = zeros((6,newn),'d')
      self.vs = zeros((6,newn),'d')
      self.ns = zeros((6,newn))
      self.parity = zeros(newn)
      self.mglevel = zeros(newn)

      self.ix[:n1] = ix
      self.iy[:n1] = iy
      self.iz[:n1] = iz
      self.xx[:n1] = xx
      self.yy[:n1] = yy
      self.zz[:n1] = zz
      self.dels[:,:n1] = dels
      self.vs[:,:n1] = vs
      self.ns[:,:n1] = ns
      self.parity[:n1] = parity
      self.mglevel[:n1] = mglevel

    self.ix[n1:n1+n2] = d.ix[:n2]
    self.iy[n1:n1+n2] = d.iy[:n2]
    self.iz[n1:n1+n2] = d.iz[:n2]
    self.xx[n1:n1+n2] = d.xx[:n2]
    self.yy[n1:n1+n2] = d.yy[:n2]
    self.zz[n1:n1+n2] = d.zz[:n2]
    self.dels[:,n1:n1+n2] = d.dels[:,:n2]
    self.vs[:,n1:n1+n2] = d.vs[:,:n2]
    self.ns[:,n1:n1+n2] = d.ns[:,:n2]
    self.parity[n1:n1+n2] = d.parity[:n2]
    self.mglevel[n1:n1+n2] = d.mglevel[:n2]
    self.ndata = n1 + n2

  def install(self):
    """
Installs the data into the WARP database
    """
    ntot = 0
    nc = f3d.ncond
    nn = sum(where(self.parity[:self.ndata] == -1,1,0))
    ntot = ntot + nn
    if nn > 0:
      if nc + nn > f3d.ncondmax:
        f3d.ncondmax = nn + nc
        gchange("Conductor3d")
      f3d.ncond = f3d.ncond + nn
      ii = compress(self.parity[:self.ndata] == -1,arange(self.ndata))
      f3d.ixcond[nc:nc+nn] = take(self.ix,ii)
      f3d.iycond[nc:nc+nn] = take(self.iy,ii)
      f3d.izcond[nc:nc+nn] = take(self.iz,ii)
      f3d.condvolt[nc:nc+nn] = take(self.vs[0,:],ii)
      f3d.condnumb[nc:nc+nn] = take(self.ns[0,:],ii)
      f3d.icondlevel[nc:nc+nn] = take(self.mglevel,ii)

    ne = f3d.necndbdy
    nn = sum(where(self.parity[:self.ndata] == 0,1,0))
    ntot = ntot + nn
    if nn > 0:
      if ne + nn > f3d.ncndmax:
        f3d.ncndmax = nn + ne
        gchange("Conductor3d")
      f3d.necndbdy = f3d.necndbdy + nn
      ii = compress(self.parity[:self.ndata] == 0,arange(self.ndata))
      f3d.iecndx[ne:ne+nn] = take(self.ix,ii)
      f3d.iecndy[ne:ne+nn] = take(self.iy,ii)
      f3d.iecndz[ne:ne+nn] = take(self.iz,ii)
      f3d.ecdelmx[ne:ne+nn] = take(self.dels[0,:],ii)
      f3d.ecdelpx[ne:ne+nn] = take(self.dels[1,:],ii)
      f3d.ecdelmy[ne:ne+nn] = take(self.dels[2,:],ii)
      f3d.ecdelpy[ne:ne+nn] = take(self.dels[3,:],ii)
      f3d.ecdelmz[ne:ne+nn] = take(self.dels[4,:],ii)
      f3d.ecdelpz[ne:ne+nn] = take(self.dels[5,:],ii)
      f3d.ecvoltmx[ne:ne+nn] = take(self.vs[0,:],ii)
      f3d.ecvoltpx[ne:ne+nn] = take(self.vs[1,:],ii)
      f3d.ecvoltmy[ne:ne+nn] = take(self.vs[2,:],ii)
      f3d.ecvoltpy[ne:ne+nn] = take(self.vs[3,:],ii)
      f3d.ecvoltmz[ne:ne+nn] = take(self.vs[4,:],ii)
      f3d.ecvoltpz[ne:ne+nn] = take(self.vs[5,:],ii)
      f3d.ecnumbmx[ne:ne+nn] = take(self.ns[0,:],ii)
      f3d.ecnumbpx[ne:ne+nn] = take(self.ns[1,:],ii)
      f3d.ecnumbmy[ne:ne+nn] = take(self.ns[2,:],ii)
      f3d.ecnumbpy[ne:ne+nn] = take(self.ns[3,:],ii)
      f3d.ecnumbmz[ne:ne+nn] = take(self.ns[4,:],ii)
      f3d.ecnumbpz[ne:ne+nn] = take(self.ns[5,:],ii)
      f3d.ecvolt[ne:ne+nn] = take(self.vs[0,:],ii)
      f3d.ecnumb[ne:ne+nn] = take(self.ns[0,:],ii)
      f3d.iecndlevel[ne:ne+nn] = take(self.mglevel,ii)

    no = f3d.nocndbdy
    nn = sum(where(self.parity[:self.ndata] == 1,1,0))
    ntot = ntot + nn
    if nn > 0:
      if no + nn > f3d.ncndmax:
        f3d.ncndmax = nn + no
        gchange("Conductor3d")
      f3d.nocndbdy = f3d.nocndbdy + nn
      ii = compress(self.parity[:self.ndata] == 1,arange(self.ndata))
      f3d.iocndx[no:no+nn] = take(self.ix,ii)
      f3d.iocndy[no:no+nn] = take(self.iy,ii)
      f3d.iocndz[no:no+nn] = take(self.iz,ii)
      f3d.ocdelmx[no:no+nn] = take(self.dels[0,:],ii)
      f3d.ocdelpx[no:no+nn] = take(self.dels[1,:],ii)
      f3d.ocdelmy[no:no+nn] = take(self.dels[2,:],ii)
      f3d.ocdelpy[no:no+nn] = take(self.dels[3,:],ii)
      f3d.ocdelmz[no:no+nn] = take(self.dels[4,:],ii)
      f3d.ocdelpz[no:no+nn] = take(self.dels[5,:],ii)
      f3d.ocvoltmx[no:no+nn] = take(self.vs[0,:],ii)
      f3d.ocvoltpx[no:no+nn] = take(self.vs[1,:],ii)
      f3d.ocvoltmy[no:no+nn] = take(self.vs[2,:],ii)
      f3d.ocvoltpy[no:no+nn] = take(self.vs[3,:],ii)
      f3d.ocvoltmz[no:no+nn] = take(self.vs[4,:],ii)
      f3d.ocvoltpz[no:no+nn] = take(self.vs[5,:],ii)
      f3d.ocnumbmx[no:no+nn] = take(self.ns[0,:],ii)
      f3d.ocnumbpx[no:no+nn] = take(self.ns[1,:],ii)
      f3d.ocnumbmy[no:no+nn] = take(self.ns[2,:],ii)
      f3d.ocnumbpy[no:no+nn] = take(self.ns[3,:],ii)
      f3d.ocnumbmz[no:no+nn] = take(self.ns[4,:],ii)
      f3d.ocnumbpz[no:no+nn] = take(self.ns[5,:],ii)
      f3d.ocvolt[no:no+nn] = take(self.vs[0,:],ii)
      f3d.ocnumb[no:no+nn] = take(self.ns[0,:],ii)
      f3d.iocndlevel[no:no+nn] = take(self.mglevel,ii)
    if ntot > 0:
      if(w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom):
        frz.install_conductors_rz()

  def __neg__(self):
    "Delta not operator."
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 -self.dels,self.vs,self.ns)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(self.dels,right.dels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,right.dels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    rdels = -right.dels
    c = less(self.dels,rdels)
    return Delta(self.ix,self.iy,self.iz,self.xx,self.yy,self.zz,
                 choose(c,(self.dels,rdels)),
                 choose(c,(self.vs  ,right.vs)),
                 choose(c,(self.ns  ,right.ns)))

  def __str__(self):
    "Prints out delta"
    return repr(self.dels)+" "+repr(self.vs)+" "+repr(self.ns)

##############################################################################

class Distance:
  """
Class to hold the distance between points and a conductor
Distances have the sign of the outward normal surface vector, i.e.
distances outside the surface are positive, inside negative.
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    distance=None,generator=None,kwlist=[]):
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = zeros(self.ndata,'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                self.distance[:]])
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.distance = distance
   
  def __neg__(self):
    "Delta not operator."
    return Distance(self.xx,self.yy,self.zz, -self.distance)

  def __mul__(self,right):
    "'and' operator, returns maximum of distances to surfaces."
    c = less(abs(self.distance),abs(right.distance))
    return Distance(self.xx,self.yy,self.zz,
                    choose(c,(self.distance,right.distance)))

  def __add__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    c = greater(abs(self.distance),abs(right.distance))
    dd = Distance(self.xx,self.yy,self.zz,
                  choose(c,(self.distance,right.distance)))
    dd.distance = where((self.distance < 0.) & (right.distance > 0.),
                        self.distance,dd.distance)
    dd.distance = where((self.distance > 0.) & (right.distance < 0.),
                        right.distance,dd.distance)
    return dd

  def __sub__(self,right):
    "'or' operator, returns minimum of distances to surfaces."
    # --- Warning - while this work for many cases, there is no
    # --- gaurantee of robustness! It should only work when the right
    # --- hand side is a cylinder.
    rdistance = -right.distance
    dd = Distance(self.xx,self.yy,self.zz,self.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance >= 0.),
                        sqrt(rdistance**2+self.distance**2),
                        dd.distance)
    dd.distance = where((rdistance >= 0.) & (self.distance <= 0.),
                        rdistance,dd.distance)
    dd.distance = where((rdistance < 0.) & (self.distance <= 0.),
                        maximum(rdistance,self.distance),dd.distance)
    return dd


  def __str__(self):
    "Prints out delta"
    return repr(self.distance)

##############################################################################

class IsInside:
  """
Class to hold flag whether or not a point is inside a conductor.
  """

  def __init__(self,xx=None,yy=None,zz=None,
                    isinside=None,generator=None,condid=1,kwlist=[]):
    self.condid = condid
    if generator is not None:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      distance = zeros(self.ndata,'d')
      fuzz = 1.e-13
      apply(generator,kwlist + [self.ndata,self.xx,self.yy,self.zz,
                                distance[:]])
      self.isinside = where(distance <= 0.,condid,0.)
    else:
      self.ndata = len(xx)
      self.xx = xx
      self.yy = yy
      self.zz = zz
      self.isinside = isinside*self.condid
   
  def __neg__(self):
    "Delta not operator."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_not(self.isinside),condid=self.condid)

  def __mul__(self,right):
    "'and' operator, returns logical and of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,right.isinside),
                    condid=self.condid)

  def __add__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_or(self.isinside,right.isinside),
                    condid=self.condid)

  def __sub__(self,right):
    "'or' operator, returns logical or of isinsides."
    return IsInside(self.xx,self.yy,self.zz,
                    logical_and(self.isinside,logical_not(right.isinside)),
                    condid=self.condid)

  def __str__(self):
    "Prints out delta"
    return repr(self.isinside)

##############################################################################
##############################################################################
##############################################################################

class Grid:
  """
Class holding the grid info.
Constructor arguments:
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of conductors. Defaults to the
    mesh size. These only need to be set for optimization, to avoid looking
    for conductors where there are none. They can also be used to crop a
    conductor.
Call getdata(a,dfill) to generate the conductor data. 'a' is a geometry object.
Call installdata() to install the data into the WARP database.
  """

  def __init__(self,xmin=None,xmax=None,ymin=None,ymax=None,
                    zmin=None,zmax=None):
    """
Creates a grid object which can generate conductor data.
    """
    self.xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    self.xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    self.ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    self.ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    self.zmin = where(zmin==None,w3d.zmmin,zmin)[0]
    self.zmax = where(zmax==None,w3d.zmmax,zmax)[0]

    self.nx = w3d.nx
    self.ny = w3d.ny
    self.nz = w3d.nz
    self.nzfull = w3d.nzfull
    self.xmmin = w3d.xmmin
    self.ymmin = w3d.ymmin
    self.zmmin = w3d.zmminglobal

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    self.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
    if(w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom):
      self.dy = self.dx
    else:
      self.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    self.dz = (w3d.zmmax - w3d.zmmin)/w3d.nz

    if top.fstype in [7,11,10]:
      if top.fstype in [7,11]:
        setmglevels(self.nx,self.ny,self.nz,self.nzfull,self.dx,self.dy,self.dz)
      if top.fstype == 10:
        init_base(self.nx,self.nz,self.dx,self.dz,0.,self.zmin,false)
      self.mglevels = f3d.mglevels
      self.mglevelsnx = f3d.mglevelsnx[:f3d.mglevels]
      self.mglevelsny = f3d.mglevelsny[:f3d.mglevels]
      self.mglevelsiz = f3d.mglevelsiz[:f3d.mglevels]
      self.mglevelsnz = f3d.mglevelsnz[:f3d.mglevels]
      self.mglevelslx = f3d.mglevelslx[:f3d.mglevels]
      self.mglevelsly = f3d.mglevelsly[:f3d.mglevels]
      self.mglevelslz = f3d.mglevelslz[:f3d.mglevels]
    else:
      self.mglevels = 1
      self.mglevelsnx = [self.nx]
      self.mglevelsny = [self.ny]
      self.mglevelsiz = [top.izfsslave[me]]
      self.mglevelsnz = [self.nz]
      self.mglevelslx = [1]
      self.mglevelsly = [1]
      self.mglevelslz = [1]

  def getmesh(self,mglevel=0):
    i = mglevel
    dx = self.dx*self.mglevelslx[i]
    dy = self.dy*self.mglevelsly[i]
    dz = self.dz*self.mglevelslz[i]
    nx = self.mglevelsnx[i]
    ny = self.mglevelsny[i]
    iz = self.mglevelsiz[i]
    nz = self.mglevelsnz[i]

    zmmin = self.zmmin + iz*dz

    xmesh = self.xmmin + dx*arange(nx+1)
    ymesh = self.ymmin + dy*arange(ny+1)
    zmesh =      zmmin + dz*arange(nz+1)
    xmesh = compress(logical_and(self.xmin-dx <= xmesh,
                                                 xmesh <= self.xmax+dx),xmesh)
    ymesh = compress(logical_and(self.ymin-dy <= ymesh,
                                                 ymesh <= self.ymax+dy),ymesh)
    zmesh = compress(logical_and(self.zmin-dz <= zmesh,
                                                 zmesh <= self.zmax+dz),zmesh)
    x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
    y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
    z = zeros(len(xmesh)*len(ymesh),'d')
    ix = nint((x - self.xmmin)/dx)
    iy = nint((y - self.ymmin)/dy)
    iz = zeros(len(xmesh)*len(ymesh))
    return ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh

  def getdata(self,a,dfill=top.largepos):
    """
Given an Assembly, accumulate the appropriate data to represent that
Assembly on this grid.
 - a: the assembly
 - dfill=top.largepos: points at a depth in the conductor greater than dfill
                       are skipped.
    """
    starttime = wtime()
    tt2 = zeros(8,'d')
    self.dall = Delta()
    for i in range(self.mglevels):
      tt1 = wtime()
      ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(i)

      tt2[0] = tt2[0] + wtime() - tt1
      if len(x) == 0: continue
      for zz in zmesh:
        tt1 = wtime()
        z[:] = zz
        iz[:] = nint((zz - zmmin)/dz)
        tt2[1] = tt2[1] + wtime() - tt1
        tt1 = wtime()
        d = a.griddistance(ix,iy,iz,x,y,z)
        tt2[2] = tt2[2] + wtime() - tt1
        tt1 = wtime()
        d.normalize(dx,dy,dz)
        tt2[3] = tt2[3] + wtime() - tt1
        tt1 = wtime()
        d.setparity(dfill)
        tt2[4] = tt2[4] + wtime() - tt1
        tt1 = wtime()
        d.clean()
        tt2[5] = tt2[5] + wtime() - tt1
        tt1 = wtime()
        d.setlevels(i)
        tt2[6] = tt2[6] + wtime() - tt1
        tt1 = wtime()
        self.dall.append(d)
        tt2[7] = tt2[7] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

  def installdata(self):
    """
Installs the conductor data into the fortran database
    """
    self.dall.install()

  def getdistances(self,a,mglevel=0):
    """
Given an Assembly, accumulate the distances between the assembly and the
grid points.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(mglevel)
    try:
      self.distances[0,0,0]
    except:
      self.distances = zeros((1+nx,1+ny,1+nz),'d')
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz
      iz[:] = nint((zz - zmmin)/dz)
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.distance(x,y,z)
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.distance
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.distances[ix1:ix2+1,iy1:iy2+1,iz[0]] = dd
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

  def getisinside(self,a,mglevel=0):
    """
Given an Assembly, set flag for each grid point whether it is inside the
assembly.
 - a: the assembly
 - mglevel=0: coarsening level to use
    """
    starttime = wtime()
    tt2 = zeros(4,'d')
    tt1 = wtime()
    ix,iy,iz,x,y,z,zmmin,dx,dy,dz,nx,ny,nz,zmesh = self.getmesh(mglevel)
    try:
      self.isinside[0,0,0]
    except:
      self.isinside = zeros((1+nx,1+ny,1+nz),'d')
    ix1 = min(ix)
    ix2 = max(ix)
    iy1 = min(iy)
    iy2 = max(iy)
    tt2[0] = tt2[0] + wtime() - tt1
    if len(x) == 0: return
    for zz in zmesh:
      tt1 = wtime()
      z[:] = zz  #####
      iz[:] = nint((zz - zmmin)/dz)  #####
      tt2[1] = tt2[1] + wtime() - tt1
      tt1 = wtime()
      d = a.isinside(x,y,z)  #####
      tt2[2] = tt2[2] + wtime() - tt1
      tt1 = wtime()
      dd = d.isinside
      dd.shape = (ix2-ix1+1,iy2-iy1+1)
      self.isinside[ix1:ix2+1,iy1:iy2+1,iz[0]] = dd
      tt2[3] = tt2[3] + wtime() - tt1
    endtime = wtime()
    self.generatetime = endtime - starttime
    #print tt2

##############################################################################
##############################################################################
##############################################################################

#============================================================================
class Box(Assembly):
  """
Box class
  - xsize,ysize,zsize: box size
  - voltage=0: box voltage
  - xcent=0.,ycent=0.,zcent=0.: center of box
  - condid=1: conductor id of box, must be integer
  """
  def __init__(self,xsize,ysize,zsize,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist=['xsize','ysize','zsize']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           f3d.boxconductorf,f3d.boxconductord)
    self.xsize = xsize
    self.ysize = ysize
    self.zsize = zsize

#============================================================================
class Cylinder(Assembly):
  """
Cylinder class
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           f3d.cylinderconductorf,f3d.cylinderconductord)
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi

#============================================================================
class Cylinders(Assembly):
  """
Cylinders class for a list of cylinders
  - radius,length: cylinder size
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,theta=0.,phi=0.,
                    voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['ncylinders','radius','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           f3d.cylindersconductorf,f3d.cylindersconductord)
    self.ncylinders = len(radius)
    self.radius = radius
    self.length = length
    self.theta  = theta
    self.phi    = phi

#============================================================================
class ZCylinder(Assembly):
  """
Cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                           f3d.zcylinderconductorf,f3d.zcylinderconductord)
    self.radius = radius
    self.length = length

#============================================================================
class ZRoundedCylinder(Assembly):
  """
Cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,radius2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zroundedcylinderconductorf,
                      f3d.zroundedcylinderconductord)
    self.radius = radius
    self.length = length
    self.radius2 = radius2


#============================================================================
class ZCylinderOut(Assembly):
  """
Outside of a cylinder aligned with z-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zcylinderoutconductorf,f3d.zcylinderoutconductord)
    self.radius = radius
    self.length = length

#============================================================================
class ZRoundedCylinderOut(Assembly):
  """
Outside of a cylinder with rounded corners aligned with z-axis
  - radius,length: cylinder size
  - radius2: radius of rounded corners
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,radius2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length','radius2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zroundedcylinderoutconductorf,
                      f3d.zroundedcylinderoutconductord)
    self.radius = radius
    self.length = length
    self.radius2 = radius2

#============================================================================
class YCylinder(Assembly):
  """
Cylinder aligned with y-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.ycylinderconductorf,f3d.ycylinderconductord)
    self.radius = radius
    self.length = length

#============================================================================
class XCylinder(Assembly):
  """
Cylinder aligned with x-axis
  - radius,length: cylinder size
  - voltage=0: cylinder voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cylinder
  - condid=1: conductor id of cylinder, must be integer
  """
  def __init__(self,radius,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.xcylinderconductorf,f3d.xcylinderconductord)
    self.radius = radius
    self.length = length

#============================================================================
class Sphere(Assembly):
  """
Sphere
  - radius: radius
  - voltage=0: sphere voltage
  - xcent=0.,ycent=0.,zcent=0.: center of sphere
  - condid=1: conductor id of sphere, must be integer
  """
  def __init__(self,radius,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['radius']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.sphereconductorf,f3d.sphereconductord)
    self.radius = radius

#============================================================================
class Cone(Assembly):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.coneconductorf,f3d.coneconductord)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.theta = theta
    self.phi = phi
    self.length = length

#============================================================================
class ConeSlope(Assembly):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.coneconductorf,f3d.coneconductord)
    self.slope = slope
    self.intercept = intercept
    self.theta = theta
    self.phi = phi
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class Cones(Assembly):
  """
Cones
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - theta=0,phi=0: angle of cylinder relative to z-axis
    theta is angle in z-x plane
    phi is angle in z-y plane
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,theta,phi,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['ncones','r_zmin','r_zmax','length','theta','phi']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.conesconductorf,f3d.conesconductord)
    self.ncones = len(r_zmin)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.theta = theta
    self.phi = phi
    self.length = length

#============================================================================
class ZCone(Assembly):
  """
Cone
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zconeconductorf,f3d.zconeconductord)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

#============================================================================
class ZConeSlope(Assembly):
  """
Cone
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,intercept,length,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zconeconductorf,f3d.zconeconductord)
    self.slope = slope
    self.intercept = intercept
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class ZConeOut(Assembly):
  """
Cone outside
  - r_zmin: radius at z min
  - r_zmax: radius at z max
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r_zmin,r_zmax,length,voltage=0.,xcent=0.,ycent=0.,zcent=0.,
                    condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zconeoutconductorf,f3d.zconeoutconductord)
    self.r_zmin = r_zmin
    self.r_zmax = r_zmax
    self.length = length

#============================================================================
class ZConeOutSlope(Assembly):
  """
Cone outside
  - slope: ratio of radius at zmax minus radius at zmin over length
  - intercept: location where line defining cone crosses the axis, relative
               to zcent
  - length: length
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,slope,intercept,length,voltage=0.,
                    xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r_zmin','r_zmax','length']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.zconeoutconductorf,f3d.zconeoutconductord)
    self.slope = slope
    self.intercept = intercept
    self.length = length
  def getkwlist(self):
    self.r_zmin = self.slope*(-self.length/2. - self.intercept)
    self.r_zmax = self.slope*(+self.length/2. - self.intercept)
    return Assembly.getkwlist(self)

#============================================================================
class ZTorus(Assembly):
  """
Torus
  - r1: toroidal radius
  - r2: poloidal radius
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,r1,r2,voltage=0.,xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['r1','r2']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.ztorusconductorf,f3d.ztorusconductord)
    self.r1 = r1
    self.r2 = r2

#============================================================================
class Beamletplate(Assembly):
  """
Plate from beamlet pre-accelerator
  - za: location of spherical center in the x-plane
  - zb: location of spherical center in the y-plane
  - z0: location of the center of the plate on the z-axis
  - thickness: thickness of the plate
  - voltage=0: cone voltage
  - xcent=0.,ycent=0.,zcent=0.: center of cone
  - condid=1: conductor id of cone, must be integer
  """
  def __init__(self,za,zb,z0,thickness,voltage=0.,
               xcent=0.,ycent=0.,zcent=0.,condid=1):
    kwlist = ['za','zb','z0','thickness']
    Assembly.__init__(self,voltage,xcent,ycent,zcent,condid,kwlist,
                      f3d.beamletplateconductorf,f3d.beamletplateconductord)
    self.za = za
    self.zb = zb
    self.z0 = z0
    self.thickness = thickness

  def visualize(self,xmin=None,xmax=None,ymin=None,ymax=None):
    xmin = where(xmin==None,w3d.xmmin,xmin)[0]
    xmax = where(xmax==None,w3d.xmmax,xmax)[0]
    ymin = where(ymin==None,w3d.ymmin,ymin)[0]
    ymax = where(ymax==None,w3d.ymmax,ymax)[0]
    nx = min(w3d.nx,20)
    ny = min(w3d.ny,20)
    nz = w3d.nz
    xmmin = w3d.xmmin
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
    xmmax = w3d.xmmax
    ymmax = w3d.ymmax
    zmmax = w3d.zmmax

    # --- Calculate dx, dy, and dz in case this is called before
    # --- the generate.
    dx = (xmmax - xmmin)/nx
    dy = (ymmax - ymmin)/ny
    dz = (zmmax - zmmin)/nz
    if(w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom):
      dy = dx

    xmesh = xmmin + dx*arange(nx+1)
    ymesh = ymmin + dy*arange(ny+1)
    xmesh = compress(logical_and(xmin-dx <= xmesh,xmesh <= xmax+dx),xmesh)
    ymesh = compress(logical_and(ymin-dy <= ymesh,ymesh <= ymax+dy),ymesh)
    x = ravel(xmesh[:,NewAxis]*ones(len(ymesh)))
    y = ravel(ymesh*ones(len(xmesh))[:,NewAxis])
    ix = nint((x - xmmin)/dx)
    iy = nint((y - ymmin)/dy)
    if len(x) == 0: return

    xx = x
    yy = y
    xx.shape = (len(xmesh),len(ymesh))
    yy.shape = (len(xmesh),len(ymesh))

    # --- Outer face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') - self.thickness
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zl = z[0] + d.dels[5,:]
    zl.shape = (len(xmesh),len(ymesh))
    ml = VPythonobjects.VisualMesh(xx,yy,zl,twoSided=true)

    # --- Inner face
    z = self.z0*ones(len(xmesh)*len(ymesh),'d') + 0.5*self.za
    iz = nint((z - zmmin)/dz)
    d = self.griddistance(ix,iy,iz,x,y,z)
    zr = z[0] - d.dels[4,:]
    zr.shape = (len(xmesh),len(ymesh))
    mr = VPythonobjects.VisualMesh(xx,yy,zr,twoSided=true)

    # --- Four sides between faces
    xside = xx[:,0]*ones(2)[:,NewAxis]
    yside = yy[:,0]*ones(2)[:,NewAxis]
    zside = array([zl[:,0],zr[:,0]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[:,-1]*ones(2)[:,NewAxis]
    yside = yy[:,-1]*ones(2)[:,NewAxis]
    zside = array([zl[:,-1],zr[:,-1]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[0,:]*ones(2)[:,NewAxis]
    yside = yy[0,:]*ones(2)[:,NewAxis]
    zside = array([zl[0,:],zr[0,:]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

    xside = xx[-1,:]*ones(2)[:,NewAxis]
    yside = yy[-1,:]*ones(2)[:,NewAxis]
    zside = array([zl[-1,:],zr[-1,:]])
    ms1 = VPythonobjects.VisualMesh(xside,yside,zside,twoSided=true)

#============================================================================
#############################################################################

