"""Implements adaptive mesh refinement in 3d
"""
from warp import *
from multigrid import MultiGrid
from fieldsolver import SubcycledPoissonSolver
from find_mgparam import find_mgparam
from pyOpenDX import Visualizable,DXCollection,viewboundingbox
import MA
import __main__
import time
#import threading
try:
  import psyco
except ImportError:
  pass

#########################################################################
# Note that MRBlock is psyco.bind at the end of the file
class MRBlock(MultiGrid,Visualizable):
  """
Implements adaptive mesh refinement in 3d
 - parent:
 - refinement=None: amount of refinement along each axis
 - lower,upper: extent of domain in relative to parent, in its own grid
                cell size, and not including any guard cells
 - dims: dimensions of the grid, only used for root block, the one with
         no parents
 - mins,maxs: locations of the grid lower and upper bounds in the beam frame
 - root: coarsest level grid
 - children: list of tuples, each containing three elements,
             (lower,upper,refinement). Children can also be added later
             using addchild.
 - lreducedpickle=true: when true, a small pickle is made by removing all of
                        the big arrays. The information can be regenerated
                        upon restart.
  """
  def __init__(self,parent=None,refinement=None,
                    lower=None,upper=None,
                    fulllower=None,fullupper=None,
                    dims=None,mins=None,maxs=None,
                    nguard=[1,1,1],
                    children=None,lreducedpickle=1,**kw):

    # --- Pass the input value of lreducedpickle into the MultiGrid class
    kw['lreducedpickle'] = lreducedpickle

    if parent is None:

      # --- No parents, so just create empty lists
      self.parents = []
      self.root = self

      # --- It is assumed that the root block will be the first one created.
      # --- So clear out the global block list and count.
      self.totalnumberofblocks = 0
      self.listofblocks = []
      self.finalized = 0
      self.my_index = me

      # --- For the root, the dimensions and extent of the grid should
      # --- be specified. If not, they will be taken from w3d.
      self.dims = dims
      self.mins = mins
      self.maxs = maxs
      self.totalrefinement = ones(3)
      self.forcesymmetries = 1
      self.nguard = nguard
      self.deltas = None

      self.refinement = None

      # --- top.nparpgrp is now set in the init of the base class.
      
    else:

      # --- Save the parent and the index number. These are saved in lists
      # --- since a block can have multiple parents.
      self.parents = [parent.blocknumber]
      self.root = parent.root
      self.my_index = 0

      # --- Make sure that refinement is an array of length three. If a scalar
      # --- is input, it is broadcast to all three axis.
      if type(refinement) in [IntType,FloatType]:
        refinement = 3*[refinement]
      self.refinement = array(refinement)
      if type(nguard) in [IntType,FloatType]:
        nguard = 3*[nguard]
      self.nguard = array(nguard)

      self.totalrefinement = parent.totalrefinement*self.refinement
      self.deltas = parent.deltas/self.refinement
      self.rootdims = self.root.dims*self.totalrefinement
      self.forcesymmetries = 0

      if lower is None and upper is None:
        # --- The grid mins and maxs are input.
        # --- The lower and upper are calculated to be an integer number of
        # --- parent grid cells times the refinement factor. The lower is
        # --- rounded down and upper rounded up to ensure that the patch
        # --- includes the entire extent specified by mins and maxs.
        self.mins = array(mins)
        self.maxs = array(maxs)
        self.mins = maximum(self.mins,self.root.mins)
        self.maxs = minimum(self.maxs,self.root.maxs)
        self.lower = (nint(floor((self.mins - self.root.minsglobal)/parent.deltas))*
                     self.refinement)
        self.upper = (nint(ceil((self.maxs - self.root.minsglobal)/parent.deltas))*
                     self.refinement)
        self.lower = maximum(self.root.lower*self.totalrefinement,self.lower)
        self.upper = minimum(self.root.upper*self.totalrefinement,self.upper)

      else:
        # --- The grid lower and upper bounds are input. The bounds are
        # --- relative to the root grid, but scaled by the total refinement.
        self.lower = array(lower)
        self.upper = array(upper)

      # --- In parallel, an extra grid cell in z can be added since the
      # --- information in the z guard planes of phi is correct.
      self.extradimslower = zeros(3)
      self.extradimsupper = zeros(3)
      if me > 0:      self.extradimslower[-1] = 1
      if me < npes-1: self.extradimsupper[-1] = 1

      rootfulllower = self.root.fulllower*self.totalrefinement - self.extradimslower*self.totalrefinement
      rootfullupper = self.root.fullupper*self.totalrefinement + self.extradimsupper*self.totalrefinement

      # --- Now, extend the domain by the given number of guard cells. Checks
      # --- are made so that the domain doesn't extend beyond the root grid.
      if fulllower is None:
        self.fulllower = maximum(rootfulllower,self.lower - self.nguard*self.refinement)
      else:
        self.fulllower = array(fulllower)
      if fullupper is None:
        self.fullupper = minimum(rootfullupper,self.upper + self.nguard*self.refinement)
      else:
        self.fullupper = array(fullupper)

      # --- Get the number of grid points along each dimension
      self.dims = self.fullupper - self.fulllower

      # --- Make sure that the number of grid points is even.
      # --- If it is odd, then enough cells are added to extend to the next
      # --- grid cell of the parent. It is then cutoff at the root grid.
      self.fulllower = where(self.dims%2==1,self.fulllower-self.refinement,
                                            self.fulllower)
      self.fulllower = maximum(rootfulllower,self.fulllower)
      self.dims = self.fullupper - self.fulllower

      # --- If it is still odd (which means that the cells added above
      # --- where cutoff at zero) then add some at the top.
      self.fullupper = where(self.dims%2==1,self.fullupper+self.refinement,
                                            self.fullupper)
      self.fullupper = minimum(rootfullupper,self.fullupper)
      self.dims = self.fullupper - self.fulllower

      # --- If it is still odd, then there is some serious problem. The number
      # --- in the base grid may be odd.
      assert alltrue(self.dims%2 == 0),\
             """The number of grid cells in one of the dimensions is odd - they
             all must be even. Check that the number of cells in the base grid
             is even."""

      # --- Now calculate the extent of the grid
      self.mins = self.root.minsglobal + self.fulllower*self.deltas
      self.maxs = self.root.minsglobal + self.fullupper*self.deltas

      # --- Make sure there is some overlap of the child with the parent
      mins = maximum(self.mins,parent.mins)
      maxs = minimum(self.maxs,parent.maxs)
      assert alltrue(maxs >= mins),\
             "The child is not within the extent of the parent"

      # --- Make sure that the patch has a finite extent in all dimensions
      assert alltrue(self.upper>self.lower),\
             "The child must have a finite extent in all dimensions"

      # --- First, just use same boundary conditions as root.
      self.bounds = self.root.bounds.copy()
      self.pbounds = self.root.pbounds.copy()

      # --- Check if the mesh doesn't reach the edge of the root grid.
      # --- If not, switch to Dirichlet boundary.
      self.bounds[::2] = where(self.fulllower > 0,0,self.bounds[::2])
      self.bounds[1::2] = where(self.fullupper < self.root.dimsglobal*self.totalrefinement,
                                0,self.bounds[1::2])
      self.l2symtry = self.root.l2symtry
      self.l4symtry = self.root.l4symtry
      self.pbounds[::2] = where(self.fulllower > 0,0,self.pbounds[::2])
      self.pbounds[1::2] = where(self.fullupper < self.root.dimsglobal*self.totalrefinement,
                                0,self.pbounds[1::2])

      # --- Create some temporaries for optimization
      self.fullloweroverrefinement = self.fulllower/self.refinement
      self.fullupperoverrefinement = self.fullupper/self.refinement

    # --- Set individual quantities based on the values in the arrays,
    # --- if they have been set.
    if self.deltas is not None:
      self.dx = self.deltas[0]
      self.dy = self.deltas[1]
      self.dz = self.deltas[2]
    if self.dims is not None:
      self.nx = self.dims[0]
      self.ny = self.dims[1]
      self.nz = self.dims[2]
    if self.mins is not None:
      self.xmmin = self.mins[0]
      self.ymmin = self.mins[1]
      self.zmmin = self.mins[2]
    if self.maxs is not None:
      self.xmmax = self.maxs[0]
      self.ymmax = self.maxs[1]
      self.zmmax = self.maxs[2]

    # --- Do some further initialization.
    MultiGrid.__init__(self,**kw)

    if parent is None:
      # --- This is only needed by the root grid in cases when the grid
      # --- parameters are obtained from w3d instead of the argument list.
      self.dims = array([self.nx,self.ny,self.nz])
      self.dimsglobal = array([self.nx,self.ny,self.nzfull])
      self.deltas = array([self.dx,self.dy,self.dz])
      self.mins = array([self.xmmin,self.ymmin,self.zmmin])
      self.maxs = array([self.xmmax,self.ymmax,self.zmmax])
      self.minsglobal = array([self.xmmin,self.ymmin,self.zmminglobal])
      self.maxsglobal = array([self.xmmax,self.ymmax,self.zmmaxglobal])
      self.lower = nint(((self.mins - self.minsglobal)/self.deltas))
      self.upper = nint(((self.maxs - self.minsglobal)/self.deltas))
      self.fulllower = self.lower.copy()
      self.fullupper = self.upper.copy()
      self.rootdims = self.dims

    # --- childdomains is the node centered grid which keeps track of which
    # --- cells are owned by which children. If there are no children,
    # --- then it is not needed.
    self.childdomains = None

    # --- Set up variable time steps
    if top.chdtspid > 0:
      top.dxpid = nextpid()
      top.dypid = nextpid()
      top.dzpid = nextpid()

    # --- Get the current global block number and increment the counter.
    # --- Also, add self to the global list of blocks.
    self.blocknumber = self.root.totalnumberofblocks
    self.root.totalnumberofblocks += 1
    self.root.listofblocks.append(self)

    # --- Note that a dictionary is used for the overlaps so that lookups
    # --- are faster, also, the values in the dictionary are list containing
    # --- the domain of the overlap. The blocks with lower and higher block
    # --- number are treated differently, so create separate lists for each.
    # --- Two separate lists is likely only a small optimization.
    self.overlapslower = {}
    self.overlapshigher = {}
    self.overlapsparallelleft = {}
    self.overlapsparallelright = {}

    # --- Now add any specified children
    self.children = []
    if children is not None:
      for l,u,r in children:
        self.addchild(l,u,refinement=r)

  def __getstate__(self):
    """
Check whether this instance is the registered solver so that upon unpickling
it knows whether to re-register itself.
    """
    # --- Make sure the the lreducedpickle option gets propagated to all
    # --- of the blocks.
    for child in self.children:
      child.lreducedpickle = self.lreducedpickle
    dict = MultiGrid.__getstate__(self)
    if self.lreducedpickle:
      # --- Remove the big objects from the dictionary. This can be
      # --- regenerated upon the restore.
      dict['childdomains'] = None
    return dict

  def __setstate__(self,dict):
    MultiGrid.__setstate__(self,dict)
   #self.makefortranordered('phi')
   #self.makefortranordered('rho')
   #self.makefortranordered('selfe')
    if (self == self.root and self.lreducedpickle and
        not self.lnorestoreonpickle):
      # --- It is assumed that at this point, all of the children have been
      # --- restored.
      # --- Regenerate childdomains
      self.initializechilddomains()
      # --- If rho and phi weren't saved, make sure that they are setup.
      # --- Though, this may not always be the right thing to do.
      # --- These can only be done at the end of the restart since only then
      # --- is it gauranteed that the particles are read in.
      installafterrestart(self.loadrho)
      installafterrestart(self.solve)

  def makefortranordered(self,vname):
    a = getattr(self,vname)
    if type(a) is ArrayType:
      setattr(self,vname,fzeros(shape(a),a.typecode()))
      getattr(self,vname)[...] = a

  def addchild(self,lower=None,upper=None,fulllower=None,fullupper=None,
                    mins=None,maxs=None,
                    refinement=[2,2,2],nguard=None,nslaves=1):
    """
Add a mesh refined block to this block.
  -lower,upper,mins,maxs,refinement: All have same meanings as for the
                                     constructor.
  -nslaves=1: defaults to one so it is not parallelized
    """
    try:
      if nguard is None: nguard = self.nguard
      child = MRBlock(parent=self,lower=lower,upper=upper,
                      fulllower=fulllower,fullupper=fullupper,
                      mins=mins,maxs=maxs,
                      refinement=refinement,nguard=nguard,
                      nslaves=nslaves)
      self.children.append(child)
      return child
    except:
      if self.root.nslaves <= 1: raise

  def resetroot(self):
    # --- No parents, so just create empty lists
    self.parents = []
    self.root = self
    self.totalnumberofblocks = 1
    self.listofblocks = [self]
    self.childdomains = None
    self.children = []

  #--------------------------------------------------------------------------
  # --- The next several methods handle conductors
  #--------------------------------------------------------------------------

  def installconductor(self,conductor,dfill=top.largepos):
    if not self.isfirstcall(): return
    MultiGrid.installconductor(self,conductor,dfill=dfill)
    for child in self.children:
      child.installconductor(conductor,dfill=dfill)

  def clearconductors(self):
    if not self.isfirstcall(): return
    MultiGrid.clearconductors(self)
    for child in self.children:
      child.clearconductors()

  def hasconductors(self):
    return (self.conductors.interior.n > 0 or
            self.conductors.evensubgrid.n > 0 or
            self.conductors.oddsubgrid.n > 0)

  def getconductors(self,alllevels=1,result=None):
    if result is None: result = []
    result.append(self.conductors)
    if alllevels:
      for child in self.children:
        child.getconductors(alllevels,result)
    return result

  def setconductorvoltage(self,voltage,condid=0,discrete=false,
                          setvinject=false):
    'Recursively calls setconductorvoltage for base and all children'
    if not self.isfirstcall(): return
    setconductorvoltage(voltage,condid,discrete,setvinject,
                        conductors=self.conductors)
    for child in self.children:
      child.setconductorvoltage(voltage,condid,discrete)

  #--------------------------------------------------------------------------
  # --- The next several methods handle initialization that is done after
  # --- all blocks have been added.
  #--------------------------------------------------------------------------

  def finalize(self):
    # --- This should only be called at the top level.
    if self != self.root or self.finalized: return
    blocklists = self.generateblocklevellists()
    self.blocklists = blocklists
    blocklistsleft,blocklistsright = self.swapblocklistswithprocessneighbors(blocklists)
    self.clearparentsandchildren()
    self.findallchildren(blocklists)
    self.initializechilddomains()
    self.findoverlappingsiblings(blocklists[1:],
                                 blocklistsleft[1:],blocklistsright[1:])
    self.finalized = 1

  def generateblocklevellists(self,blocklists=None):
    if blocklists is None:
      # --- This will only happen at the top level.
      # --- Create a list of empty lists. Each empty list will get the blocks
      # --- at the appropriate levels appended to it. Note that 100 is
      # --- assumed to be a large enough number - there almost certainly
      # --- will never be 100 levels of refinement.
      blocklists = [[] for i in range(100)]
    # --- Add this instance to the top level of the list and pass the rest
    # --- of it to the children
    if self not in blocklists[0]:
      blocklists[0].append(self)
      for child in self.children:
        b = child.generateblocklevellists(blocklists[1:])
    return blocklists

  def swapblocklistswithprocessneighbors(self,blocklists):
    if not lparallel:
      blocklistsleft = [[] for i in range(100)]
      blocklistsright = [[] for i in range(100)]
    else:
      # --- First, set so amount data sent to neighbors will be small
      lreducedpicklesave = self.lreducedpickle
      for block in self.listofblocks:
        self.lreducedpickle = 1
        self.lnorestoreonpickle = 1
      if (me > 0     ): mpi.send(blocklists,me-1)
      if (me < npes-1): blocklistsright = mpi.recv(me+1)[0]
      else:             blocklistsright = [[] for i in range(100)]
      if (me < npes-1): mpi.send(blocklists,me+1)
      if (me > 0     ): blocklistsleft = mpi.recv(me-1)[0]
      else:             blocklistsleft = [[] for i in range(100)]
      # --- Restore the flags
      for block in self.listofblocks:
        self.lreducedpickle = lreducedpicklesave
        self.lnorestoreonpickle = 0
    return blocklistsleft,blocklistsright
      
  def clearparentsandchildren(self):
    self.parents = []
    for child in self.children:
      child.clearparentsandchildren()
    self.children = []

  def findallchildren(self,blocklists):
    for block in blocklists[1]:
      # --- Get extent of possible overlapping domain
      l = maximum(block.fullloweroverrefinement,self.fulllower)
      u = minimum(block.fullupperoverrefinement,self.fullupper)
      #if alltrue(u >= l):
      if (u[0] >= l[0] and u[1] >= l[1] and u[2] >= l[2]):
        self.children.append(block)
        block.parents.append(self.blocknumber)

    # --- Only the first block in the list makes the call for the next level.
    # --- This guarantees that this method is called only once for each block.
    if blocklists[0][0] == self:
      for block in blocklists[1]:
        block.findallchildren(blocklists[1:])

  def initializechilddomains(self):
    """
Sets the regions that are covered by the children.
    """
    if not self.isfirstcall(): return
    # --- Loop over the children, first calling each, then setting
    # --- childdomain appropriately.
    for child in self.children:
      child.initializechilddomains()

      # --- Set full domain to negative of child number first.
      l = maximum(self.fulllower,child.fullloweroverrefinement)
      u = child.fullupperoverrefinement
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims all unclaimed areas.
      ii = self.getchilddomains(l,u)
      ii[...] = where(ii==self.blocknumber,-child.blocknumber,ii)

      # --- Set interior to positive child number.
      l = maximum(self.fulllower,child.lower/child.refinement)
      u = child.upper/child.refinement
      # --- Check against the case where only the guard cells of the child
      # --- overlap the parent.
      if (u[0] < l[0] or u[1] < l[1] or u[2] < l[2]): continue
      # --- If the child extends to the edge of the parent mesh, it claims the
      # --- grid points on the upper edges.
      u = u + where(u == self.fullupper,1,0)
      # --- The child claims its full interior area
      ii = self.getchilddomains(l,u)
      ii[...] = +child.blocknumber

  def findoverlappingsiblings(self,blocklists,blocklistsleft,blocklistsright):
    """
Recursive routine to find, at each level of refinement, all overlapping
siblings.
This is faster than the original routine since each pair of blocks is checked
only once, rather than twice for each parent as in the original.
    """
    # --- When the list is empty, there are no blocks, so just return.
    if (len(blocklists[0]) == 0 and
        len(blocklistsleft[0]) == 0 and
        len(blocklistsright[0]) == 0): return
    # --- Make the call for the next level.
    self.findoverlappingsiblings(blocklists[1:],
                                 blocklistsleft[1:],blocklistsright[1:])

    # --- Get a copy of the list (which will be mangled below).
    blocklistscopy = copy.copy(blocklists[0])

    # --- Loop over all blocks.
    for block in blocklists[0]:

      # --- Loop only over blocks that havn't been checked yet.
      del blocklistscopy[0]
      for sibling in blocklistscopy:

        # --- Get the area common to the block and its sibling.
        # --- Don't do anything if there is no overlap.
        sl = maximum(sibling.fulllower,block.fulllower)
        su = minimum(sibling.fullupper,block.fullupper)
        if sl[0] > su[0] or sl[1] > su[1] or sl[2] > su[2]: continue

        # --- The ordering ofthe blocks in the list matter since lower blocks
        # --- get precedence, and the rho is accumualated there.
        if block.blocknumber < sibling.blocknumber:
          block.overlapshigher[sibling.blocknumber] = [sl,su]
          sibling.overlapslower[block.blocknumber] = [sl,su]
        else:
          block.overlapslower[sibling.blocknumber] = [sl,su]
          sibling.overlapshigher[block.blocknumber] = [sl,su]

    # --- Now find overlapping blocks in neighboring processors. 
    # --- For the precedence, even numbered processors get precdence
    # --- over odd numbered processors. Otherwise, the logic is the same
    # --- as above.
    for block in blocklists[0]:
      for neighborlists,pe in zip([blocklistsleft[0],blocklistsright[0]],[me-1,me+1]):
        if pe < 0 or pe == npes: continue
        for neighborblock in neighborlists:
          sl = maximum(neighborblock.fulllower,block.fulllower)
          su = minimum(neighborblock.fullupper,block.fullupper)
          if sl[0] >= su[0] or sl[1] >= su[1] or sl[2] >= su[2]: continue
          if pe < me:
            block.overlapsparallelleft[neighborblock.blocknumber] = [sl,su,pe]
          else:
            block.overlapsparallelright[neighborblock.blocknumber] = [sl,su,pe]

  def clearinactiveregions(self,nbcells,parent=None,level=1):
    """
For regions which should not be refined but are included because of the
coalescing of blocks, the childdomains is set so that the E-fields will
not be fetched from there (it is set negative).
 - nbcells: array containing the refinement level of the grid cells.
            The shape will be the same as the intersection of self and the
            calling parent.
 - parent: the calling parent
 - level: the total amount refinement
    """
    if parent is None:
      # --- If there is no parent, it is assumed that this is the root block
      l = self.fulllower
      u = self.fullupper
    else:
      # --- Find intersection of parent and self.
      l = maximum(parent.fulllower*self.refinement,self.fulllower)
      u = minimum(parent.fullupper*self.refinement,self.fullupper)

    for child in self.children:
      # --- Find intersection of parent, self, and child
      cl = maximum(child.fullloweroverrefinement,l)
      cu = minimum(child.fullupperoverrefinement,u)
      #if sometrue(cl > cu): continue
      if cl[0] > cu[0] or cl[1] > cu[1] or cl[2] > cu[2]: continue

      # --- Get childdomains in the intersection region, Wherever the
      # --- the refinement level is lower than the childs, force childdomains
      # --- to be negative.
      ii = self.getchilddomains(cl,cu,1)
      nbc = self.getlocalarray(nbcells,cl,cu,fulllower=l)
      # --- Note that if the upper edge of the child does not extend to the
      # --- upper edge of self, then the childdomains at cu will have the
      # --- blocknumber of self, and so should  not be reset. The check
      # --- for (ii == child.blocknumber) prevents this. The check also
      # --- prevents one child from clearing out the domain of another, though
      # --- that wouldn't break anything since that other child would be
      # --- clearing out the same region itself. Also, with the check,
      # --- the second argument of the where can just be -ii instead
      # --- of -abs(ii) since only positive values of ii will have the sign
      # --- changed.
      r = child.refinement
      ii[...] = where((nbc<max(level*r)) & (ii == child.blocknumber),-ii,ii)

      # --- Stretch out the array so it has the refined cell size of the child
      nbcstretched = zeros(1+(cu-cl)*r)
      for k in range(r[2]):
        if k == 0: ksl = slice(None)
        else:      ksl = slice(-1)
        for j in range(r[1]):
          if j == 0: jsl = slice(None)
          else:      jsl = slice(-1)
          for i in range(r[0]):
            if i == 0: isl = slice(None)
            else:      isl = slice(-1)
            nbcstretched[i::r[0],j::r[1],k::r[2]] = nbc[isl,jsl,ksl]

      child.clearinactiveregions(nbcstretched,self,level*r)

  #--------------------------------------------------------------------------
  # --- The next several methods handle the charge density calculation.
  #--------------------------------------------------------------------------

  def setrhopforparticles(self,*args):
    for block in self.listofblocks:
      SubcycledPoissonSolver.setrhopforparticles(block,*args)

  def allocatedataarrays(self):
    # --- If not root, than only allocate the arrays of this patch
    if self != self.root:
      SubcycledPoissonSolver.allocatedataarrays(self)
      return
    # --- Otherwise, do all patches.
    # --- Make sure that the final setup was done. This is put here
    # --- since this routine is called by loadrho and solve, which
    # --- call finalize anyway.
    self.finalize()
    # --- Now loop over blocks, calling allocatedataarrays
    for block in self.listofblocks:
      SubcycledPoissonSolver.allocatedataarrays(block)

  def zerorhop(self):
    for block in self.listofblocks:
      SubcycledPoissonSolver.zerorhop(block)

  def averagerhopwithsubcycling(self):
    for block in self.listofblocks:
      SubcycledPoissonSolver.averagerhopwithsubcycling(block)

  def setrhop(self,x,y,z,uz,q,w,zgrid,lrootonly=0):
    """
Given the list of particles, a charge and a weight, deposits the charge
density of the mesh structure.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
The sort still takes up about 40% of the time. Note that this depends on
having ichilddomains filled with the blocknumber rather than the child number
relative to the parent.
    """
    if len(self.children) > 0 and not lrootonly:

      ichild = zeros(len(x))
      self.getichild(x,y,z,ichild)

      x,y,z,uz,nperchild = self.sortbyichild(ichild,x,y,z,uz)

    else:
      nperchild = [len(x)]

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MultiGrid.setrhop(block,x[i:i+n],y[i:i+n],z[i:i+n],uz[i:i+n],q,w,zgrid)
      i = i + n

  def aftersetrhop(self,lzero):
    # --- distribute charge density among blocks
    if lzero:
      # --- propagate rhop between patches
      self.getrhopfromoverlaps()
      self.gatherrhopfromchildren()
      self.exchangerhopwithneighbors()
      self.restorerhopinoverlaps()

  def exchangerhopwithneighbors(self):
    """
Exchange rhop in blocks overlapping blocks on neighboring processors.
    """
    assert self is self.root,"This should only be called by the root block"
    if npes == 0: return

    # --- Make a list of nonblocking sends to wait for to finish.
    parallelwaitlist = []

    # --- All blocks first send the overlapping rhop to the neighbors.
    # --- Note that the precedences for rhop is handled locally so
    # --- here, all exchanges are made. Note that all data to be
    # --- sent is gathered into a dictionary and sent together.
    # --- Send to the left first
    if me > 0:
      senddict = {}
      for block in self.listofblocks:
        for othernumber,data in block.overlapsparallelleft.items():
          l,u,pe = data
          rhop = block.getrhop(l,u)
          senddict.setdefault(othernumber,[]).append((l,u,rhop))
      parallelwaitlist.append(mpi.isend(senddict,me-1))
    # --- Then to the right
    if me < npes-1:
      senddict = {}
      for block in self.listofblocks:
        for othernumber,data in block.overlapsparallelright.items():
          l,u,pe = data
          rhop = block.getrhop(l,u)
          senddict.setdefault(othernumber,[]).append((l,u,rhop))
      parallelwaitlist.append(mpi.isend(senddict,me+1))

    # --- The parallel receives are done afterward since they are blocking. At
    # --- this point, all of the sends have been made so there should be little
    # --- waiting and no contention.
    # --- First receive the data from the right
    if me < npes-1:
      recvdict = mpi.recv(me+1)[0]
      for blocknumber,data in recvdict.items():
        block = self.getblockfromnumber(blocknumber)
        for l,u,orhop in data:
          srhop = block.getrhop(l,u)
          add(srhop,orhop,srhop)
    # --- The from the left
    if me > 0:
      recvdict = mpi.recv(me-1)[0]
      for blocknumber,data in recvdict.items():
        block = self.getblockfromnumber(blocknumber)
        for l,u,orhop in data:
          srhop = block.getrhop(l,u)
          add(srhop,orhop,srhop)

    # --- This is not really necessary, but helps keeps things clean.
    if len(parallelwaitlist) > 0: mpi.waitall(parallelwaitlist)

  def getrhopfromoverlaps(self):
    """
Add in the rhop from overlaping areas. The rhop is gathered into the block with
the lowerest number. Later on, the rhop will be copied back to the higher
numbered blocks. Note that overlaps from neighboring processors has already
been taken care of. This should only ever be called by the root block.
    """
    assert self is self.root,"This should only be called by the root block"

    # --- This loops over the blocks in ascending order to ensure that in any
    # --- area with overlap, the block with the lowest number is the one that
    # --- gets the rhop. This avoids problems of double counting rhop. This
    # --- could also be done be zeroing out orhop, but that is extra
    # --- (unecessary) computational work, since it already will be done
    # --- at the end of gatherrhopfromchildren.
    for block in self.listofblocks:
      for othernumber,overlapdomain in block.overlapshigher.items():
        other = block.getblockfromnumber(othernumber)
        l,u = overlapdomain
        srhop = block.getrhop(l,u)
        orhop = other.getrhop(l,u)
        add(srhop,orhop,srhop)
        orhop[...] = 0.

  def sortbyichild(self,ichild,x,y,z,uz):
    xout,yout,zout,uzout = zeros((4,len(x)),'d')
    nperchild = zeros(self.root.totalnumberofblocks)
    sortparticlesbyindex(len(x),ichild,x,y,z,uz,self.root.totalnumberofblocks,
                         xout,yout,zout,uzout,nperchild)
    return xout,yout,zout,uzout,nperchild

  def getichild(self,x,y,z,ichild):
    """
Gathers the ichild for the setrhop.
    """
    # --- This must wait until all of the parents have have set ichild
    # --- so that the value in the children takes precedence.
    if not self.islastcall(): return
    if len(x) == 0: return
    if len(self.children) > 0:
      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      getichild(self.blocknumber,len(x),x,y,z,ichild,
                self.nx,self.ny,self.nz,self.childdomains,
                self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                self.zmmin,self.zmmax,top.zgrid,
                self.l2symtry,self.l4symtry)
      for child in self.children:
        child.getichild(x,y,z,ichild)

  def zerorhopinoverlap(self):
    """
This zeros out rhop in overlapping regions for higher numbered blocks.  When
rhop is passed from child to parent, in any overlapping regions only one child
needs to pass the data to the parent.  The choice as to which does the
passing is determined by the blocknumber - the lower gets to do the passing.
For the others, the rhop in the overlapping region is cleared out. That rhop
will be restored later by a call to restorerhopinoverlaps.
Note that this is not recursive, since it is called separately by each block
from gatherrhopfromchildren.
    """
    for othernumber,overlapdomain in self.overlapslower.items():
      l,u = overlapdomain
      srhop = self.getrhop(l,u)
      srhop[...] = 0.

  def gatherrhopfromchildren(self):
    """
Fortran version
    """
    # --- Do this only the first time this is called. This should only be
    # --- done once and since each parent requires that this be done
    # --- before it can get its rhop from here, it must be done on the
    # --- first call.
    if not self.isfirstcall(): return

    # --- Loop over the children
    for child in self.children:

      # --- Make sure that the child has gathered rhop from its children.
      child.gatherrhopfromchildren()

      # --- Get coordinates of child relative to this domain
      l = maximum(child.fullloweroverrefinement,self.fulllower)
      u = minimum(child.fullupperoverrefinement,self.fullupper)

      # --- Check for any Nuemann boundaries
     #dopbounds = (sometrue(child.pbounds == 1) and
     #            (sometrue(l == 0) or
     #             sometrue(u == self.rootdims)))
      dopbounds = ((child.pbounds[0] == 1 or
                    child.pbounds[1] == 1 or
                    child.pbounds[2] == 1) and
                  (l[0] == 0 or l[1] == 0 or l[2] == 0 or
                   u[0] == self.rootdims[0] or
                   u[1] == self.rootdims[1] or
                   u[2] == self.rootdims[2]))

      w = self.getwarrayforrhop(child.refinement)
      gatherrhofromchild(self.rhop,self.dims,child.rhop,child.dims,
                         l,u,self.fulllower,child.fulllower,child.fullupper,
                         child.refinement,w,
                         dopbounds,child.pbounds,self.rootdims)

    # --- zerorhopinoverlap is call here so that any contribution from
    # --- the children in the overlap regions will get zeroed as necessary.
    self.zerorhopinoverlap()

  def restorerhopinoverlaps(self):
    """
Restore rhop in overlapping areas for blocks which had the rhop zeroed out, the
higher numbered blocks. This should only ever be called by the root block.
    """
    assert self is self.root,"This should only be called by the root block"
    # --- The loop does not need to be in ascending order, but this just
    # --- matches the getrhopfromoverlaps routine.
    for block in self.listofblocks:
      for othernumber,overlapdomain in block.overlapslower.items():
        other = block.getblockfromnumber(othernumber)
        l,u = overlapdomain
        srhop = block.getrhop(l,u)
        orhop = other.getrhop(l,u)
        srhop[...] = orhop

  #--------------------------------------------------------------------------
  # --- Methods to carry out the field solve
  #--------------------------------------------------------------------------

  def dosolve(self,iwhich=0,*args):

    # --- Make sure that the final setup was done.
    self.finalize()

    # --- Wait until all of the parents have called here until actually
    # --- doing the solve. This ensures that the phi in all of the parents
    # --- which is needed on the boundaries will be up to date.
    if not self.islastcall(): return

    # --- solve on phi, first getting phi from the parents - both the
    # --- boundary conditions and the interior values as the initial
    # --- value.
    self.setphifromparents()
    MultiGrid.dosolve(self,iwhich)

    # --- solve for children, using the routine which does the correct
    # --- referencing for subcycling and self-B correction
    for child in self.children:
      child.dosolveonphi(iwhich,*args)
    """
    if self == self.root:
      t = threading.Thread(target=MultiGrid.dosolve,args=tuple([self,iwhich]))
      t.start()
      t.join()
      for blocklist in self.blocklists[1:]:
        if len(blocklist) == 0: break
        tlist = []
        for block in blocklist:
          t = threading.Thread(target=block.dosolveonphi,args=tuple([iwhich]+list(args)))
          t.start()
          tlist.append(t)
          print self.blocknumber,len(tlist)
        for t in tlist:
          t.join()
    else:
      self.setphifromparents()
      MultiGrid.dosolve(self,iwhich)
    """


  def setphifromparents(self):
    """
Sets phi, using the values from the parent grid. Setting the full phi array
gives a better initial guess for the field solver.
    """
    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      # --- Coordinates of mesh relative to parent's mesh location
      # --- and refinement. The minimum and maximum are needed in case
      # --- this mesh extends beyond the parent's.
      plower = parent.fulllower*self.refinement - self.extradimslower*self.refinement
      pupper = parent.fullupper*self.refinement + self.extradimsupper*self.refinement
      l = maximum(plower,self.fulllower)
      u = minimum(pupper,self.fullupper)
      # --- The full phi arrays are passed in to avoid copying the subsets
      # --- since the fortran needs contiguous arrays.
      gatherphifromparents(self.phi,self.dims,l,u,self.fulllower,
                           parent.phi,parent.dims,parent.fulllower,
                           self.refinement)

  def optimizeconvergence(self,resetpasses=1):
    if not self.isfirstcall(): return
    MultiGrid.optimizeconvergence(self,resetpasses=resetpasses)
    for child in self.children:
      child.optimizeconvergence(resetpasses=resetpasses)
    
  #--------------------------------------------------------------------------
  # --- Methods to fetch E-fields and potential
  #--------------------------------------------------------------------------

  def fetchefrompositions(self,x,y,z,ex,ey,ez,pgroup=None):
    if pgroup is None:
      self.fetchefrompositionswithoutpgroup(x,y,z,ex,ey,ez)
    else:
      self.fetchefrompositionswithpgroup(x,y,z,ex,ey,ez,pgroup)

  def fetchefrompositionswithpgroup(self,x,y,z,ex,ey,ez,pgroup=None):
    """
Given the list of particles, fetch the E fields.
This first gets the blocknumber of the block where each of the particles are
to be deposited. This is then sorted once. The loop is then over the list
of blocks, rather than walking through the tree structure.
The sort takes up about 40% of the time. It is significantly faster
using the fortran sort.
Note that this depends on having ichilddomains filled with the
blocknumber rather than the child number relative to the parent.
Also, this ends up with the input data remaining sorted.
    """
    if len(self.children) > 0:

      ichild = zeros(len(x))
      # --- This assumes that the root block has blocknumber zero.
      self.getichild_positiveonly(x,y,z,ichild)

      # --- This sorts the particle data in place, including
      # --- the velocities, gaminv, and pid.
      nn = self.root.totalnumberofblocks
      nperchild = zeros(nn)
      particlesortbyindex(pgroup,ichild,0,w3d.ipminfsapi,w3d.npfsapi,
                          nn,nperchild)

    else:
      nperchild = [len(x)]

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MultiGrid.fetchefrompositions(block,x[i:i+n],y[i:i+n],z[i:i+n],
                                          ex[i:i+n],ey[i:i+n],ez[i:i+n])
      if top.chdtspid > 0:
        pgroup.pid[i:i+n,top.dxpid-1] = block.dx
        pgroup.pid[i:i+n,top.dypid-1] = block.dy
        pgroup.pid[i:i+n,top.dzpid-1] = block.dz
      i = i + n

  def fetchefrompositionswithoutpgroup(self,x,y,z,ex,ey,ez):
    """
This is the old version of fetchefrompositions that doesn't rely on having
access to the particle group and does not sort the input data.
    """
    if len(self.children) > 0:

      ichild = zeros(len(x))
      # --- This assumes that the root block has blocknumber zero.
      self.getichild_positiveonly(x,y,z,ichild)

      x,y,z,isort,nperchild = self.sortbyichildgetisort(ichild,x,y,z)

      # --- Create temporary arrays to hold the E field
      tex,tey,tez = zeros((3,len(x)),'d')

    else:
      isort = None
      nperchild = [len(x)]
      tex,tey,tez = ex,ey,ez

    # --- For each block, pass to it the particles in it's domain.
    i = 0
    for block,n in zip(self.root.listofblocks,nperchild):
      MultiGrid.fetchefrompositions(block,x[i:i+n],y[i:i+n],z[i:i+n],
                                          tex[i:i+n],tey[i:i+n],tez[i:i+n])
      i = i + n

    # --- Now, put the E fields back into the original arrays, unsorting
    # --- the data
    if isort is not None:
      n = len(x)
      putsortedefield(len(tex),isort,tex,tey,tez,ex[:n],ey[:n],ez[:n])

  def fetchphifrompositions(self,x,y,z,phi):
    """
Fetches the potential, given a list of positions
    """
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's. It is assumed at first to be from the local
      # --- domain, and is only set to one of the childs domains where
      # --- childdomains is positive (which does not include any guard cells).
      ichild = zeros(len(x))
      add(ichild,self.blocknumber,ichild)
      getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                            self.nx,self.ny,self.nz,self.childdomains,
                            self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                            self.zmmin,self.zmmax,
                            top.zgridprv,self.l2symtry,self.l4symtry)

      for block in [self]+self.children:
        # --- Get list of particles within the ith domain
        # --- Note that when block==self, the particles selected are the
        # --- ones from this instances domain.
        ii = nonzero(ichild==block.blocknumber)
        if len(ii) == 0: continue
        # --- Get positions of those particles.
        xc = take(x,ii)
        yc = take(y,ii)
        zc = take(z,ii)
        # --- Create temporary arrays to hold the potential
        tphi = zeros(len(xc),'d')
        # --- Now get the field
        if block == self:
          MultiGrid.fetchphifrompositions(self,xc,yc,zc,tphi)
        else:
          block.fetchphifrompositions(xc,yc,zc,tphi)
        # --- Put the potential into the passed in arrays
        put(phi,ii,tphi)

    else:

      # --- Get phi from this domain
      MultiGrid.fetchphifrompositions(self,x,y,z,phi)

  def setphipforparticles(self,*args):
    for block in self.listofblocks:
      SubcycledPoissonSolver.setphipforparticles(block,*args)

  def sortbyichildgetisort(self,ichild,x,y,z):
    xout,yout,zout = zeros((3,len(x)),'d')
    isort = zeros(len(x))
    nperchild = zeros(self.root.totalnumberofblocks)
    sortparticlesbyindexgetisort(len(x),ichild,x,y,z,
                                 self.root.totalnumberofblocks,
                                 xout,yout,zout,isort,nperchild)
    return xout,yout,zout,isort,nperchild

  def getichild_positiveonly(self,x,y,z,ichild):
    """
Gathers the ichild for the fetche_allsort.
    """
    # --- This must wait until all of the parents have have set ichild
    # --- so that the value in the children takes precedence.
    if not self.islastcall(): return
    if len(x) == 0: return
    if len(self.children) > 0:

      # --- Find out whether the particles are in the local domain or one of
      # --- the children's.
      getichildpositiveonly(self.blocknumber,len(x),x,y,z,ichild,
                            self.nx,self.ny,self.nz,self.childdomains,
                            self.xmmin,self.xmmax,self.ymmin,self.ymmax,
                            self.zmmin,self.zmmax,top.zgridprv,
                            self.l2symtry,self.l4symtry)
      for child in self.children:
        child.getichild_positiveonly(x,y,z,ichild)

  #--------------------------------------------------------------------------
  # --- Utility methods
  #--------------------------------------------------------------------------

  def getblockfromnumber(self,number):
    return self.root.listofblocks[number]

  def islastcall(self):
    "Returns true when last parent has called"
    try:                   self.ncallsfromparents
    except AttributeError: self.ncallsfromparents = 0
    self.ncallsfromparents += 1
    if self.ncallsfromparents < len(self.parents): return 0
    self.ncallsfromparents = 0
    return 1

  def isfirstcall(self):
    "Returns true when first parent has called"
    try:                   self.ncallsfromparents
    except AttributeError: self.ncallsfromparents = 0
    self.ncallsfromparents += 1
    if self.ncallsfromparents > 1:
      if self.ncallsfromparents == len(self.parents):
        self.ncallsfromparents = 0
      return 0
    # --- This extra check is needed in case there is one or no parent.
    if self.ncallsfromparents >= len(self.parents):
      self.ncallsfromparents = 0
    return 1

  def setname(self,name='c',ichild=None):
    if not self.isfirstcall(): return
    import __main__
    if ichild is not None:
      name = name + '%d'%ichild
      __main__.__dict__[name] = self
    self.mainname = name
    for child,ichild in zip(self.children,range(1,1+len(self.children))):
      child.setname(name,ichild)

  def getmem(self):
    if not self.isfirstcall(): return
    memtot = product(self.dims + 1)
    for child in self.children:
      memtot = memtot + child.getmem()
    return memtot
      
  def getwarrayforrhop(self,r):
    # --- Create weight array needed for rhop deposition.
    # --- Is linear falloff in the weights correct for r > 2?
    wi = [0,0,0]
    for i in range(3):
      wi[i] = r[i] - abs(1.*iota(-r[i]+1,+r[i]-1))
      wi[i] = wi[i]/sum(wi[i])
    # --- Expand into 3 dimensions
    w = outerproduct(wi[0],outerproduct(wi[1],wi[2]))
    w.shape = (2*r[0]-1,2*r[1]-1,2*r[2]-1)
    result = fzeros(2*r-1,'d')
    result[...] = w
    return result

  def getphip(self,lower=None,upper=None,**kw):
    if len(kw) > 0: return SubcycledPoissonSolver.getphip(self,**kw)
#   if lower is None: lower = self.fulllower - array([0,0,1])
#   if upper is None: upper = self.fullupper + array([0,0,1])
    # --- Note that this takes into account the guard cells in z.
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1 
    iz1 = iz1 + 1
    iz2 = iz2 + 1
    return self.phip[ix1:ix2,iy1:iy2,iz1:iz2,...]
  def getrhop(self,lower=None,upper=None,r=[1,1,1],**kw):
    if len(kw) > 0: return SubcycledPoissonSolver.getrhop(self,**kw)
#   if lower is None: lower = self.fulllower
#   if upper is None: upper = self.fullupper
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.rhop[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2],...]
  def getphi(self,lower,upper):
    # --- Note that this takes into account the guard cells in z.
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1 
    iz1 = iz1 + 1
    iz2 = iz2 + 1
    return self.phi[ix1:ix2,iy1:iy2,iz1:iz2]
  def getrho(self,lower,upper,r=[1,1,1]):
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    return self.rho[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getselfe(self,lower=None,upper=None,comp=slice(None),r=[1,1,1]):
    if lower is None: lower = self.lower
    if upper is None: upper = self.upper
    if type(comp) == StringType:
      comp = ['x','y','z'].index(comp)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + 1
    selfe = MultiGrid.getselfe(self,recalculate=0)
    return selfe[comp,ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]
  def getchilddomains(self,lower,upper,upperedge=0):
    if self.childdomains is None:
      #self.childdomains = fzeros(1+self.dims)  + self.blocknumber
      self.childdomains = fzeros(1+self.dims)
      add(self.childdomains,self.blocknumber,self.childdomains)
    ix1,iy1,iz1 = lower - self.fulllower
    ix2,iy2,iz2 = upper - self.fulllower + upperedge
    return self.childdomains[ix1:ix2,iy1:iy2,iz1:iz2]
  def getlocalarray(self,array,lower,upper,r=[1,1,1],fulllower=None,
                    upperedge=1):
    if fulllower is None: fulllower = self.fulllower
    ix1,iy1,iz1 = lower - fulllower
    ix2,iy2,iz2 = upper - fulllower + upperedge
    return array[ix1:ix2:r[0],iy1:iy2:r[1],iz1:iz2:r[2]]

  def setmgtol(self,mgtol=None):
    """
Sets the convergence tolerance for all blocks. If mgtol is not given, it uses
f3d.mgtol.
    """
    if mgtol is None: mgtol = f3d.mgtol
    self.mgtol = mgtol
    for child in self.children:
      child.setmgtol(mgtol)

  def setmgmaxiters(self,mgmaxiters=None):
    """
Sets the maximum number of iterations for all blocks. If mgmaxiters is
not given, it uses f3d.mgmaxiters.
    """
    if mgmaxiters is None: mgmaxiters = f3d.mgmaxiters
    self.mgmaxiters = mgmaxiters
    for child in self.children:
      child.setmgmaxiters(mgmaxiters)

  def find_mgparam(self,lsavephi=false,resetpasses=0):
    for block in self.listofblocks:
      print "Finding mgparam for block number ",block.blocknumber,me
      # --- Temporarily remove the children so the solve is only done
      # --- on this patch
      childrensave = block.children
      block.children = []
      MultiGrid.find_mgparam(block,lsavephi=lsavephi,resetpasses=resetpasses)
      # --- Restore the children
      block.children = childrensave

  def arraysliceoperation(self,ip,idim,arraystring,op,opnd,null,comp=None):
    """
Applies the operator to the array at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    # --- Each block only needs to check once
    # --- XXX This call breaks something
    #if not self.islastcall(): return null
    # --- Don't do anything if the ip is outside the block
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return null
    # --- Get the appropriate slice of phi and the childdomains array
    ii = [slice(None),slice(None),slice(None)]
    ii[idim] = ip - self.fulllower[idim]
    ix,iy,iz = ii
    if arraystring == 'phi': getarray = self.getphi
    elif arraystring == 'rho': getarray = self.getrho
    elif arraystring == 'selfe': getarray = self.getselfe
    if comp is None: array = getarray(self.fulllower,self.fullupper)
    else:            array = getarray(self.fulllower,self.fullupper,comp)
    if len(self.children) > 0:
      # --- Skip points that don't self doesn't own
      c = self.getchilddomains(self.fulllower,self.fullupper,1)
      array = where(c[ix,iy,iz]==self.blocknumber,array[ix,iy,iz],null)
    else:
      # --- The transpose is taken so that the array is in C ordering so
      # --- the opnd will be faster.
      array = transpose(array[ix,iy,iz])
    # --- Find the max of self's and the children's phi
    result = opnd(array)
    for child in self.children:
      ipc = ip*child.refinement[idim]
      cresult = child.arraysliceoperation(ipc,idim,arraystring,op,opnd,null,
                                          comp)
      result = op(result,cresult)
    return result

# ---------------------------------------------------------------------------
# --- Routines used for plotting
  def getphislicemin(self,ip,idim):
    """
Finds the minimum value of phi at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'phi',min,minnd,+largepos)

  def getphislicemax(self,ip,idim):
    """
Finds the maximum value of phi at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'phi',max,maxnd,-largepos)

  def getrhoslicemin(self,ip,idim):
    """
Finds the minimum value of rho at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'rho',min,minnd,+largepos)

  def getrhoslicemax(self,ip,idim):
    """
Finds the maximum value of rho at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'rho',max,maxnd,-largepos)

  def getselfeslicemin(self,ip,idim,comp):
    """
Finds the minimum value of selfe at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'selfe',min,minnd,
                                    +largepos,comp)

  def getselfeslicemax(self,ip,idim,comp):
    """
Finds the maximum value of selfe at the specified plane. The blocks only
contribute within their domains of ownership.
    """
    return self.arraysliceoperation(ip,idim,'selfe',max,maxnd,
                                    -largepos,comp)

  #--------------------------------------------------------------------------
  # --- The following are used for plotting.
  #--------------------------------------------------------------------------

  def genericpf(self,kw,idim,pffunc,ip=None):
    """
Generic plotting routine. This plots only the local domain. Domains of the
children are also skipped, but the same call is made for them so they will
be plotted.
    """
    # --- Wait until all parents have called so that the child's domain
    # --- if not overlapped by a parent. This only affects the cellaray plots.
    # --- This also avoids the child plotting multiple times.
    if not self.islastcall(): return

    # --- Get the plane to be plotted
    if ip is None:
      ip = kw.get(('ix','iy','iz')[idim],None)
      if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    else:
      ip = ip*self.refinement[idim]
      kw[('ix','iy','iz')[idim]] = ip - self.fulllower[idim]

    # --- Set the values of cmin and cmax for all levels. This must be
    # --- done by the root level.
    if self is self.root:
      cmin = kw.get('cmin',None)
      cmax = kw.get('cmax',None)

      if kw.get('plotselfe',0):
        comp = kw.get('comp',2)
        if cmin is None: cmin = self.getselfeslicemin(ip,idim,comp)
        if cmax is None: cmax = self.getselfeslicemax(ip,idim,comp)
      elif kw.get('plotrho',0):
        if cmin is None: cmin = self.getrhoslicemin(ip,idim)
        if cmax is None: cmax = self.getrhoslicemax(ip,idim)
      else:
        if cmin is None: cmin = self.getphislicemin(ip,idim)
        if cmax is None: cmax = self.getphislicemax(ip,idim)

      cmin = globalmin(cmin)
      cmax = globalmax(cmax)
      kw['cmin'] = cmin
      kw['cmax'] = cmax
      kw['local'] = 1

      accumulateplotlists()

    try:
      # --- Only make the plot if the plane is included in the domain.
      # --- Even if there is no overlap, the children must be called since
      # --- they may overlap (in the domain of a different parent).
      # --- Note that the full extent is used so that the childdomains slice
      # --- will have the same shape as ireg.
      if self.fulllower[idim] <= ip and ip <= self.fullupper[idim]:
        if self.childdomains is not None:
          # --- Create the ireg array, which will be set to zero in the domain
          # --- of the children.
          ss = list(shape(self.rho))
          del ss[idim]
          ireg = zeros(ss)
          ii = [slice(-1),slice(-1),slice(-1)]
          ii[idim] = ip - self.fulllower[idim]
          ix,iy,iz = ii
          ireg[1:,1:]=equal(self.childdomains[ix,iy,iz],self.blocknumber)
          if idim != 2: ireg = transpose(ireg)
          kw['ireg'] = ireg
        else:
          kw['ireg'] = None
        MultiGrid.genericpf(self,kw,pffunc)
        kw['titles'] = 0
        kw['lcolorbar'] = 0
  
      for child in self.children:
        child.genericpf(kw,idim,pffunc,ip)

    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def pfxy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,2,pfxy)
  def pfzx(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,1,pfzx)
  def pfzy(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,0,pfzy)
  def pfxyg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,2,pfxyg)
  def pfzxg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,1,pfzxg)
  def pfzyg(self,kwdict=None,**kw):
    if kwdict is None: kwdict = {}
    kwdict.update(kw)
    self.genericpf(kwdict,0,pfzyg)


  def plphiz(self,ix=None,iy=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.phi[ix-self.fulllower[0],iy-self.fulllower[1],1:-1],self.zmesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plphiz(ix*child.refinement[0],iy*child.refinement[1],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plphix(self,iy=None,iz=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.phi[:,iy-self.fulllower[1],iz-self.fulllower[2]+1],self.xmesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plphix(iy*child.refinement[1],iz*child.refinement[2],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plphiy(self,ix=None,iz=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.phi[ix-self.fulllower[0],:,iz-self.fulllower[2]+1],self.ymesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plphiy(ix*child.refinement[0],iz*child.refinement[2],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plrhoz(self,ix=None,iy=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if ix < self.fulllower[0]: return
    if iy < self.fulllower[1]: return
    if ix > self.fullupper[0]: return
    if iy > self.fullupper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.rho[ix-self.fulllower[0],iy-self.fulllower[1],:],self.zmesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plrhoz(ix*child.refinement[0],iy*child.refinement[1],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plrhox(self,iy=None,iz=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if iy < self.fulllower[1]: return
    if iz < self.fulllower[2]: return
    if iy > self.fullupper[1]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.rho[:,iy-self.fulllower[1],iz-self.fulllower[2]],self.xmesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plrhox(iy*child.refinement[1],iz*child.refinement[2],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plrhoy(self,ix=None,iz=None,colors=None,selfonly=0):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if ix < self.fulllower[0]: return
    if iz < self.fulllower[2]: return
    if ix > self.fullupper[0]: return
    if iz > self.fullupper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.rho[ix-self.fulllower[0],:,iz-self.fulllower[2]],self.ymesh,
          color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plrhoy(ix*child.refinement[0],iz*child.refinement[2],colors=colors)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plselfez(self,comp=2,ix=None,iy=None,colors=None,selfonly=0,withguard=1):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      iz = slice(None)
    else:
      lower,upper = self.lower,self.upper
      iz = slice(self.lower[2] - self.fulllower[2],
                 self.upper[2] - self.fulllower[2] + 1)
    if ix < lower[0]: return
    if iy < lower[1]: return
    if ix > upper[0]: return
    if iy > upper[1]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.selfe[comp,ix-self.fulllower[0],iy-self.fulllower[1],iz],
          self.zmesh[iz],color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plselfez(comp,ix*child.refinement[0],iy*child.refinement[1],
                         colors=colors,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plselfex(self,comp=2,iy=None,iz=None,colors=None,selfonly=0,withguard=1):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      ix = slice(None)
    else:
      lower,upper = self.lower,self.upper
      ix = slice(self.lower[0] - self.fulllower[0],
                 self.upper[0] - self.fulllower[0] + 1)
    if iy < lower[1]: return
    if iz < lower[2]: return
    if iy > upper[1]: return
    if iz > upper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.selfe[comp,ix,iy-self.fulllower[1],iz-self.fulllower[2]],
                     self.xmesh[ix],color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plselfex(comp,iy*child.refinement[1],iz*child.refinement[2],
                         colors=colors,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def plselfey(self,comp=2,ix=None,iz=None,colors=None,selfonly=0,withguard=1):
    if colors is None: colors = color
    elif not operator.isSequenceType(colors): colors = list([colors])
    if withguard:
      lower,upper = self.fulllower,self.fullupper
      iy = slice(None)
    else:
      lower,upper = self.lower,self.upper
      iy = slice(self.lower[1] - self.fulllower[1],
                 self.upper[1] - self.fulllower[1] + 1)
    if ix < lower[0]: return
    if iz < lower[2]: return
    if ix > upper[0]: return
    if iz > upper[2]: return
    if self is self.root: accumulateplotlists()
    try:
      plg(self.selfe[comp,ix-self.fulllower[0],iy,iz-self.fulllower[2]],
          self.ymesh[iy],color=colors[self.blocknumber%len(colors)])
      if not selfonly:
        for child in self.children:
          child.plselfey(comp,ix*child.refinement[0],iz*child.refinement[2],
                         colors=colors,withguard=withguard)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def drawbox(self,ip=None,idim=2,withguards=1,color=[],selfonly=0):
    if len(color)==0: color=['red', 'green', 'blue', 'cyan', 'magenta','yellow']
    if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    if withguards:
      i01 = self.mins[ii[0]]
      i02 = self.maxs[ii[0]]
      i11 = self.mins[ii[1]]
      i12 = self.maxs[ii[1]]
    else:
      i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
      i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
      i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
      i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
    yy = [i01,i01,i02,i02,i01]
    xx = [i11,i12,i12,i11,i11]
    if idim==2:
      yy,xx = xx,yy
    else:
      xx = array(xx) + top.zbeam
    if self is self.root: accumulateplotlists()
    try:
      plg(yy,xx,color=color[0])
      if not selfonly:
        for child in self.children:
          child.drawbox(ip=ip*child.refinement[idim],idim=idim,
                        withguards=withguards,color=color[1:])
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def drawboxzy(self,ix=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=ix,idim=0,withguards=withguards,color=color,
                 selfonly=selfonly)
  def drawboxzx(self,iy=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=iy,idim=1,withguards=withguards,color=color,
                 selfonly=selfonly)
  def drawboxxy(self,iz=None,withguards=1,color=[],selfonly=0):
    self.drawbox(ip=iz,idim=2,withguards=withguards,color=color,
                 selfonly=selfonly)

  def drawfilledbox(self,ip=None,idim=2,withguards=1,ibox=None,selfonly=0):
    if ip is None: ip = nint(-self.mins[idim]/self.deltas[idim])
    if ip < self.fulllower[idim] or ip > self.fullupper[idim]: return
    ii = [0,1,2]
    del ii[idim]
    if withguards:
      i01 = self.mins[ii[0]]
      i02 = self.maxs[ii[0]]
      i11 = self.mins[ii[1]]
      i12 = self.maxs[ii[1]]
    else:
      i01 = self.root.mins[ii[0]] + self.lower[ii[0]]*self.deltas[ii[0]]
      i02 = self.root.mins[ii[0]] + self.upper[ii[0]]*self.deltas[ii[0]]
      i11 = self.root.mins[ii[1]] + self.lower[ii[1]]*self.deltas[ii[1]]
      i12 = self.root.mins[ii[1]] + self.upper[ii[1]]*self.deltas[ii[1]]
    xx = [i01,i01,i02,i02,i01]
    yy = [i11,i12,i12,i11,i11]
    if idim==2: xx,yy = yy,xx
    if ibox is None: ibox = ones(1,'b')
    else:            ibox = (ibox+1).astype('b')
    if self is self.root: accumulateplotlists()
    try:
      plfp(ibox,yy,xx,[5])
      if not selfonly:
        for child in self.children:
          child.drawfilledbox(ip=ip*child.refinement[idim],idim=idim,
                              withguards=withguards,ibox=ibox)
    finally:
      if self is self.root: plotlistofthings(lturnofflist=1)

  def createdxobject(self,kwdict={},**kw):
    """
Create DX object drawing the object.
  - withguards=1: when true, the guard cells are included in the bounding box
    """
    kw.update(kwdict)
    withguards = kw.get('withguards',1)
    xmin,xmax = self.xmmin,self.xmmax
    ymin,ymax = self.ymmin,self.ymmax
    zmin,zmax = self.zmmin,self.zmmax
    if not withguards:
      ng = self.nguard*self.refinement
      xmin,xmax = xmin+ng[0]*self.dx, xmax-ng[0]*self.dx
      ymin,ymax = ymin+ng[1]*self.dy, ymax-ng[1]*self.dy
      zmin,zmax = zmin+ng[2]*self.dz, zmax-ng[2]*self.dz
    dxlist = [viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax)]
    for child in self.children:
      dxlist.append(child.getdxobject(kwdict=kw))
    self.dxobject = DXCollection(*dxlist)
















  #===========================================================================
  def solve2down(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    self.phisave[:,:,:] = self.phi
    cond_potmg(self.conductors.interior,
               self.nx,self.ny,self.nz,self.phisave,0,false,
               2,true)
    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phisave,self.rhosave,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,true,self.lcndbndy,self.icndbndy,self.conductors)
    self.rho[:,:,:] = self.res[:,:,1:-1]
    self.phi[:,:,:] = 0.
    print 1,self.res[10,10,10]

    for child in self.children:
      child.solve2down()

    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      restrict3d(sx2-sx1,sy2-sy1,sz2-sz1,pz2-pz1,sz2-sz1,
                 self.res[sx1:sx2+1,sy1:sy2+1,sz1:sz2+3],
                 parent.res[px1:px2+1,py1:py2+1,pz1:pz2+3],
                 self.boundxy,
                 self.bound0,self.boundnz,self.bound0,self.boundnz,
                 0,0,self.l2symtry,self.l4symtry)

  def solve2up(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      expand3d(px2-px1,py2-py1,pz2-pz1,sz1-sz2,pz2-pz1,
               parent.phi[px1:px2+1,py1:py2+1,pz1:pz2+3],
               self.phi[sx1:sx2+1,sy1:sy2+1,sz1:sz2+3],
               self.boundxy,self.bound0,self.boundnz,0,0)

    #   for i in range(self.uppasses):
    #     self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
    #                    self.phi,self.rho,self.rstar,
    #                    dxsqi,dysqi,dzsqi,self.linbend,
    #                    self.l2symtry,self.l4symtry,self.bendx,
    #                    self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
    #                    self.lcndbndy,self.icndbndy,self.conductors)

    childrenserror = 0.
    for child in self.children:
      childerror = child.solve2up()
      childrenserror = max(childrenserror,childerror)

    add(self.phi,self.phisave,self.phi)
    print 2,self.phi[10,10,10]

    # --- When using residual correction form, the other planes do need
    # --- to be set when using other than Dirichlet boundaries since
    # --- those planes are only set with the error of phi.
    if self.bound0  == 1: self.phi[:,:,0] = self.phi[:,:,2]
    if self.boundnz == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
    if self.bound0  == 2: self.phi[:,:,0] = self.phi[:,:,-3]
    if self.boundnz == 2: self.phi[:,:,-1] = self.phi[:,:,2]

    # --- Calculate the change in phi.
    subtract(self.phisave,self.phi,self.phisave)
    absolute(self.phisave,self.phisave)
    self.mgerror = MA.maximum(self.phisave)
    print self.mgerror,childrenserror
    print 'err = ',self.mgerror
    return max(childrenserror,self.mgerror)

  #===========================================================================
  def solve2init(self):
    # --- Create temp arrays
    self.phisave = fzeros(shape(self.phi),'d')
    self.bendx = fzeros(((self.nx+1)*(self.ny+1)),'d')

    # --- Initialize temporaries
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    checkconductors(self.nx,self.ny,self.nz,self.nzfull,
                    self.dx,self.dy,self.dz,self.conductors,
                    top.my_index,top.nslaves,top.izfsslave,top.nzfsslave)

    # --- Preset rho to increase performance (reducing the number of
    # --- multiplies in the main SOR sweep loop).
    if not self.linbend:
      # --- Do the operation in place (to avoid temp arrays)
      multiply(self.rho,reps0c,self.rho)
    else:
      raise "Bends not yet supported"

    # --- Since using residual correction form, need to save the original rho.
    self.rhosave = self.rho + 0.
    self.res = fzeros(shape(self.phi),'d')

    for child in self.children:
      child.solve2init()

  #===========================================================================
  def solve2(self,iwhich=0):
    # --- No initialization needed
    if iwhich == 1: return

    self.solve2init()

    # --- Initialize temporaries
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    reps0c = self.mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
    rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

    # --- Main multigrid v-cycle loop. Calculate error each iteration since
    # --- very few iterations are done.
    self.mgiters = 0
    self.mgerror = 2.*self.mgtol + 1.
    while (self.mgerror > self.mgtol and self.mgiters < self.mgmaxiters):
      self.mgiters = self.mgiters + 1
 
      self.solve2down()

      # --- Do one vcycle.
      self.vcycle(0,self.nx,self.ny,self.nz,self.nzfull,
                  self.dx,self.dy,self.dz,self.phi,self.rho,
                  self.rstar,self.linbend,self.l2symtry,self.l4symtry,
                  self.bendx,
                  self.boundxy,self.bound0,self.boundnz,
                  self.mgparam,self.mgform,self.mgmaxlevels,
                  self.downpasses,self.uppasses,self.lcndbndy,
                  self.icndbndy,self.conductors)

      self.mgerror = self.solve2up()

      #else
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      # mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
      #                my_index,nslaves,izfsslave,nzfsslave,
      #                whosendingleft,izsendingleft,
      #                whosendingright,izsendingright)
      #endif

    # --- For Dirichlet boundary conditions, copy data into guard planes
    # --- For other boundary conditions, the guard planes are used during
    # --- the solve are so are already set.
    if (self.bound0 == 0): self.phi[:,:,0] = self.phi[:,:,1]
    if (self.boundnz == 0): self.phi[:,:,-1] = self.phi[:,:,-2]

    # --- Make a print out.
    if (self.mgerror > self.mgtol):
      print "MultiGrid: Maximum number of iterations reached"
    print ("MultiGrid: Error converged to %11.3e in %4d v-cycles"%
           (self.mgerror,self.mgiters))

    # --- If using residual correction form, restore saved rho
    self.rho[:,:,:] = self.rhosave

    # --- Restore rho
    if (not self.linbend):
      multiply(self.rho,1./reps0c,self.rho)


  #===========================================================================
  def solve2down1(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    self.phisave[:,:,:] = self.phi
    cond_potmg(self.conductors.interior,
               self.nx,self.ny,self.nz,self.phisave,0,false,
               2,true)
    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phisave,self.rhosave,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,true,self.lcndbndy,self.icndbndy,self.conductors)
    self.rho[:,:,:] = self.res[:,:,1:-1]
    self.phi[:,:,:] = 0.

    for i in range(self.downpasses):
      self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
                     self.phi,self.rho,self.rstar,
                     dxsqi,dysqi,dzsqi,self.linbend,
                     self.l2symtry,self.l4symtry,self.bendx,
                     self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
                     self.lcndbndy,self.icndbndy,self.conductors)

    residual(self.nx,self.ny,self.nz,self.nzfull,dxsqi,dysqi,dzsqi,
             self.phi,self.rho,self.res,
             0,self.bound0,self.boundnz,self.boundxy,
             self.l2symtry,self.l4symtry,
             self.mgparam,2,false,
             self.lcndbndy,self.icndbndy,self.conductors)

    for child in self.children:
      child.solve2down()

    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      restrict3d(sx2-sx1,sy2-sy1,sz2-sz1,pz2-pz1,sz2-sz1,
                 self.res[sx1:sx2+1,sy1:sy2+1,sz1:sz2+1],
                 parent.rho[px1:px2+1,py1:py2+1,pz1:pz2+1],
                 self.boundxy,
                 self.bound0,self.boundnz,self.bound0,self.boundnz,
                 0,0,self.l2symtry,self.l4symtry)

  def solve2up1(self):
    dxsqi  = 1./self.dx**2
    dysqi  = 1./self.dy**2
    dzsqi  = 1./self.dz**2
    for parentnumber in self.parents:
      parent = self.getblockfromnumber(parentnumber)
      s1 = maximum(self.fulllower,parent.fulllower*self.refinement)
      s2 = minimum(self.fullupper,parent.fullupper*self.refinement)
      sx1,sy1,sz1 = s1 - self.fulllower
      sx2,sy2,sz2 = s2 - self.fulllower
      px1,py1,pz1 = s1/self.refinement - parent.fulllower
      px2,py2,pz2 = s2/self.refinement - parent.fulllower
      expand3d(px2-px1,py2-py1,pz2-pz1,sz2-sz2,pz2-pz1,
               parent.phi[px1:px2+1,py1:py2+1,pz1:pz2+1],
               self.phi[sx1:sx2+1,sy1:sy2+1,sz1:sz2+1],
               self.boundxy,self.bound0,self.boundnz,0,0)

    for i in range(self.uppasses):
      self.sorpass3d(0,self.nx,self.ny,self.nz,self.nzfull,
                     self.phi,self.rho,self.rstar,
                     dxsqi,dysqi,dzsqi,self.linbend,
                     self.l2symtry,self.l4symtry,self.bendx,
                     self.bound0,self.boundnz,self.boundxy,self.mgparam,2,
                     self.lcndbndy,self.icndbndy,self.conductors)

    childrenserror = 0.
    for child in self.children:
      childerror = child.solve2up()
      childrenserror = max(childrenserror,childerror)

    add(self.phi,self.phisave,self.phi)

    # --- When using residual correction form, the other planes do need
    # --- to be set when using other than Dirichlet boundaries since
    # --- those planes are only set with the error of phi.
    if self.bound0  == 1: self.phi[:,:,0] = self.phi[:,:,2]
    if self.boundnz == 1: self.phi[:,:,-1] = self.phi[:,:,-3]
    if self.bound0  == 2: self.phi[:,:,0] = self.phi[:,:,-3]
    if self.boundnz == 2: self.phi[:,:,-1] = self.phi[:,:,2]

    # --- Calculate the change in phi.
    subtract(self.phisave,self.phi,self.phisave)
    absolute(self.phisave,self.phisave)
    self.mgerror = MA.maximum(self.phisave)
    print 'err= ',self.mgerror
    return max(childrenserror,self.mgerror)

# --- This can only be done after MRBlock is defined.
try:
  psyco.bind(MRBlock)
except NameError:
  pass

