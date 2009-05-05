"""Utility and convenience routines used in Warp

iota(): returns a array of sequential integers
remark(): same as print
sign(): emulation of sign function
arrayappend(): appends to multi-dimensional array
exists(): checks if a variable exists
ave(): averages an array of numbers
maxnd(): finds max of multi-dimensional array
minnd(): finds min of multi-dimensional array
span(): returns an array of evenly spaced numbers
getnextfilename(): finds next available file name in a numeric sequence
getmesh2d(): returns tuple of 2 2-d arrays holding the coordinates of the
             two dimensions
getmesh3d(): returns tuple of 3 3-d arrays holding the coordinates of the
             three dimensions
averagezdata(): Does local averaging over the first dimension of the input
                array. Can also do down select on the data.
getdatafromtextfile(): Reads in table data from a text file, returning an array
                       that holds that data.
RandomStream: Class the allows creating independent streams of random
              numbers.
"""
from __future__ import generators # needed for yield statement for P2.2
from warp import *
import struct # needed for makefortranordered
import appendablearray

warputils_version = "$Id: warputils.py,v 1.27 2009/05/05 00:55:30 dave Exp $"

def warputilsdoc():
  import warputils
  print warputils.__doc__

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
  if high is None:
    if step > 0:
      return arange(1,low+1,step)
    else:
      return arange(low,0,step)
  else:
    if step > 0:
      return arange(low,high+1,step)
    else:
      return arange(low,high-1,step)

# --- Convenience function to do printing
def remark(s):
  print s

# --- Replicate the sign function with two arguments. If only one is
# --- given, return the value from the numpy sign function.
# --- This should realy be removed.
numpysign = sign
def sign(x,y=None):
  if y is None: return numpysign(x)
  if isinstance(x,ndarray):
    result = where(greater(y,0.),abs(x),-abs(x))
    result = where(equal(y,0.),0.,result)
    return result
  else:
    if y > 0:
      return abs(x)
    elif y < 0:
      return -abs(x)
    else:
      return 0

# --- Convenience function which returns meshes filled with the coordinates.
def getmeshcoordinates(mins,dds,nns):
  """
getmeshcoordinates(mins,dds,nns)
Returns arrays holding the coordinates of the mesh points.
Lenght of list of inputs determines number of dimensions.
  """
  nns = tuple(array(nns) + 1)
  cc = indices(nns,'d')
  for i in xrange(len(mins)): cc[i] = mins[i] + cc[i]*dds[i]
  clist = []
  for i in xrange(len(mins)): clist.append(cc[i])
  return tuple(clist)

def getmesh2d(xmin,dx,nx,ymin,dy,ny):
  """
getmesh2d(xmin,dx,nx,ymin,dy,ny)
Returns 2 2-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin],[dx,dy],[nx,ny])

def getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz):
  """
getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dx,nz)
Returns 3 3-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin,zmin],[dx,dy,dz],[nx,ny,nz])

# --- This function appends a new element to the end of an array.
# --- It is not very efficient since it creates a whole new array each time.
def arrayappend(x,a):
  xshape = list(shape(x))
  if isinstance(a,ndarray):
    pass
  elif isinstance(a,ListType):
    a = array(a)
  else:
    a = array([a])
  ashape = list(shape(a))
  if len(xshape)==1 and len(ashape)==1:
    xshape[0] = xshape[0] + ashape[0]
    y = fzeros(xshape,gettypecode(x))
    y[0:xshape[0]-ashape[0]] = x
    y[xshape[0]-ashape[0]:] = a
  elif len(xshape) == len(ashape)+1 and xshape[:-1] == ashape:
    xshape[-1] = xshape[-1] + 1
    y = fzeros(xshape,gettypecode(x))
    y[...,0:-1] = x
    y[...,-1] = a
  return y

# Convenience function which returns true if variable exists
def exists(x):
  """
Checks whether or not the variable whose name is specified exists in the
main dictionary.
 - x: Name of variable - must be a string.
  """
  import __main__
  if x in __main__.__dict__.keys(): return true
  if x in locals().keys(): return true
  if x in globals().keys(): return true
  return false

# Returns the average of the input array
def ave(x,index=0):
  if len(shape(x)) == 0: return x
  if shape(x)[index] > 0:
    return sum(x,index)/shape(x)[index]
  else:
    return 0.

def averagezdata(qty,navg=0,nlines=100,n1=None,n2=None,istep=None,
                 includezeros=false):
  """
Averages data over local region. It also can down select data in the other
dimension.
  - qty: Data to be smoothed. Can be either a 1-D or 2-D array.
  - navg=0: number of data points to average over
  - nlines=100: number of lines from second dimension to choose.
  - n1=shape(qty)[0]-1:
  - n2=shape(qty)[1]-1:
  - istep=max(1,n2/nlines):
  - includezeros=false: by default, only non-zero data is averaged over.
  """
  if shape(qty)[0] < navg+navg+1: return qty
  if navg == 0 or nlines == 0: return qty
  if len(shape(qty)) == 1:
    fixqty = 1
    qty.shape = (len(qty),1)
  else:
    fixqty = 0
  if not n1: n1 = shape(qty)[0] - 1
  if not n2: n2 = shape(qty)[1] - 1
  if istep is None: istep = max(1,n2/nlines)
  hl = qty[:,::istep] + 0.
  hl[navg,:] = sum(qty[navg-navg:navg+navg+1,::istep])
  nn = 2*navg+1 + zeros(shape(hl),'l')
  if not includezeros:
    nn[navg,:] = sum(where(qty[navg-navg:navg+navg+1,::istep]==0.,0,1),0)
  for j in range(navg+1,n1-navg-1):
    hl[j,:] = hl[j-1,:] + (qty[j+navg,::istep] - qty[j-navg-1,::istep])
    nn[j,:] = nn[j-1,:] + (+ where(qty[j+navg,::istep]==0,0,1)
                           - where(qty[j-navg-1,::istep]==0,0,1))
  nn = where(nn==0,1,nn)
  hl = where(qty[:,::istep]==0.,0.,hl)
  hl[navg:n1-navg-1,:] = hl[navg:n1-navg-1,:]/nn[navg:n1-navg-1,:]
  if fixqty: qty.shape = (shape(qty)[0],)
  if n2 > 1: return hl
  else: return hl[:,0]

# --- Returns the max of the multiarray
def maxnd(x,defaultval=-1.e36):
  """Return the max element of an array of any dimension"""
  if isinstance(x,ndarray) and x.size > 0: return x.max()
  xtemp = reshape(x,tuple([product(array(shape(x)))]))
  if len(xtemp) == 0: return defaultval
  return max(xtemp)
# --- Returns the min of the multiarray
def minnd(x,defaultval=+1.e36):
  """Return the min element of an array of any dimension"""
  if isinstance(x,ndarray) and x.size > 0: return x.min()
  xtemp = reshape(x,tuple([product(array(shape(x)))]))
  if len(xtemp) == 0: return defaultval
  return min(xtemp)
# --- Returns the sum of the multiarray
def sumnd(x,defaultval=0.):
  """Return the total sum of an array of any dimension"""
  if isinstance(x,ndarray) and x.size > 0: return x.sum()
  xtemp = reshape(x,tuple([product(array(shape(x)))]))
  if len(xtemp) == 0: return defaultval
  return sum(xtemp)
# --- Returns the sum of the multiarray
def avend(x,defaultval=0.):
  """Return the average of an array of any dimension"""
  if isinstance(x,ndarray) and x.size > 0: return x.sum()/x.size
  xtemp = reshape(x,tuple([product(array(shape(x)))]))
  if len(xtemp) == 0: return defaultval
  return sum(xtemp)/len(xtemp)

def span(lo, hi, num):
  """Returns an array of num equally spaced numbers starting with lo and
ending with hi.
  """
  return lo + (hi - lo)*arange(num)/(num-1.)

def makefortranordered(x):
  """Given an array, returns the same data but with fortran ordering. 
If the array already has the correct ordering, the array is just
returned as is. Otherwise, a new array is created and the data copied."""
  assert isinstance(x,ndarray),"Input value must be an array."
  # --- Pick any package, since all have the getstrides method
  pkg = packageobject(getcurrpkg())
  # --- An array is fortran ordered if the strides are increasing,
  # --- starting at the elemental size. This means that the first
  # --- index varies the fastest in memory.
  strides = pkg.getstrides(x)
  fordered = 1
  ss = struct.calcsize(gettypecode(x))
  for id in range(len(x.shape)):
    if strides[id] != ss: fordered = 0
    ss = ss*x.shape[id]
  
  if fordered:
    # --- No change needed.
    return x
  else:
    xf = fzeros(x.shape,gettypecode(x))
    xf[...] = x
    return xf

# Gets next available filename with the format 'root.nnn.suffix'.
def getnextfilename(root,suffix):
  dir = string.join(os.listdir('.'))
  i = 0
  name = root+('.%03d.'%i)+suffix
  # --- escape adds a backslash to all non-alphanumeric characters. This is
  # --- needed in case special re search characters are used as part of the
  # --- file name, for example '+' or '.'.
  while re.search(re.escape(name),dir):
    i = i + 1
    name = root+('.%03d.'%i)+suffix
  return name

# --- Another Example profiler - this prints out the function name when the
# --- function is called and when it returns. Indentation is used to signify
# --- the call level.
# --- To use, execute the command sys.setprofile(warpprofile).
def warpprofilesample(frame,event,arg):
  try:    warpprofile.level
  except: warpprofile.level = 0
  if event == 'return': warpprofile.level = warpprofile.level - 1
  print "%s %s %s"%(warpprofile.level*'  ',event,frame.f_code.co_name)
  if event == 'call': warpprofile.level = warpprofile.level + 1

# --- Create enumerate function, which is defined in python2.3 but not 2.2
try:
  enumerate
except:
  def enumerate(ll):
    tt = []
    for i in range(len(ll)):
      tt.append((i,ll[i]))
    return tt

# --- Convenience function to read in data from a text file
def getdatafromtextfileold(filename,nskip=0,dims=[],nquantities=1,dtype='d',
                        fortranordering=1,converter=float,mode='r',get_header=false):
  """
Reads data in from a text file. The data is assumed to be laid out on a
logically Cartesian mesh.
 - filename: must be supplied
 - nskip=0: numbers of lines at the beginning of the file to skip
            e.g. lines with file info and comments
 - dims=[]: must be supplied - the size of each dimension
 - nquantities=1: number of quantities in the file
 - fortranordering=1: when true, the data will be in fortran ordering, where
                      the index that varies that fastest in the file will be
                      the first index. Otherwise use C ordering.
 - converter=float: Function which converts the strings into numbers. This should
                    be the type of result desired. It should only be float or int,
                    unless you know what you are doing.

Here's an example data file called 'testdata':
----------------------------------------------
this line is skipped
1 0.0379
2 0.0583
3 0.0768
4 0.1201
----------------------------------------------
This can be read in with either
dd = getdatafromtextfile('testdata',nskip=1,dims=[2,4])
or
dd = getdatafromtextfile('testdata',nskip=1,dims=[4],nquantities=2)
Both produce an array of shape (2,4) that looks like
>>> print dd
[[ 1.    , 2.    , 3.    , 4.    ,]
 [ 0.0379, 0.0583, 0.0768, 0.1201,]]

  """

  # --- Get total number of data values and make an array big enough to hold
  # --- them.
  ntot = nquantities*product(dims)
  data = zeros(ntot,dtype)
  if get_header:header=[]
  ff = open(filename,mode)

  # --- Skip the number of lines at the top of the file as specified
  for i in range(nskip):
    if get_header:
      header.append(ff.readline())
    else:
      ff.readline()

  # --- Loop over the file, reading in one line at a time.
  # --- Each whole line is put into data at once.
  # --- For the conversion of the strings into numbers, it is faster to use
  # --- float or int rather than eval.
  i = 0
  while i < ntot:
    dataline = map(converter,string.split(ff.readline()))
    nd = len(dataline)
    data[i:i+nd] = dataline
    i = i + nd

  if not fortranordering:
    # --- reverse the order of the dims to prepare the transpose that follows
    dims = list(dims)
    dims.reverse()

  # --- If nquantities, then add another dimension to the data array
  if nquantities > 1:
    dims0 = dims
    dims = [nquantities] + list(dims)
  else:
    dims = list(dims)

  # --- Set array to have proper shape.
  dims.reverse()
  data.shape = tuple(dims)
  if fortranordering:
    data = transpose(data)
  else:
    data = transpose(data,[len(dims0)]+range(len(dims0)))

  if get_header:
    return data,header
  else:
    return data

def getdatafromtextfile(filename,nskip=0,dims=[],dtype='d',fortranordering=1,
                        converter=float,mode='r',get_header=false):
  """
Reads data in from a text file. The data is assumed to be laid out on a
logically Cartesian mesh.
 - filename: must be supplied
 - nskip=0: numbers of lines at the beginning of the file to skip
            e.g. lines with file info and comments
 - dims=[]: must be supplied - the size of each dimension
            The last dimension can be None - it will be calculated
            automatically based on the amount of data in the file.
 - fortranordering=1: when true, the data will be in fortran ordering, where
                      the index that varies that fastest in the file will be
                      the first index. Otherwise use C ordering.
 - converter=float: Function which converts the strings into numbers. This should
                    be the type of result desired. It should only be float or int,
                    unless you know what you are doing.

Here's an example data file called 'testdata':
----------------------------------------------
this line is skipped
1 0.0379
2 0.0583
3 0.0768
4 0.1201
----------------------------------------------
This can be read in with
>>> dd = getdatafromtextfile('testdata',nskip=1,dims=[2,4])
to produce an array of shape (2,4) that looks like
>>> print dd
[[ 1.    , 2.    , 3.    , 4.    ,]
 [ 0.0379, 0.0583, 0.0768, 0.1201,]]

The second dimension of dims can be None, meaning that all of the data
will be read in and the size of that dimension will be determined
automatically.
>>> dd = getdatafromtextfile('testdata',nskip=1,dims=[2,None])
>>> print dd
[[ 1.    , 2.    , 3.    , 4.    ,]
 [ 0.0379, 0.0583, 0.0768, 0.1201,]]

  """
  ff = open(filename,mode)

  # --- Skip the number of lines at the top of the file as specified
  header = []
  for i in range(nskip):
    header.append(ff.readline())

  # --- If last dimension is given, calculate the total amount of data to
  # --- be read in to give an exit condition for the loop below.
  # --- This allows less than the whole file to be read in.
  if dims[-1] is not None:
    ntot = product(dims)
  else:
    ntot = None

  # --- Loop over the file, reading in one line at a time.
  # --- Each whole line is put into data at once.
  # --- For the conversion of the strings into numbers, it is faster to use
  # --- float or int rather than eval.
  data = appendablearray.AppendableArray()
  for line in ff.readlines():
    dataline = map(converter,string.split(line))
    data.append(dataline)
    if ntot is not None and len(data) >= ntot: break

  # --- Get the data out of the appendable array.
  data = data[...]

  # --- If last dimension is None, calculate it from the amount of data.
  if dims[-1] is None:
    ndata = len(data)
    ndims = product(dims[:-1])
    nlast = int(ndata/ndims)
    assert nlast*ndims == ndata,"Amount of data does not conform to dims"
    dims[-1] = nlast

  if fortranordering:
    # --- reverse the order of the dims
    dims.reverse()

  # --- Set array to have proper shape.
  data.shape = tuple(dims)
  if fortranordering:
    data = transpose(data)
  else:
    data = transpose(data,[len(dims)-1]+range(len(dims)-1))

  if get_header:
    return data,header
  else:
    return data

class RandomStream(object):
  """Create an independent stream of random numbers.
>>> r1 = RandomStream()
>>> r2 = RandomStream()
>>> r1.seed([1,2,3])
>>> r2.seed([1,2,3])
>>> print r1.random()
0.609861272287
>>> print r2.random()
0.609861272287
  """
  def __init__(self):
    # --- Save the initial random number state
    self._state = random.get_state()
  def __getattr__(self,name):
    # --- For all attribute access, save the attribute name and return the
    # --- _wrapper method.
    self._savedname = name
    return self._wrapper
  def _wrapper(self,*args,**kw):
    # --- This method calls the function from the random module with the saved
    # --- name, though first it switches to the stream's state, and then
    # --- switches back afterward.
    savedstate = random.get_state()
    random.set_state(self._state)
    result = getattr(random,self._savedname)(*args,**kw)
    self._state = random.get_state()
    random.set_state(savedstate)
    return result

