"""
Controller operations
=====================

These are the functions which allow installing user created functions so that
they are called at various places along the time step.  

For each controller, the following three functions are defined.
 - install___: Installs a function to be called at that specified time
 - uninstall___: Uninstalls the function (so it won't be called anymore)
 - isinstalled___: Checks if the function is installed

The functions all take a function or instance method as an argument. Note that
if an instance method is used, an extra reference to the method's object is saved.

Functions can be called at the following times:
 - :py:func:`aftergenerate <installaftergenerate>`: immediately after the generate is complete
 - :py:func:`beforefs <installbeforefs>`: before the field solve
 - :py:func:`afterfs <installafterfs>`: after the field solve
 - :py:func:`beforeloadrho <installbeforeloadrho>`: before the rho is deposited, at the beginning of loadrho
 - :py:func:`afterloadrho <installafterloadrho>`: after the rho is deposited, at the end of loadrho
 - :py:func:`othereuser <installothereuser>`: during execution of electric fields gathering
 - :py:func:`beforestep <installbeforestep>`: before the time step
 - :py:func:`afterstep <installafterstep>`: after the time step
 - :py:func:`beforescraper <installbeforescraper>`: just before the particle boundary conditions are applied
 - :py:func:`particlescraper <installparticlescraper>`: just after the particle boundary conditions are applied
                    but before lost particles are processed
 - :py:func:`afterscraper <installafterscraper>`: just after the particle boundary conditions are applied
 - :py:func:`particleloader <installparticleloader>`: at the time that the standard particle loader is called
 - :py:func:`addconductor <installaddconductor>`: at the start of the multigrid solver (to installconductors)
 - :py:func:`beforeplot <installbeforeplot>`: before a plot (actually after a frame advance)
 - :py:func:`afterplot <installafterplot>`: after a plot (acutally just before a frame advance)
 - :py:func:`plseldom <installplseldom>`: during a special time step, when position and velocity are
             synchronized, specified by itplseldom or zzplseldom
 - :py:func:`plalways <installplalways>`: during a special time step, when position and velocity are
             synchronized, specified by itplalways or zzplalways
 - :py:func:`userinjection <installuserinjection>`: called when particles injection happens, after the position
                  advance and before loadrho is called, allowing a user defined
                  particle distribution to be injected each time step
 - :py:func:`userparticlesinjection <installuserparticlesinjection>`: allows directly specifying the particles to be injected

"""
from __future__ import generators
controllers_version = "$Id: controllers.py,v 1.32 2010/12/02 23:13:50 dave Exp $"
def controllersdoc():
  import controllers
  print controllers.__doc__

#from warp import *
import warp
from types import *
import copy
import time

class ControllerFunction:
  """
Class to handle the function lists.

Note that for functions passed in that are methods of a class instance,
a full reference of the instance is saved. This extra reference means
that the object will not actually deleted if the user deletes the
original reference.  This is good since the user does not need to keep
the reference to the object (for example it can be created using a local
variable in a function). It is also good since this allows the installed
method to transfer across a dump and restart. It may be bad if the user
thinks an object was deleted, but it actually isn't since it had (unkown
to the user) installed a method in one of the controllers.

This class also provides what is effectively a picklable function
reference. Though it is not complete, since in some cases, functions
won't be restorable. For example, functions typed in interactively cannot
be restored since the source is not saved anywhere.
  """

  def __init__(self,name=None,lcallonce=0):
    self.funcs = []
    self.time = 0.
    self.timers = {}
    self.name = name
    self.lcallonce = lcallonce

  def __call__(self):
    "Call all of the functions in the list"
    tt = self.callfuncsinlist()
    self.time = self.time + tt
    if self.lcallonce: self.funcs = []

  def __getstate__(self):
    """
The instance is picklable. Only the names of functions are saved. A full
reference to a method's object is saved. The names of functions replace
the funcs attribute in the dictionary returned. Note that nothing
special is needed on a restore since the function names will
automatically be converted back into functions the first time they are
called (so there is no __setstate__). For methods, the name is saved since
instancemethods cannot be pickled.
The ControllerFunctionContainer class below ensures that top level
controllers are restored properly.
    """
    dict = self.__dict__.copy()
    del dict['funcs']
    funcnamelist = []
    for f in self.controllerfuncnames():
      funcnamelist.append(f)
    dict['funcs'] = funcnamelist
    return dict

  def hasfuncsinstalled(self):
    "Checks if there are any functions installed"
    return len(self.funcs) > 0

  def getmethodobject(self,func):
    return func[0]

  def controllerfuncnames(self):
    """Returns the names of the functions in the list, and any methods
       (which are stored in lists)"""
    for f in self.funcs:
      if type(f) == ListType:
        result = f
      elif type(f) == StringType:
        import __main__
        if f in __main__.__dict__:
          result = f
        else:
          continue
      else:
        result = f.__name__
      yield result

  def controllerfunclist(self):
    funclistcopy = copy.copy(self.funcs)
    for f in funclistcopy:
      if type(f) == ListType:
        object = self.getmethodobject(f)
        if object is None:
          self.funcs.remove(f)
          continue
        result = getattr(object,f[1])
      elif type(f) == StringType:
        import __main__
        if f in __main__.__dict__:
          result = __main__.__dict__[f]
          # --- If the function with the name is found, then replace the
          # --- name in the list with the function.
          self.funcs[self.funcs.index(f)] = result
        else:
          continue
      else:
        result = f
      if not callable(result):
        print "\n\nWarning: a controller was found that is not callable."
        print "Only callable objects can be installed."
        print "It is possible that the callable's name has been overwritten"
        print "by something not callable. This can happen during restart"
        print "if a function name had later been used as a variable name."
        if type(f) == StringType:
          print "The name of the controller is ",f
        print "\n\n"
        continue
      yield result

  def installfuncinlist(self,f):
    if type(f) == MethodType:
      # --- If the function is a method of a class instance, then save a full
      # --- reference to that instance and the method name.
      finstance = f.im_self
      fname = f.__name__
      self.funcs.append([finstance,fname])
    else:
      self.funcs.append(f)

  def uninstallfuncinlist(self,f):
    # --- An element by element search is needed
    # --- f can be a function or method object, or a name (string).
    # --- Note that method objects can not be removed by name.
    funclistcopy = copy.copy(self.funcs)
    for func in funclistcopy:
      if f == func:
        self.funcs.remove(f)
        return
      elif type(func) == ListType and type(f) == MethodType:
        object = self.getmethodobject(func)
        if f.im_self is object and f.__name__ == func[1]:
          self.funcs.remove(func)
          return
      elif type(func) == StringType:
        if f.__name__ == func:
          self.funcs.remove(func)
          return
      elif type(f) == StringType:
        if type(func) == StringType: funcname = func
        elif type(func) == ListType: funcname = None
        else:                        funcname = func.__name__
        if f == funcname:
          self.funcs.remove(func)
          return
    raise 'Warning: no such function had been installed'

  def isinstalledfuncinlist(self,f):
    # --- An element by element search is needed
    funclistcopy = copy.copy(self.funcs)
    for func in funclistcopy:
      if f == func:
        return 1
      elif type(func) == ListType and type(f) == MethodType:
        object = self.getmethodobject(func)
        if f.im_self is object and f.__name__ == func[1]:
          return 1
      elif type(func) == StringType:
        if f.__name__ == func:
          return 1
    return 0

  def callfuncsinlist(self,*args,**kw):
    bb = time.time()
    for f in self.controllerfunclist():
      t1 = time.time()
      f(*args,**kw)
      t2 = time.time()
      # --- For the timers, use the function (or method) name as the key.
      # --- This is done since instancemethods cannot be pickled.
      self.timers[f.__name__] = self.timers.get(f.__name__,0.) + (t2 - t1)
    aa = time.time()
    return aa - bb

#=============================================================================

# --- Now create the actual instances.
aftergenerate = ControllerFunction('aftergenerate')
beforefs = ControllerFunction('beforefs')
afterfs = ControllerFunction('afterfs')
beforeloadrho = ControllerFunction('beforeloadrho')
afterloadrho = ControllerFunction('afterloadrho')
othereuser = ControllerFunction('othereuser')
beforescraper = ControllerFunction('beforescraper')
afterscraper = ControllerFunction('afterscraper')
callscraper = ControllerFunction('callscraper')
callparticleloader = ControllerFunction('callparticleloader')
calladdconductor = ControllerFunction('calladdconductor')
callbeforestepfuncs = ControllerFunction('callbeforestepfuncs')
callafterstepfuncs = ControllerFunction('callafterstepfuncs')
callbeforeplotfuncs = ControllerFunction('callbeforeplotfuncs')
callafterplotfuncs = ControllerFunction('callafterplotfuncs')
callplseldomfuncs = ControllerFunction('callplseldomfuncs')
callplalwaysfuncs = ControllerFunction('callplalwaysfuncs')
callafterrestartfuncs = ControllerFunction('callafterrestartfuncs',lcallonce=1)
userinjection = ControllerFunction('userinjection')
generateuserparticlesforinjection = ControllerFunction('generateuserparticlesforinjection')

#=============================================================================
class ControllerFunctionContainer:
  """
This is a somewhat kludgy fix to how to get any saved functions restored.
A single instance of this class is created and this instance is what is save
in a dump. This instance will have a list of the controllers, so the
controllers will be saved, but not as top level python variables.
Upon restoration, this container will go through each of the saved controllers
and reinstall the functions saved therein. This installs the functions in the
original set of controllers created when this module was first imported.
Anything that may have already been installed will therefore be unaffected.
  """
  def __init__(self,clist):
    self.clist = clist
  def __setstate__(self,dict):
    import controllers
    import __main__
    self.__dict__.update(dict)
    for c in self.clist:
      for f in c.funcs:
        if type(f) is StringType:
          # --- Check if f is already in the original list of functions,
          # --- and skip it if it is. Both the function name (f) and the
          # --- actual function in main are checked.
          # --- This will be the case if, for example, the user execs the
          # --- original input file, which sets up some functions, before
          # --- doing the restart.
          origfuncs = controllers.__dict__[c.name].funcs
          try:
            ffunc = __main__.__dict__[f]
          except KeyError:
            ffunc = None
          if (f not in origfuncs and ffunc not in origfuncs):
            controllers.__dict__[c.name].installfuncinlist(f)
        else:
          # --- Otherwise, f is a method, so it can be directly installed.
          # --- A check is still made to ensure it isn't installed twice.
          # --- The check is only needed temporarily until the classes
          # --- are fixed to not resinstall in the getstate.
          ffunc = getattr(f[0],f[1])
          if not controllers.__dict__[c.name].isinstalledfuncinlist(ffunc):
            controllers.__dict__[c.name].installfuncinlist(ffunc)
    # --- The clist is obtained from the original instance in the controllers
    # --- module so that the list contains references to the original
    # --- controller instances. This is needed, since in the next dump,
    # --- this instance will be written out and must contain an updated
    # --- list of controllers.
    self.clist = controllers.__dict__['controllerfunctioncontainer'].clist

  def printtimers(self,tmin=1.,lminmax=0.,ff=None):
    """Prints timings of install functions.
 - tmin=1.: only functions with time greater than tmin will be printed
    """
    if ff is None: ff = warp.sys.stdout
    for c in self.clist:
      for f in c.controllerfunclist():
        fname = f.__name__
        try:
          vlist = warp.array(warp.gather(c.timers[fname]))
        except KeyError:
          # --- If the function fname had never been called, then there
          # --- would be no data in timers for it.
          continue
        if warp.me > 0: continue
        vsum = warp.sum(vlist)
        if vsum <= tmin: continue
        vrms = warp.sqrt(max(0.,warp.ave(vlist**2) - warp.ave(vlist)**2))
        ff.write('%20s %s %10.4f  %10.4f %10.4f'%(c.name,fname,vsum,vsum/warp.npes,vrms))
        if lminmax:
          vmin = min(vlist)
          vmax = max(vlist)
          ff.write('  %10.4f  %10.4f'%(vmin,vmax))
        if warp.top.it > 0:
          ff.write('   %10.4f'%(vsum/warp.npes/(warp.top.it)))
        ff.write('\n')

# --- This is primarily needed by warp.py so that these objects can be removed
# --- from the list of python objects which are not written out.
controllerfunctioncontainer = ControllerFunctionContainer(
                               [aftergenerate,beforefs,afterfs,
                                beforeloadrho,afterloadrho,othereuser,
                                beforescraper,afterscraper,callscraper,
                                callparticleloader,calladdconductor,
                                callbeforestepfuncs,callafterstepfuncs,
                                callbeforeplotfuncs,callafterplotfuncs,
                                callplseldomfuncs,callplalwaysfuncs,
                                callafterrestartfuncs,
                                userinjection,
                                generateuserparticlesforinjection])


#=============================================================================
# ----------------------------------------------------------------------------
def installaftergenerate(f):
  "Adds a function to the list of functions called after a generate"
  aftergenerate.installfuncinlist(f)
def uninstallaftergenerate(f):
  "Removes the function from the list of functions called after a generate"
  aftergenerate.uninstallfuncinlist(f)
def isinstalledaftergenerate(f):
  "Checks if the function is called after a generate"
  return aftergenerate.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  beforefs.installfuncinlist(f)
  warp.w3d.lbeforefs = warp.true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  beforefs.uninstallfuncinlist(f)
  if not beforefs.hasfuncsinstalled(): warp.w3d.lbeforefs = warp.false
def isinstalledbeforefs(f):
  "Checks if the function is called before a field-solve"
  return beforefs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  afterfs.installfuncinlist(f)
  warp.w3d.lafterfs = warp.true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  afterfs.uninstallfuncinlist(f)
  if not afterfs.hasfuncsinstalled(): warp.w3d.lafterfs = warp.false
def isinstalledafterfs(f):
  "Checks if the function is called after a field-solve"
  return afterfs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforeloadrho(f):
  "Adds a function to the list of functions called before a load rho"
  beforeloadrho.installfuncinlist(f)
  warp.w3d.lbeforelr = warp.true
def uninstallbeforeloadrho(f):
  "Removes the function from the list of functions called before a load rho"
  beforeloadrho.uninstallfuncinlist(f)
  if not beforeloadrho.hasfuncsinstalled(): warp.w3d.lbeforelr = warp.false
def isinstalledbeforeloadrho(f):
  "Checks if the function is called before a load rho"
  return beforeloadrho.isinstalledfuncinlist(f)

# --- This are defined for backwards compatibility
installbeforelr = installbeforeloadrho
uninstallbeforelr = uninstallbeforeloadrho
isinstalledbeforelr = isinstalledbeforeloadrho

# ----------------------------------------------------------------------------
def installafterloadrho(f):
  "Adds a function to the list of functions called after a load rho"
  afterloadrho.installfuncinlist(f)
  warp.w3d.lafterloadrho = warp.true
def uninstallafterloadrho(f):
  "Removes the function from the list of functions called after a load rho"
  afterloadrho.uninstallfuncinlist(f)
  if not afterloadrho.hasfuncsinstalled(): warp.w3d.lafterloadrho = warp.false
def isinstalledafterloadrho(f):
  "Checks if the function is called after a load rho"
  return afterloadrho.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installothereuser(f):
  "Adds a function to the list of functions called during the electric fields gathering"
  othereuser.installfuncinlist(f)
  warp.w3d.lothereuser = warp.true
def uninstallothereuser(f):
  "Removes the function from the list of functions called during the electric fields gathering"
  othereuser.uninstallfuncinlist(f)
  if not othereuser.hasfuncsinstalled(): warp.w3d.lothereuser = warp.false
def isinstalledothereuser(f):
  "Checks if the function is called during the electric fields gathering"
  return othereuser.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforescraper(f):
  "Adds a function to the list of functions called before scraping particles"
  beforescraper.installfuncinlist(f)
  warp.w3d.lbeforescraper = warp.true
def uninstallbeforescraper(f):
  "Removes the function from the list of functions called before scraping particles"
  beforescraper.uninstallfuncinlist(f)
  if not beforescraper.hasfuncsinstalled(): warp.w3d.lbeforescraper = warp.false
def isinstalledbeforescraper(f):
  "Checks if the function is called before scraping particles"
  return beforescraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterscraper(f):
  "Adds a function to the list of functions called after scraping particles"
  afterscraper.installfuncinlist(f)
  warp.w3d.lafterscraper = warp.true
def uninstallafterscraper(f):
  "Removes the function from the list of functions called after scraping particles"
  afterscraper.uninstallfuncinlist(f)
  if not afterscraper.hasfuncsinstalled(): warp.w3d.lafterscraper = warp.false
def isinstalledafterscraper(f):
  "Checks if the function is called after scraping particles"
  return afterscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installparticlescraper(f):
  "Adds a function to the list of functions called to scrape particles"
  callscraper.installfuncinlist(f)
  warp.w3d.lcallscraper = warp.true
def uninstallparticlescraper(f):
  "Removes the function from the list of functions called to scrape particles"
  callscraper.uninstallfuncinlist(f)
  if not callscraper.hasfuncsinstalled(): warp.w3d.lcallscraper = warp.false
def isinstalledparticlescraper(f):
  "Checks if the function is called to scrape particles"
  return callscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installparticleloader(f):
  "Adds a function to the list of functions called to load particles"
  callparticleloader.installfuncinlist(f)
  warp.w3d.lcallparticleloader = warp.true
def uninstallparticleloader(f):
  "Removes the function from the list of functions called to load particles"
  callparticleloader.uninstallfuncinlist(f)
  if not callparticleloader.hasfuncsinstalled():
    warp.w3d.lcallparticleloader = warp.false
def isinstalledparticleloader(f):
  "Checks if the function is called to load particles"
  return callparticleloader.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installaddconductor(f):
  "Adds a function to the list of functions called to add conductors"
  calladdconductor.installfuncinlist(f)
  warp.f3d.laddconductor = warp.true
def uninstalladdconductor(f):
  "Removes the function from the list of functions called to add conductors"
  calladdconductor.uninstallfuncinlist(f)
  if not calladdconductor.hasfuncsinstalled(): warp.f3d.laddconductor = warp.false
def isinstalledaddconductor(f):
  "Checks if the function is called to add conductors"
  return calladdconductor.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  callbeforestepfuncs.installfuncinlist(f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  callbeforestepfuncs.uninstallfuncinlist(f)
def isinstalledbeforestep(f):
  "Checks if the function is called before a step"
  return callbeforestepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  callafterstepfuncs.installfuncinlist(f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  callafterstepfuncs.uninstallfuncinlist(f)
def isinstalledafterstep(f):
  "Checks if the function is called after a step"
  return callafterstepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforeplot(f):
  "Adds a function to the list of functions called before a plot"
  beforeplotfuncs.installfuncinlist(f)
def uninstallbeforeplot(f):
  "Removes the function from the list of functions called before a plot"
  beforeplotfuncs.uninstallfuncinlist(f)
def isinstalledbeforeplot(f):
  "Checks if the function is called before a plot"
  return beforeplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterplot(f):
  "Adds a function to the list of functions called after a plot"
  callafterplotfuncs.installfuncinlist(f)
def uninstallafterplot(f):
  "Removes the function from the list of functions called after a plot"
  callafterplotfuncs.uninstallfuncinlist(f)
def isinstalledafterplot(f):
  "Checks if the function is called after a plot"
  return callafterplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installplseldom(f):
  "Adds a function to the list of functions controlled by itplseldom and zzplseldom"
  callplseldomfuncs.installfuncinlist(f)
def uninstallplseldom(f):
  "Removes the function from the list of functions controlled by itplseldom and zzplseldom"
  callplseldomfuncs.uninstallfuncinlist(f)
def isinstalledplseldom(f):
  "Checks if the function is controlled by itplseldom and zzplseldom"
  return callplseldomfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installplalways(f):
  "Adds a function to the list of functions controlled by itplalways and zzplalways"
  callplalwaysfuncs.installfuncinlist(f)
def uninstallplalways(f):
  "Removes the function from the list of functions controlled by itplalways and zzplalways"
  callplalwaysfuncs.uninstallfuncinlist(f)
def isinstalledplalways(f):
  "Checks if the function is controlled by itplalways and zzplalways"
  return callplalwaysfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterrestart(f):
  "Adds a function to the list of functions called immediately after a restart"
  callafterrestartfuncs.installfuncinlist(f)
def uninstallafterrestart(f):
  "Removes the function from the list of functions called immediately after a restart"
  callafterrestartfuncs.uninstallfuncinlist(f)
def isinstalledafterrestart(f):
  "Checks if the function is called immediately after a restart"
  return callafterrestartfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installuserinjection(f):
  """
Adds a user defined function that is to be called when particles
injection happens, after the position advance and before loadrho is
called, allowing a user defined particle distribution to be injected
each time step"""
  userinjection.installfuncinlist(f)
def uninstalluserinjection(f):
  "Removes the function installed by installuserinjection"
  userinjection.uninstallfuncinlist(f)
def isinstalleduserinjection(f):
  "Checks if the function is called when particles injection happens"
  return userinjection.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installuserparticlesinjection(f):
  """
Adds a user defined function that is to be called during injection which
allows the user to specify the distribution of particles on the emitting
surface. For expert use only.
To use, the installed function should set w3d.npgrp to the number of
particles to inject, call gchange("Setpwork3d") to allocate the arrays, and
fill the arrays w3d.xt, yt, uxt, uyt, and uzt with the particle data. The
particles start on the emitting surface and the code will advance them away
from the surface. The function will be called once for each species each time
step, with the variable w3d.inj_js set to the species being injected. Note
that if no particles are to be injected, set w3d.npgrp=0 to avoid injection
of bad particles."""
  warp.w3d.l_inj_user_particles = warp.true
  generateuserparticlesforinjection.installfuncinlist(f)
def uninstalluserparticlesinjection(f):
  "Removes the function installed by installuserparticlesinjection"
  generateuserparticlesforinjection.uninstallfuncinlist(f)
  if not generateuserparticlesforinjection.hasfuncsinstalled():
    warp.w3d.l_inj_user_particles = warp.false
def isinstalleduserparticlesinjection(f):
  "Checks if the function is called during injection"
  return generateuserparticlesforinjection.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def fixcontrollersfromolddump():
  import __main__
  import controllers
  controllernames = ['aftergenerate','beforefs','afterfs','callscraper',
                     'calladdconductor','callbeforestepfuncs',
                     'callafterstepfuncs','callbeforeplotfuncs',
                     'callafterplotfuncs','callplseldomfuncs',
                     'callplalwaysfuncs']
  for cname in controllernames:
    if cname in __main__.__dict__:
      controller = __main__.__dict__[cname]
      if 'funcnamelist' in controller.__dict__:
        controllers.__dict__[controller.name] = controller
        controller.funcs = controller.funcnamelist
        del controller.funcnamelist
    else:
      print "Controller ",cname," not found"

