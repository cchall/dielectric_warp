"""Creates a class for handling extrapolated particle windows, ExtPart.
Type doc(ExtPart) for more help.
Two functions are available for saving the object in a file.
 - dumpExtPart(object,filename)
 - restoreExtPart(object,filename)
"""
from warp import *
from appendablearray import *
import cPickle
import string
extpart_version = "$Id: extpart.py,v 1.33 2004/04/15 18:39:08 dave Exp $"

def extpartdoc():
  import extpart
  print extpart.__doc__

_extforcenorestore = 0
def extforcenorestore():
  global _extforcenorestore
  _extforcenorestore = 1
def extnoforcenorestore():
  global _extforcenorestore
  _extforcenorestore = 0

############################################################################
class ExtPart:
  """This class defines a container to setup and keep track of extropolated
particle data. It can optionally accumulate the data over multiple time steps.
The creator options are:
 - iz: grid location where the extrapolated data is saved.
 - zz: lab location where data is saved.
 - wz: width of lab window
 - nepmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
           otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.
 - name=None: descriptive name for location
 - lautodump=0: when true, after the grid moves beyond the z location,
                automatically dump the data to a file, clear the arrays and
                disable itself. Also must have name set.
 - dumptofile=0: when true, the particle data is always dumped to a file
                 and not saved in memory. Name must be set. Setting this
                 to true implies that the data is accumulated.

One of iz or zz must be specified.

Available methods:
 - setaccumulate(v=1): Turns accumulation on or off. If turned on, all
                       currently saved data is deleted.
 - clear(): Will clear the existing data.
 - disable(): Turns off collecting of data
 - enable(): Turns on collecting of data (only needed after disable)

The follow all take an optional argument to specify species number.
 - getn: Get number of particles
 - gett: Get time at which particle was saved
 - getx, y, ux, uy, uz, vx, vy, vz, xp, yp, r, theta, rp: Get the various
     coordinates or velocity of particles

The following are available plot routines. All take an optional argument to
specify species number. Additional arguments are the same as the 'pp' plotting
routines (such as ppxxp).
 - pxy, pxxp, pyyp, pxpyp, prrp, ptx, pty, ptxp, ptyp, ptux, ptuy, ptuz, ptvx
 - ptvy, ptvz, ptrace
  """

  def __init__(self,iz=-1,zz=0.,wz=None,nepmax=None,laccumulate=0,
               name=None,lautodump=1,dumptofile=0):
    # --- Save input values, getting default values when needed
    assert iz >= 0 or zz is not None,"Either iz or zz must be specified"
    self.iz = iz
    self.zz = zz
    if wz is None: self.wz = w3d.dz
    else:          self.wz = wz
    self.laccumulate = laccumulate
    self.lautodump = lautodump
    self.name = name
    self.dumptofile = dumptofile
    self.dt = top.dt
    if nepmax is None:
      self.nepmax = 10000
      if top.allocated("pnumz") and 0 <= self.getiz() <= top.nzmmnt:
        if top.pnumz[self.getiz()] > 0: self.nepmax = top.pnumz[self.getiz()]*3
    else:
      self.nepmax = nepmax
    # --- Add this new window to the ExtPart group in top
    self.enabled = 0
    self.enable()
    # --- Setup empty arrays for accumulation if laccumulate if true.
    # --- Otherwise, the arrays will just point to the data in ExtPart.
    self.setuparrays(top.ns)

  def getiz(self):
    if self.iz >= 0:
      return self.iz
    else:
      return int((self.zz - top.zmmntmin)*top.dzmi)

  def setuparrays(self,ns,bump=None):
    if self.laccumulate and not self.dumptofile:
      if bump is None: bump = self.nepmax
      self.tep = []
      self.xep = []
      self.yep = []
      self.uxep = []
      self.uyep = []
      self.uzep = []
      for js in range(ns):
        self.tep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.xep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.yep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uxep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uyep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uzep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
    else:
      self.tep = ns*[zeros(0,'d')]
      self.xep = ns*[zeros(0,'d')]
      self.yep = ns*[zeros(0,'d')]
      self.uxep = ns*[zeros(0,'d')]
      self.uyep = ns*[zeros(0,'d')]
      self.uzep = ns*[zeros(0,'d')]

  def addspecies(self):
    if self.laccumulate and not self.dumptofile:
      for js in range(len(self.tep),top.ns):
        bump = self.nepmax
        self.tep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.xep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.yep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uxep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uyep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
        self.uzep.append(AppendableArray(self.nepmax,type='d',autobump=bump))
    else:
      self.tep = top.ns*[zeros(0,'d')]
      self.xep = top.ns*[zeros(0,'d')]
      self.yep = top.ns*[zeros(0,'d')]
      self.uxep = top.ns*[zeros(0,'d')]
      self.uyep = top.ns*[zeros(0,'d')]
      self.uzep = top.ns*[zeros(0,'d')]

  def clear(self):
    self.setuparrays(top.ns)

  def getid(self,safe=0):
    'If safe, then return None is id is not found rather than raising error'
    assert self.enabled,"This window is disabled and there is no associated id"
    for i in range(top.nepwin):
      if top.izepwin[i] == self.iz and self.iz >= 0: return i
      if top.zzepwin[i] == self.zz and self.iz == -1: return i
    if not safe:
      raise "Uh Ooh! Somehow the window was deleted! I can't continue! "+self.titleright(None,None)
    else:
      return None

  def setupid(self):
    top.nepwin = top.nepwin + 1
    if top.nepmax < self.nepmax: top.nepmax = self.nepmax
    err = gchange("ExtPart")
    top.izepwin[-1] = self.iz
    top.zzepwin[-1] = self.zz
    top.wzepwin[-1] = self.wz

  def enable(self):
    # --- Add this window to the list
    # --- Only add this location to the list if it is not already there.
    # --- Note that it is not an error to have more than one instance
    # --- have the same location. For example one could be accumulating
    # --- while another isn't or the widths could be different.
    if self.enabled: return
    self.setupid()
    # --- Set so accumulate method is called after time steps
    installafterstep(self.accumulate)
    self.enabled = 1

  def disable(self):
    if not self.enabled: return
    # --- Set so accumulate method is not called after time steps
    uninstallafterstep(self.accumulate)
    # --- Remove this window from the list. Turn safe on when gettin
    # --- the id, since it may for some reason not be consistent.
    id = self.getid(safe=1)
    if id is not None:
      for i in range(id,top.nepwin-1):
        top.izepwin[i] = top.izepwin[i+1]
        top.zzepwin[i] = top.zzepwin[i+1]
        top.wzepwin[i] = top.wzepwin[i+1]
      top.nepwin = top.nepwin - 1
      gchange("ExtPart")
    self.enabled = 0

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    if not self.enabled: return
    id = self.getid(safe=1)
    if id is None: self.setupid()
    self.restoredata()
    if not isinstalledafterstep(self.accumulate):
      installafterstep(self.accumulate)

  def autodump(self):
    if not self.lautodump or self.name is None: return
    if not self.laccumulate and not self.dumptofile: return
    if self.iz >= 0: return
    if self.zz+self.wz > w3d.zmminglobal+top.zbeam: return
    # --- Check if there is any data. If there is none, then don't make
    # --- a dump.
    ntot = 0
    for js in range(self.getns()):
      ntot = ntot + self.getn(js=js)
    if ntot > 0 and me == 0:
      ff = None
#     try:
#       ff = PWpyt.PW(self.name+'_epdump.pyt')
#       dumpsmode = 1
#     except:
      ff = PW.PW(self.name+'_epdump.pdb')
      dumpsmode = 0
      if ff is None:
         print "ExtPart: %s unable to dump data to file."%self.name
         return
      ff.write(self.name+'@pickle',cPickle.dumps(self,dumpsmode))
      ff.close()
    self.nepmax = 1
    self.clear()
    # --- Disable is done last so that the object written out to the
    # --- file is still enabled. That flag is used in restoredata to
    # --- determine whether or not to restore the data. The logic is set
    # --- so that the object in an autodump file will restore the data
    # --- but one is a generic dump file won't (unless it was not auto
    # --- dumped, in which case the object in the generic dump is the only
    # --- copy). There will of course be exceptions, so restoredata takes
    # --- and option argument to force restoration of data, and the
    # --- extforcenorestore function turns any restores off.
    self.disable()

  def dodumptofile(self):
    if me != 0: return
    ff = None
#   try:
#     --- For now, pytables doesn't work since it has a limit of the number
#     --- of arrays that can be written out.
#     ff = PWpyt.PW(self.name+'_ep.pyt','a',verbose=0)
#   except:
    ff = PW.PW(self.name+'_ep.pdb','a',verbose=0)
    if ff is None:
       print "ExtPart: %s unable to dump data to file."%self.name
       return
    for js in range(top.ns):
      suffix = "_%d_%d"%(top.it,js)
      if self.getn(js=js) > 0:
        ff.write('n'+suffix,self.getn(js=js))
        ff.write('t'+suffix,self.gett(js=js))
        ff.write('x'+suffix,self.getx(js=js))
        ff.write('y'+suffix,self.gety(js=js))
        ff.write('ux'+suffix,self.getux(js=js))
        ff.write('uy'+suffix,self.getuy(js=js))
        ff.write('uz'+suffix,self.getuz(js=js))
    ff.close()

  def accumulate(self):
    # --- If top.nepwin is 0 then something is really wrong - this routine
    # --- should never be called if top.nepwin is zero.
    if top.nepwin == 0: return
    # --- Check if the number of species has changed. This is done to ensure
    # --- crashes don't happen.
    if top.ns > self.getns(): self.addspecies()
    # --- If this windows is outside of the grid, then just return.
    if (self.iz == -1 and 
        (self.zz+self.wz < w3d.zmminglobal+top.zbeam or
         self.zz-self.wz > w3d.zmmaxglobal+top.zbeam)):
      self.autodump()
      return
    # --- Make sure the arrays didn't overflow. Note that this is not an
    # --- issue for lab frame windows.
    if globalmax(maxnd(top.nep)) == top.nepmax:
      print "************* WARNING *************"
      print "**** Not enough space was allocated for the ExtPart arrays."
      print "**** The data will not be correct and will not be saved."
      print "**** A guess will be made as to how much to increase the size"
      print "**** of the arrays. Please run another timestep to accumulate new"
      print "**** data"
      if top.allocated("pnumz") and 0 <= self.iz <= top.nzmmnt:
        guess = 3*top.pnumz[self.iz]
      else:
        guess = 0
      # --- Only do this on if the relation is true. This avoids unnecessarily
      # --- increasing the size of the arrays on processors where no data
      # --- is gathered.
      if maxnd(top.nep) == top.nepmax:
        top.nepmax = max(2*top.nepmax,guess)
        err = gchange("ExtPart")
      return
    id = self.getid()
    # --- Loop over species, collecting only ones where some particles
    # --- were saved.
    for js in range(top.ns):
      # --- Gather the data.
      # --- In parallel, the data is gathered in PE0, return empty arrays
      # --- on other processors. In serial, the arrays are just returned as is.
      nn = top.nep[id,js]
      ntot = globalsum(nn)
      if ntot == 0: continue
      t = gatherarray(top.tep[:nn,id,js],othersempty=1)
      x = gatherarray(top.xep[:nn,id,js],othersempty=1)
      y = gatherarray(top.yep[:nn,id,js],othersempty=1)
      ux = gatherarray(top.uxep[:nn,id,js],othersempty=1)
      uy = gatherarray(top.uyep[:nn,id,js],othersempty=1)
      uz = gatherarray(top.uzep[:nn,id,js],othersempty=1)
      if self.laccumulate and not self.dumptofile:
        self.tep[js].append(t+0.)
        self.xep[js].append(x+0.)
        self.yep[js].append(y+0.)
        self.uxep[js].append(ux+0.)
        self.uyep[js].append(uy+0.)
        self.uzep[js].append(uz+0.)
      else:
        self.tep[js] = t
        self.xep[js] = x
        self.yep[js] = y
        self.uxep[js] = ux
        self.uyep[js] = uy
        self.uzep[js] = uz
    if self.dumptofile: self.dodumptofile()
    # --- Force nep to zero to ensure that particles are not saved twice.
    top.nep[id,:] = 0

  def setaccumulate(self,v=1):
    self.laccumulate = v
    if self.laccumulate: self.setuparrays(top.ns)

  ############################################################################
  def restoredata(self,lforce=0):
    """
Restores data dumped to a file. Note that this turns off the dumptofile
feature.
  - lforce=0: if true, force a restore, despite the value of enabled.
    """
    if not self.dumptofile: return
    if not lforce and (not self.enabled or _extforcenorestore): return
    self.dumptofile = 0
    self.laccumulate = 1
    try:
      ff = PRpyt.PR(self.name+'_ep.pyt','a',verbose=0)
    except:
      ff = PR.PR(self.name+'_ep.pdb','a',verbose=0)
    # --- Get total number of particles
    ntot = []
    jsmax = 0
    varlist = list(ff.inquire_names())
    varlist.sort()
    for var in varlist:
      if var[0] == 'n':
        name,ii,js = string.split(var,'_')
        jsmax = max(jsmax,eval(js))
        while jsmax >= len(ntot): ntot.append(0)
        ntot[jsmax] = ntot[jsmax] + ff.read(var)
    self.setuparrays(jsmax+1,bump=max(array(ntot))+1)
    for var in varlist:
      if var[0] == 'n':
        name,iis,jss = string.split(var,'_')
        nn = ff.read(var)
        ii = eval(iis)
        js = eval(jss)
        self.tep[js].append(ff.read('t_%d_%d'%(ii,js)))
        self.xep[js].append(ff.read('x_%d_%d'%(ii,js)))
        self.yep[js].append(ff.read('y_%d_%d'%(ii,js)))
        self.uxep[js].append(ff.read('ux_%d_%d'%(ii,js)))
        self.uyep[js].append(ff.read('uy_%d_%d'%(ii,js)))
        self.uzep[js].append(ff.read('uz_%d_%d'%(ii,js)))
    ff.close()

  ############################################################################
  def selectparticles(self,val,js=0,tc=None,wt=None,tp=None):
    if tc is None: return val[js][:]
    if wt is None: wt = self.dt
    if tp is None: tp = self.tep[js][:]
    ii = compress((tc-wt<tp)&(tp<tc+wt),arange(len(tp)))
    return take(val[js][:],ii)

  def getns(self): return len(self.tep)
  def gett(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.tep,js,tc,wt,tp)
  def getx(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.xep,js,tc,wt,tp)
  def gety(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.yep,js,tc,wt,tp)
  def getux(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uxep,js,tc,wt,tp)
  def getuy(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uyep,js,tc,wt,tp)
  def getuz(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uzep,js,tc,wt,tp)
  def getvx(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uxep,js,tc,wt,tp)
  def getvy(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uyep,js,tc,wt,tp)
  def getvz(self,js=0,tc=None,wt=None,tp=None):
    return self.selectparticles(self.uzep,js,tc,wt,tp)

  def getxp(self,js=0,tc=None,wt=None,tp=None):
    return self.getux(js,tc,wt,tp)/self.getuz(js,tc,wt,tp)
  def getyp(self,js=0,tc=None,wt=None,tp=None):
    return self.getuy(js,tc,wt,tp)/self.getuz(js,tc,wt,tp)
  def getr(self,js=0,tc=None,wt=None,tp=None):
    return sqrt(self.getx(js,tc,wt,tp)**2 + self.gety(js,tc,wt,tp)**2)
  def gettheta(self,js=0,tc=None,wt=None,tp=None):
    return arctan2(self.gety(js,tc,wt,tp),self.getx(js,tc,wt,tp))
  def getrp(self,js=0,tc=None,wt=None,tp=None):
    return (self.getxp(js,tc,wt,tp)*cos(self.gettheta(js,tc,wt,tp)) +
            self.getyp(js,tc,wt,tp)*sin(self.gettheta(js,tc,wt,tp)))
  def getn(self,js=0,tc=None,wt=None,tp=None):
    return len(self.gett(js,tc,wt,tp))

  def xxpslope(self,js=0,tc=None,wt=None,tp=None):
    if self.getn(js,tc,wt,tp) == 0:
      return 0.
    else:
      return ((ave(self.getx(js,tc,wt,tp)*self.getxp(js,tc,wt,tp)) -
               ave(self.getx(js,tc,wt,tp))*ave(self.getxp(js,tc,wt,tp)))/
              (ave(self.getx(js,tc,wt,tp)*self.getx(js,tc,wt,tp)) -
               ave(self.getx(js,tc,wt,tp))*ave(self.getx(js,tc,wt,tp))))
  def yypslope(self,js=0,tc=None,wt=None,tp=None):
    if self.getn(js,tc,wt,tp) == 0:
      return 0.
    else:
      return ((ave(self.gety(js,tc,wt,tp)*self.getyp(js,tc,wt,tp)) -
               ave(self.gety(js,tc,wt,tp))*ave(self.getyp(js,tc,wt,tp)))/
              (ave(self.gety(js,tc,wt,tp)*self.gety(js,tc,wt,tp)) -
               ave(self.gety(js,tc,wt,tp))*ave(self.gety(js,tc,wt,tp))))
  def rrpslope(self,js=0,tc=None,wt=None,tp=None):
    if self.getn(js,tc,wt,tp) == 0:
      return 0.
    else:
      return (ave(self.getr(js,tc,wt,tp)*self.getrp(js,tc,wt,tp))/
              ave(self.getr(js,tc,wt,tp)**2))

  ############################################################################
  ############################################################################
  # --- Define plotting routines for the extrapolated particles.

  def checkplotargs(self,kw):
    """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
    """
    badargs = ppgeneric(checkargs=1,kwdict=kw)
    kw['allowbadargs'] = 1
    if badargs: raise "bad arguments ",string.join(badargs.keys())

  def titleright(self,tc,wt):
    if tc is None:
      ttext = ''
    else:
      if wt is None: wt = self.dt
      ttext = "  time = %e ^+_-%e"%(tc,wt)
    if self.iz >= 0:
      ztext =  "iz = %d (z = %f m)"%(self.iz,w3d.zmminglobal+self.iz*w3d.dz)
    else:
      ztext =  "z = %f m"%self.zz
    return ztext + ttext

  ############################################################################
  def pxy(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots X-Y for extraploated particles"""
    self.checkplotargs(kw)
    x = self.getx(js,tc,wt,tp)
    y = self.gety(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
    settitles("Y vs X","X","Y",self.titleright(tc,wt))
    ppgeneric(y,x,kwdict=kw)

  ############################################################################
  def pxxp(self,js=0,tc=None,wt=None,tp=None,slope=0.,offset=0.,
           **kw):
    """Plots X-X' for extraploated particles"""
    self.checkplotargs(kw)
    x = self.getx(js,tc,wt,tp)
    xp = self.getxp(js,tc,wt,tp)
    if type(slope) == type(''):
      if len(x) > 0:
        slope = (ave(x*xp)-ave(x)*ave(xp))/(ave(x*x) - ave(x)**2)
        offset = ave(xp)-slope*ave(x)
      else:
        slope = 0.
        offset = 0.
    kw['slope'] = slope
    kw['offset'] = offset
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
    settitles("X' vs X","X","X'",self.titleright(tc,wt))
    ppgeneric(xp,x,kwdict=kw)

  ############################################################################
  def pyyp(self,js=0,tc=None,wt=None,tp=None,slope=0.,offset=0.,
           **kw):
    """Plots Y-Y' for extraploated particles"""
    self.checkplotargs(kw)
    y = self.gety(js,tc,wt,tp)
    yp = self.getyp(js,tc,wt,tp)
    if type(slope) == type(''):
      if len(y) > 0:
        slope = (ave(y*yp)-ave(y)*ave(yp))/(ave(y*y) - ave(y)**2)
        offset = ave(yp)-slope*ave(y)
      else:
        slope = 0.
        offset = 0.
    kw['slope'] = slope
    kw['offset'] = offset
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
    settitles("Y' vs Y","Y","Y'",self.titleright(tc,wt))
    ppgeneric(yp,y,kwdict=kw)

  ############################################################################
  def pxpyp(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots X'-Y' for extraploated particles"""
    self.checkplotargs(kw)
    xp = self.getxp(js,tc,wt,tp)
    yp = self.getyp(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
    settitles("Y' vs X'","X'","Y'",self.titleright(tc,wt))
    ppgeneric(yp,xp,kwdict=kw)

  ############################################################################
  def prrp(self,js=0,tc=None,wt=None,tp=None,scale=0.,slope=0.,offset=0.,
           **kw):
    """Plots R-R' for extraploated particles"""
    self.checkplotargs(kw)
    x = self.getx(js,tc,wt,tp)
    y = self.gety(js,tc,wt,tp)
    xp = self.getxp(js,tc,wt,tp)
    yp = self.getyp(js,tc,wt,tp)
    xscale = 1.
    yscale = 1.
    xpscale = 1.
    ypscale = 1.
    if scale:
      xscale = 2.*sqrt(ave(x*x) - ave(x)**2)
      yscale = 2.*sqrt(ave(y*y) - ave(y)**2)
      xpscale = 2.*sqrt(ave(xp*xp) - ave(xp)**2)
      ypscale = 2.*sqrt(ave(yp*yp) - ave(yp)**2)
    x = x/xscale
    y = y/yscale
    xp = xp/xpscale
    yp = yp/ypscale
    r = sqrt(x**2 + y**2)
    t = arctan2(y,x)
    rp = xp*cos(t) + yp*sin(t)
    if type(slope) == type(''):
      if len(r) > 0:
        slope = ave(r*rp)/ave(r*r)
      else:
        slope = 0.
    kw['slope'] = slope
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                        top.xpplmin/xpscale,top.xpplmax/ypscale)
    settitles("R' vs R","R","R'",self.titleright(tc,wt))
    ppgeneric(rp,r,kwdict=kw)

  ############################################################################
  def ptx(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-X for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    x = self.getx(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.xplmin,top.xplmax)
    settitles("X vs time","time","X",self.titleright(tc,wt))
    ppgeneric(x,t,kwdict=kw)

  ############################################################################
  def pty(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-Y for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    y = self.gety(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.yplmin,top.yplmax)
    settitles("Y vs time","time","Y",self.titleright(tc,wt))
    ppgeneric(y,t,kwdict=kw)

  ############################################################################
  def ptxp(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-X' for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    xp = self.getxp(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.xpplmin,top.xpplmax)
    settitles("X' vs time","time","X'",self.titleright(tc,wt))
    ppgeneric(xp,t,kwdict=kw)

  ############################################################################
  def ptyp(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-Y' for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    yp = self.getyp(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.ypplmin,top.ypplmax)
    settitles("Y' vs time","time","Y'",self.titleright(tc,wt))
    ppgeneric(yp,t,kwdict=kw)

  ############################################################################
  def ptux(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-ux for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    ux = self.getux(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("ux vs time","time","ux",self.titleright(tc,wt))
    ppgeneric(ux,t,kwdict=kw)

  ############################################################################
  def ptuy(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-uy for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    uy = self.getuy(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("uy vs time","time","uy",self.titleright(tc,wt))
    ppgeneric(uy,t,kwdict=kw)

  ############################################################################
  def ptuz(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-uz for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    uz = self.getuz(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("uz vs time","time","uz",self.titleright(tc,wt))
    ppgeneric(uz,t,kwdict=kw)

  ############################################################################
  def ptvx(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-Vx for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    vx = self.getvx(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vx vs time","time","Vx",self.titleright(tc,wt))
    ppgeneric(vx,t,kwdict=kw)

  ############################################################################
  def ptvy(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-Vy for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    vy = self.getvy(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vy vs time","time","Vy",self.titleright(tc,wt))
    ppgeneric(vy,t,kwdict=kw)

  ############################################################################
  def ptvz(self,js=0,tc=None,wt=None,tp=None,**kw):
    """Plots time-Vz for extraploated particles"""
    self.checkplotargs(kw)
    t = self.gett(js,tc,wt,tp)
    vz = self.getvz(js,tc,wt,tp)
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vz vs time","time","Vz",self.titleright(tc,wt))
    ppgeneric(vz,t,kwdict=kw)

  ############################################################################
  def ptrace(self,js=0,tc=None,wt=None,tp=None,slope=0.,
             pplimits=None,**kw):
    """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments for X-X' and Y-Y' plots.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
    """
    self.checkplotargs(kw)
    x = self.getx(js,tc,wt,tp)
    y = self.gety(js,tc,wt,tp)
    xp = self.getxp(js,tc,wt,tp)
    yp = self.getyp(js,tc,wt,tp)
    titler = self.titleright(tc,wt)
    defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                       (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                       (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                       (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
    if pplimits is None:
      pplimits = defaultpplimits
    else:
      kw['lframe'] = 1
      if type(pplimits[0]) != type(()):
        pplimits = 4*[pplimits]
      else:
        for i in range(4):
          if i == len(pplimits): pplimits.append(defaultpplimits[i])
          if not pplimits[i]: pplimits[i] = defaultpplimits[i]
 
    kw['view'] = 3
    kw['pplimits'] = pplimits[0]
    if type(slope)==type(''): kw['slope'] = 0.
    settitles("Y vs X","X","Y",titler)
    ppgeneric(y,x,kwdict=kw)
 
    kw['view'] = 4
    kw['pplimits'] = pplimits[1]
    if type(slope)==type(''):
      kw['slope'] = (ave(y*yp)-ave(y)*ave(yp))/dvnz(ave(y*y) - ave(y)**2)
    settitles("Y' vs Y","Y","Y'",titler)
    ppgeneric(yp,y,kwdict=kw)

    kw['view'] = 5
    kw['pplimits'] = pplimits[2]
    if type(slope)==type(''):
      kw['slope'] = (ave(x*xp)-ave(x)*ave(xp))/dvnz(ave(x*x) - ave(x)**2)
    settitles("X' vs X","X","X'",titler)
    ppgeneric(xp,x,kwdict=kw)
 
    kw['view'] = 6
    kw['pplimits'] = pplimits[3]
    if type(slope)==type(''): kw['slope'] = 0.
    settitles("X' vs Y'","Y'","X'",titler)
    ppgeneric(xp,yp,kwdict=kw)

##############################################################################
def dumpExtPart(object,filename):
  """Dump the saved extrapolated data to a file
 - filename: The name of the file to save the data in"""
  if me == 0:
    # --- Only PE0 writes the object to the file since it is the processor
    # --- where the data is gathered.
    ff = open(filename,'w')
    cPickle.dump(object,ff,1)
    ff.close()

def restoreExtPart(object,filename):
  """Restore extrapolated data from the given file"""
  if me == 0:
    # --- Only PE0 wrote the object to the file since it is the processor
    # --- where the data was gathered.
    ff = open(filename,'r')
    result = cPickle.load(ff)
    ff.close()
    result.enable()
    # --- Get the value of iz
    iz = result.iz
  else:
    # --- Create temp iz
    iz = 0
  # --- PE0 broadcasts its value of iz to all of the other processors
  # --- which create new instances of the ExtPart class.
  iz = broadcast(iz)
  if me > 0: result = ExtPart(iz)
  return result






