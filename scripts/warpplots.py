from warp import *
import RandomArray
import re
import os
import __main__
if me == 0:
  try:
    import pl3d
    import plwf
  except ImportError:
    pass
warpplots_version = "$Id: warpplots.py,v 1.106 2003/08/04 19:49:17 dave Exp $"

##########################################################################
# This setups the plot handling for warp.
##########################################################################

##########################################################################
warpplotsdocbasic = """
Basic graphics commands
winon(): creates X graphic windows
hcp(): send current plot to hard-copy file
fma(): do a frame advance
plg(): basic plotting routine
plp(): plots markers (dots) instead of lines
pla(): plots multi-dimensional array as a series of lines
plotc(): contour plots 2-D data
plotfc(): contour plots 2-D data with colored contour levels
limits(): sets plot limits in order left, right, bottom, top
zoombox(): when called, use the mouse left button (press and hold) to draw a
           box around the area to be zoomed to.
mouse commmands: left button zoom in, middle shifts, right zoom out

These return or set a slice out of the rho or phi array.
getrho(), getphi(), setrho(), setphi()

The following plot various particles projections.
ppzxy(), ppzx(), ppzy(), ppzr(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz(), ppxy()
ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy(), ppxvz(), ppyvz(), pprrp()
pprvz(), pptrace()

The following plot various particles projections using color.
ppzxco(), ppzyco(), ppzxyco(), ppzvzco()

Plots arbitrary particle projections using color
ppco()

Plots various quantities versus z in beam frame (see pzplotsdoc())

Run histplotsdoc() for a list of history plots.

Plots solution from envelope code.
penv()

Plots contours of charge density (rho) or electrostatic potential (phi) in
various planes.
pcrhozy(), pcrhozx(), pcrhoxy()
pcphizy(), pcphizx(), pcphixy()
pcselfezy(), pcselfezx(), pcselfexy()

Dynamically view any gist 3-D surface plot
viewsurface

Remove extra surface plots
hidesurfaces()
"""
##########################################################################
warpplotsdocmore = """
setup(): does the work needed to start writing plots to a file automatically
plotruninfo(): plots run info at bottom of plots (called by fma and hcp)
beforeplot[] list of functions called after each frame advance
afterplot[] list of functions called just before each frame advance
nf() = fma() Emulation of Basis command
sf() = redraw() Emulation of Basis command
warpplg(): for plotting data distributed across processors
warpplp(): for plotting particles distributed across processors
plothist(): convenience routine for plotting generic history data

Define variable names for various colors
fg, bg, white, black, red, green, blue, cyan, magenta, yellow 

Define variable names for plot markers
point, plus, star, circle

settitles(): set plot titles
ptitles(): draw plot titles on the current frame

plps[]: list of interpreter defined plots controlled by itplps
plfreq[]: list of interpreter defined plots controlled by itplfreq
plseldom[]: list of interpreter defined plots controlled by itplseldom
plalways[]: list of interpreter defined plots controlled by itplalways
pltfld3d(): makes fields plots which have been turned on
onedplts(): makes 1-D plots which have been turned on
psplots(): makes particle phase space plots which have been turned on
"""
##########################################################################

def warpplotsdoc():
  print warpplotsdocbasic
  print warpplotsdocmore

##########################################################################

##########################################################################
# The setup routine does the work needed to start writing plots to a file
# automatically.
top.lpsplots = true
always = top.always
seldom = top.seldom
never = top.never
cgmlogfile = None
numframes = 0
if me == 0: pldefault(marks=0) # --- Set plot defaults, no line marks
def setup(makepsfile=0,prefix=None,cgmlog=1,runcomments='',
          cgmfilesize=100000):
  """
Does the work needed to start writing plots to a file automatically
  - makepsfile=0: allows the specification of a ps file instead of cgm
  - prefix=None: optional prefix to use for plotfile name instead of runid
  - cgmlog=1: Set to 0 to inhibit cgmlog file creation
  - runcomments='': Additional comments to print on the first plot page
  - cgmfilesize=100000: Max cgmfilesize in units of MBytes.
  """
  # --- cgmlogfile is needed elsewhere
  global cgmlogfile
  # --- Only PE0 (or serial processor) should run this routine.
  if me > 0: return
  # --- Set cgmfilesize
  try:
    pldefault(cgmfilesize=cgmfilesize)
  except:
    pass
  # --- Get next available plot file name.
  if not prefix: prefix = arraytostr(top.runid)
  if makepsfile:
    pname = getnextfilename(prefix,'ps')
  else:
    pname = getnextfilename(prefix,'cgm')
  # --- Create window(0), but have it only dump to the file pname for now.
  # --- Note that only plots made to window(0) are dumped to the file.
  window(0,display='',hcp=pname,dump=1)
  print "Plot file name",pname
  # --- Set so all fma's dump plot to file.
  hcpon()
  if cgmlog:
    # --- Create plot log file and write headint to it.
    plogname = getnextfilename(prefix,'cgmlog')
    cgmlogfile = open(plogname,"w")
    cgmlogfile.write("CGMLOG file for "+pname+"\n\n")
  # --- Print the versions to the plot file.
  plt(versionstext()+'\n'+runcomments,0.15,0.88,justify="LT")
  fma()

# --- Convenience function to open a window with default value specilized to
# --- WARP. By default, this opens up a window on the current display. If
# --- setup has been called, this just creates a window which is attached to
# --- the already created device. Otherwise, open a window attached to a
# --- new device.
def winon(winnum=0,dpi=100):
  """
Opens up an X window
  - winnum=0 is the window number
  - dpi=100 is the dots per inch (either 100 or 75)
  """
  if winnum==0:
    window(winnum,dpi=dpi,display=os.environ['DISPLAY'])
  else:
    window(winnum,dpi=dpi)

##########################################################################
# Plot run info to the current plot and plot info to the log file.
framet=''
frameb=''
framel=''
framer=''
def plotruninfo():
  "Plot run info to the current plot and plot info to the log file"
  global numframes
  ss = (arraytostr(top.pline3)+'\n'+
        arraytostr(top.pline2)+'\n'+
        arraytostr(top.pline1))
  ss = re.sub(r'x10\|S2\|','e',ss)
  plt(ss,0.12,0.28)
  runmaker = arraytostr(top.runmaker)
  codeid = arraytostr(top.codeid)
  rundate = arraytostr(top.rundate)
  runtime = arraytostr(top.runtime)
  runid = arraytostr(top.runid)
  ss = '%-28s  %-8s  %-8s  %-9s  %-8s'%(runmaker,codeid,rundate,runtime,runid)
  plt(ss,0.12,0.24)
  if current_window()==0:
    # --- Only increment and print frame number and log if the active
    # --- device is window(0).
    numframes = numframes + 1
    plt(repr(numframes),0.68,0.9,justify='RA')
    if cgmlogfile:
      cgmlogfile.write('%d Step %d %s %s %s %s\n' %
                       (numframes,top.it,framet,frameb,framel,framer))

##########################################################################
# Frame advance and redraw routines. The fma routine from gist is replaced
# with one that prints informative text at the bottom of each frame just
# before the normal gist fma is called. Also created are alternate (Basis
# like) names for fma and redraw.
gistfma = fma
gisthcp = hcp
afterplot = []
beforeplot = []
def fma(legend=1):
  """
Frame advance - plots run info on the bottom of the frame, gets graphics window
ready for next plot and sends image to hard copy file if one is opened. Checks
for before and after plot commands.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  if legend: plotruninfo()
  for f in afterplot: f()
  gistfma()
  for f in beforeplot: f()
  oldlimits = limits()
def hcp(legend=1):
  """
Hardcopy - plots run info on the bottom of the frame and sends image to hard
copy file.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  for f in afterplot: f()
  if legend: plotruninfo()
  for f in beforeplot: f()
  gisthcp()

nf = fma
sf = redraw

##########################################################################
# Create the plotting routines. It is different in the serial and parallel
# versions.
if not lparallel:
  def warpplp(y,x,color="fg",linetype="none",marker="\1",msize=1.0):
    "Plots particles, same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=linetype,marker=marker,msize=msize)
  def warpplg(y,x,color="fg",linetype="solid",marks=0,marker=None,msize=1.0,
              width=1.0):
    "Same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=linetype,marks=marks,marker=marker,msize=msize,
          width=width)
else:
  warpplp = plotpart
  warpplg = plotarray

# --- Plot particles
circle = '\4'
star = '\3'
plus = '\2'
point = '\1'
def plp(y,x=None,color="fg",marker="\1",msize=1.0):
  """Plots particles, same as plg but with different defaults so it plots
markers instead of lines"""
  if len(y) == 0: return
  if x is not None:
    plg(y,x,type="none",marker=marker,color=color,msize=msize)
  else:
    plg(y,type="none",marker=marker,color=color,msize=msize)

# --- Plot history data. Convenience function that is only needed until
# --- the 'limited' capability is implemented.
def plothist(v,iw):
  """
Plots any history versus z
   - v is the history array
   - iw is the window number to plot
  """
  plg(v[iw,0:top.jhist+1],top.hzbeam[0:top.jhist+1])

# --- Simple interface to contour plotting. Only requires the 2-D array
# --- to be plotted.
def plotc(zz,xx=None,yy=None,ireg=None,color='fg',levs=None,contours=8,
          filled=0,width=1.,linetype='solid',cmin=None,cmax=None):
  """
Simple interface to contour plotting, same arguments as plc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours=8 Optional number of levels or list of levels
  - filled=0 When 1, draws filled contours
  - cmin, cmax: min and max of contours levels
  """
  s = shape(zz)
  if len(s) != 2:
    print 'First argument must be a 2-Dimensional array'
    return
  if not xx:
    xx = arange(s[0])[:,NewAxis]*ones(s[1],'d')
  elif len(shape(xx))==1:
    xx = xx[:,NewAxis]*ones(s[1],'d')
  if not yy:
    yy = arange(s[1])*ones(s[0],'d')[:,NewAxis]
  elif len(shape(yy))==1:
    yy = yy*ones(s[0],'d')[:,NewAxis]
  if not ireg:
    ireg = ones(s)
  if levs is not None: contours = levs
  if type(contours) == ListType: contours = array(contours)
  if type(contours) == TupleType: contours = array(contours)
  if type(contours) == type(1):
    if cmin is None: cmin = minnd(zz)
    if cmax is None: cmax = maxnd(zz)
    contours = 1.*iota(0,contours)*(cmax-cmin)/contours + cmin
  if filled:
    plfc(zz,xx,yy,ireg,contours=contours)
  else:
    plc(zz,xx,yy,color=color,levs=contours,width=width,type=linetype)
def plotfc(zz,xx=None,yy=None,ireg=None,contours=8):
  """
Simple interface to filled contour plotting, same arguments as plfc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours Optional number of levels or list of levels
  """
  plotc(zz,xx=xx,yy=yy,ireg=ireg,color=color,contours=contours,filled=1)

# --- Define variables names for the allowed colors
fg = 'fg'
bg = 'bg'
white = 'white'
black = 'black'
red = 'red'
green = 'green'
blue = 'blue'
cyan = 'cyan'
magenta = 'magenta'
yellow = 'yellow'

########################################################################
########################################################################
########################################################################
# The next part of this file contains Python routines which emulate compiled
# plot routines.
#
# Here are the plots available so far
#
# ppzx(), ppzy(), ppzr(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz()
# ppxy(), ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy()
# ppxvz(), ppyvz(), pprvz(), ppzxy(), pptrace()
# ppco(y,x,z;uz,xmin,xmax,ymin,ymax,zmin,zmax)
#
# The following only work properly serially
#
# ppzxco(), ppzyco(), ppzvzco(), ppzxyco()
#
##########################################################################
# List of available named colors.
color = ["red","green","blue","cyan","magenta","yellow"]

########################################################################
# Note: Subtracted off 0.0337 from X position of titlel (10/21/99)
ptitle_placement = [
  [[0.3950, 0.8863], [0.3950, 0.3927], [0.1200, 0.6800], [0.397, 0.37]],
  [[0.3950, 0.8863], [0.3950, 0.3927], [0.1200, 0.6800], [0.397, 0.37]],
  [[0.2634, 0.8863], [0.2634, 0.6559], [0.1300, 0.8006], [0.397, 0.37]],
  [[0.5266, 0.8863], [0.5266, 0.6559], [0.3932, 0.8006], [0.397, 0.37]],
  [[0.2634, 0.6231], [0.2634, 0.3927], [0.1300, 0.5374], [0.397, 0.37]],
  [[0.5266, 0.6231], [0.5266, 0.3927], [0.3932, 0.5374], [0.397, 0.37]],
  [[0.2634, 0.8863], [0.2634, 0.3927], [0.1300, 0.6800], [0.397, 0.37]],
  [[0.5266, 0.8863], [0.5266, 0.3927], [0.3932, 0.6800], [0.397, 0.37]],
  [[0.3950, 0.8863], [0.3950, 0.6559], [0.1300, 0.8006], [0.397, 0.37]],
  [[0.3950, 0.6231], [0.3950, 0.3927], [0.1300, 0.5374], [0.397, 0.37]]]
default_titlet=""
default_titleb=""
default_titlel=""
default_titler=""
def settitles(titlet="",titleb="",titlel="",titler=""):
  "Sets titles which are plotted by ptitles"
  global default_titlet,default_titleb,default_titlel,default_titler
  default_titlet = titlet
  default_titleb = titleb
  default_titlel = titlel
  default_titler = titler
def ptitles(titlet="",titleb="",titlel="",titler="",v=1):
  "Plots titles, either uses input or titles set by settitles"
  global framet,frameb,framel,framer
  if titlet=="" and default_titlet: titlet = default_titlet
  if titleb=="" and default_titleb: titleb = default_titleb
  if titlel=="" and default_titlel: titlel = default_titlel
  if titler=="" and default_titler: titler = default_titler
  framet=titlet
  frameb=titleb
  framel=titlel
  framer=titler
  if titlet:
    plt(titlet,ptitle_placement[v-1][0][0],ptitle_placement[v-1][0][1],
        justify="CC",orient=0)
  if titleb:
    plt(titleb,ptitle_placement[v-1][1][0],ptitle_placement[v-1][1][1],
        justify="CC",orient=0)
  if titlel:
    plt(titlel,ptitle_placement[v-1][2][0],ptitle_placement[v-1][2][1],
        justify="CC",orient=1)
  if titler:
    plt(titler,ptitle_placement[v-1][3][0],ptitle_placement[v-1][3][1],
        justify="CC",orient=0)
  settitles()
def ptitlebottom(text=""):
  plt(text,0.3950,0.37,justify="CC")

##########################################################################
##########################   UTILITY ROUTINES  ###########################
##########################################################################
##########################################################################
def checkarguments(input,arglist):
  "Compare inputs against and argument list and return list of bad arguments"
  inputcopy = input.copy()
  for i in inputcopy.keys():
    if i in arglist.keys(): del inputcopy[i]
  return inputcopy

##########################################################################
def pptitleright(iw=0,kwdict={},**kw):
  "Returns right plot title. Takes same arguments as selectparticles"
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,"slope":0,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  if badargs and not allowbadargs:
    raise "bad argument",string.join(badargs.keys())

  # --- Return appropriate right title
  if zl is not None or zu is not None:
    if z is None: prefix = ""
    else: prefix = "z "
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    result = prefix+"range (%9.4e, %9.4e)"%(zl,zu)
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    result = "ix = %d, x range (%9.4e, %9.4e)"%(ix,xl,xu)
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    result = "iy = %d, y range (%9.4e, %9.4e)"%(iy,yl,yu)
  elif iz is not None:
    zl = w3d.zmminglobal + iz*w3d.dz - wz*w3d.dz + top.zbeam
    zu = w3d.zmminglobal + iz*w3d.dz + wz*w3d.dz + top.zbeam
    result = "iz = %d, z range (%9.4e, %9.4e)"%(iz,zl,zu)
  elif iw < 0:
    if psubset==[]: setup_subsets()
    result = "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles"
  else:
    if win is None:
      win = top.zwindows[:,iw] + top.zbeam
      prefix = "z "
    else:
      prefix = ""
    if len(shape(win)) == 2: win = win[:,iw]
    result = prefix+"window%d = %9.4e, %9.4e"%(iw,win[0],win[1])
  if slope != 0:
    result = result + ", slope=%7.4f"%slope
  return result

#-------------------------------------------------------------------------
def ppmoments(text):
  "Plots text in upper right hand corner of the plot"
  plt(text,0.61,.855,justify="RT",height=12,font="courierB")

#############################################################################
#############################################################################
vpythonscenelist = []
#############################################################################
def ppgeneric_doc(x,y):
  doc = selectparticles.__doc__ + """
  - zz: optional third particle data quantity - when supplied, it is deposited
       on a grid and that is used for contour levels.
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - nx, ny: grid size, defaults to 20x20
  - slope=0.: slope to subtract from %(y)s coordinate (%(y)s-slope*%(x)s)
  - xoffset=0.: average %(x)s of particles
  - yoffset=0.: average %(y)s of particles
  - xscale=1.: scaling factor applied to x data
  - yscale=1.: scaling factor applied to y data
  - titles=1: when true, plot the titles
  - lframe=0: when true, the plot limits are set to the plmin and plmax input
              arguments, which default to the plmin and plmax variables from
              the group InDiag
  - pplimits=None: a tuple of (xplmin, xplmax, yplmin, yplmax), limits of plot
                   range (used when lframe=1)
  - xmin, xmax, ymin, ymax: extrema of density grid, defaults to particle
                            extrema (x for %(x)s and y for %(y)s)
  - cmin=min(grid), cmax=max(grid): min and max of data for coloration
  - xbound=dirichlet: sets boundary condition on gridded data for x
  - ybound=dirichlet: sets boundary condition on gridded data for y
  - particles=0: when true, plot particles
  - uselog=0: when true, logarithmic levels of the number density are used
  - color='fg': color of particles, when=='density', color by number density
  - ncolor=None: when plotting particle color by number density, number of
                 colors to use, defaults to top.ncolor
  - denmin, denmax: thresholds for removing particles, only particles located
                    where density is between denmin and denmax are plotted
  - chopped=None: only particles where r < chopped*maxdensity/density
                  are plotted, where r is a random number between 0 and 1
                  and density is the density at the particle location
  - marker=dot: particle marker to plot
  - msize=1.: scaled size of marker
  - hash=0: flag to turn on or off the hash plot
  - line_scale=.9: scaling factor on line length
  - hcolor='fg': color of hash marks
  - width=1.0: width of hash marks
  - contours=None: number of countours to plot
  - filled=0: when true, plot filled contours
  - ccolor='fg': contour color (when not filled)
  - cellarray=0: when true, plot grid as cell array
  - centering='node': centering of cells with cellarray, other option are 'cell'                      and 'old' (for old incorrect scaling)
  - ctop=199: max color index for cellarray plot
  - ldensityscale=0: when true, scale the density by its max.
  - view=1: view window to use (experts only)
  - lcolorbar=1: when plotting colorized data, include a colorbar
  - colbarunitless=0: when true, color-bar scale is unitless
  - colbarlinear=1: when true, the color-bar is laid out linearly in density,
                    otherwise each contour level gets an equal sized area.
                    Only in effect when a list of colorbars is specified.
  - surface=0: when true, a 3-d surface plot is made of the gridded data
               Note: to remove window, use the hidesurfaces() command
                     rather than closing the window.
  - returngrid=0: when true, and when particle data is passed in and a plot
                  which requires a grid is requested (such as a contour
                  plot), no plotting is done and the grid and extrema
                  are returned in a tuple
  """
  return doc%vars()
#-------------------------------------------------------------------------
def ppgeneric(y=None,x=None,kwdict={},**kw):
  """
Generic particle plotting routine. Allows plotting of particle points, density
contours, and/or density hash marks.
Note that either the x and y coordinates or the grid must be passed in.
  - y, x: optional particle data (instead of using inputted grid)
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {'zz':None,'weights':None,'grid':None,
                'nx':20,'ny':20,'slope':0.,
                'xoffset':0.,'yoffset':0.,'offset':0.,
                'xscale':1.,'yscale':1.,'titles':1,'lframe':0,
                'xmin':None,'xmax':None,'ymin':None,'ymax':None,
                'pplimits':('e','e','e','e'),
                'particles':0,'uselog':0,'color':'fg','ncolor':top.ncolor,
                'usepalette':1,'marker':'\1','msize':1.0,
                'denmin':None,'denmax':None,'chopped':None,
                'hash':0,'line_scale':.9,'hcolor':'fg','width':1.0,
                'contours':None,'filled':0,'ccolor':'fg',
                'cellarray':0,'centering':'node','ctop':199,
                'cmin':None,'cmax':None,
                'xbound':dirichlet,'ybound':dirichlet,
                'ldensityscale':0,'view':1,
                'lcolorbar':1,'colbarunitless':0,'colbarlinear':1,'surface':0,
                'xmesh':None,'ymesh':None,
                'returngrid':0,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  assert (not badargs or allowbadargs), \
         "bad argument: %s"%string.join(badargs.keys())

  # --- Do some error checking on the consistency of the input
  assert (type(grid) == ArrayType or \
          (type(x) == ArrayType and type(y) == ArrayType)), \
         "either the grid and/or both x and y must be specified"
  assert (not particles or (type(x) == ArrayType and type(y) == ArrayType)), \
         "both x and y must be specified if particles are to be plotted"
  assert ((type(x) != ArrayType and type(y) != ArrayType) or len(x) == len(y)),\
         "both x and y must be of the same length"
  assert (zz is None) or (type(zz) == ArrayType and len(zz) == len(x)),\
         "zz must be the same length as x"
  assert (type(slope) != StringType),"slope must be a number"
  assert (zz is None) or (grid is None),\
         "only one of zz and grid can be specified"
  assert (centering == 'node' or centering == 'cell' or centering == 'old'),\
         "centering must take one of the values 'node', 'cell', or 'old'"

  # --- If there are no particles and no grid to plot, just return
  if type(x) == ArrayType and type(y) == ArrayType: np = globalsum(len(x))
  else: np = 0
  if np == 0 and grid is None: return

  # --- If filled is turned on, but contours is not set, set it to the
  # --- default value of 8.
  if filled and contours is None: contours = 8

  # --- If particle data was passed in and no specific plots were requested,
  # --- just plot the particles.
  if y is not None and \
     (not hash and not contours and not surface and not cellarray):
    particles = 1

  # --- If a grid is passed in and no specific plots were requested,
  # --- make a cellarray plot.
  if grid is not None and \
     (not hash and not contours and not surface and not cellarray):
    cellarray = 1

  # --- Make sure that nothing is not plotted over a surface plot
  if surface:
    particles = 0
    contours = 0
    cellarray = 0
    hash = 0
    lframe = 0
    titles = 0

  # -- Set the plotting view window
  plsys(view)

  # --- Make sure that the grid size nx and ny are consistent with grid
  # --- is one is input
  if type(grid) == ArrayType:
    nx = shape(grid)[0] - 1
    ny = shape(grid)[1] - 1

  # --- Calculate extrema of the particles
  if type(x) == ArrayType and type(y) == ArrayType:
    # --- Get slope subtracted value of y
    yms = y - (x-xoffset)*slope - yoffset - offset
    # --- Get mins and maxs of particles that were not supplied by the user.
    if lparallel:
      if xmin is None: xmintemp = globalmin(x)
      if xmax is None: xmaxtemp = globalmax(x)
      if ymin is None: ymintemp = globalmin(yms)
      if ymax is None: ymaxtemp = globalmax(yms)
    else:
      xmintemp = 0.
      xmaxtemp = 0.
      ymintemp = 0.
      ymaxtemp = 0.
      if xmin is None and len(x) > 0: xmintemp = min(x)
      if xmax is None and len(x) > 0: xmaxtemp = max(x)
      if ymin is None and len(yms) > 0: ymintemp = min(yms)
      if ymax is None and len(yms) > 0: ymaxtemp = max(yms)
    # --- When neither the min or max are supplied by the user, extend
    # --- extrema by one grid cell so that all particles are within the
    # --- limits of the grid. This is the most common case.
    if xmin is None and xmax is None:
      xmintemp = xmintemp - (xmaxtemp-xmintemp)/(nx-2)
      xmaxtemp = xmaxtemp + (xmaxtemp-xmintemp)/(nx-2)
    if ymin is None and ymax is None:
      ymintemp = ymintemp - (ymaxtemp-ymintemp)/(ny-2)
      ymaxtemp = ymaxtemp + (ymaxtemp-ymintemp)/(ny-2)
    # --- Now set main versions of min and max
    if xmin is None: xmin = xmintemp
    if xmax is None: xmax = xmaxtemp
    if ymin is None: ymin = ymintemp
    if ymax is None: ymax = ymaxtemp

    # --- Scale the data
    x = x*xscale
    yms = yms*yscale
  else:
    # --- If no particles are inputted and the extrema are not set, then
    # --- can only make a guess.
    if xmin is None: xmin = 0
    if xmax is None: xmax = nx
    if ymin is None: ymin = 0
    if ymax is None: ymax = ny

  # --- Scale the extrema
  xmin = xmin*xscale
  xmax = xmax*xscale
  ymin = ymin*yscale
  ymax = ymax*yscale

  # --- Get grid cell sizes
  if nx != 0: dx = (xmax-xmin)/nx
  else:       dx = 1.
  if ny != 0: dy = (ymax-ymin)/ny
  else:       dy = 1.

  # --- If the grid is needed for the plot and it was not passed in, generate
  # --- it from the inputted particle data (if there was any)
  if type(grid) != ArrayType and \
     (hash or contours or color=='density' or chopped or denmin or denmax or surface or cellarray):
    if zz is None:
      densitygrid = 1

      # --- Create space for data
      grid = fzeros((1+nx,1+ny),'d')

      # --- Deposit the density onto the grid.
      if(weights is None):
        setgrid2d(len(x),x,yms,nx,ny,grid,xmin,xmax,ymin,ymax)
      else:
        setgrid2dw(len(x),x,yms,weights,nx,ny,grid,xmin,xmax,ymin,ymax)
      # --- If parallel, do a reduction on the grid
      grid = parallelsum(grid)

    else:
      densitygrid = 0

      # --- Create space for data
      grid = fzeros((1+nx,1+ny),'d')
      gridcount = fzeros((1+nx,1+ny),'d')

      # --- Deposit the data onto the grid. itask is 1 so that the parallel
      # --- version can be done properly.
      if(weights is None):
        deposgrid2d(1,len(x),x,yms,zz,nx,ny,grid,gridcount,xmin,xmax,ymin,ymax)
      else:
        deposgrid2dw(1,len(x),x,yms,zz,weights,nx,ny,grid,gridcount,xmin,xmax,ymin,ymax)

      # --- If parallel, do a reduction on the grid
      grid = parallelsum(grid)
      gridcount = parallelsum(gridcount)

      # --- Divide out the particle counts by hand.
      grid = grid/where(greater(gridcount,0.),gridcount,1.)

    # --- Enforce boundary conditions on the grid
    if xbound == neumann:
      grid[0,:] = 2.*grid[0,:]
      grid[-1,:] = 2.*grid[-1,:]
    elif xbound == periodic:
      grid[0,:] = grid[0,:] + grid[-1,:]
      grid[-1,:] = grid[0,:]
    if ybound == neumann:
      grid[:,0] = 2.*grid[:,0]
      grid[:,-1] = 2.*grid[:,-1]
    elif ybound == periodic:
      grid[:,0] = grid[:,0] + grid[:,-1]
      grid[:,-1] = grid[:,0]

    # --- If requested, return the grid and extrema, doing no plotting
    if returngrid: return (grid,xmin,xmax,ymin,ymax)

  elif (hash or contours or color=='density' or chopped or denmin or denmax
        or surface):
    densitygrid = 0
 
  # --- Scale the grid by its maximum if requested.
  if ldensityscale and grid is not None:
    gridmax = maxnd(abs(grid))
    if gridmax != 0.:
      grid[:,:] = grid/gridmax

  # --- If using logarithmic number density levels, take the log of the grid
  # --- data. The original grid is left unchanged since that is still needed
  # --- by some operations below.
  if uselog:
    if densitygrid:
      # --- Take the log, raising all values below 0.1 to 0.1. The
      # --- threshold is used so that none of the elements are zero.
      # --- That value 0.1 is used since values any smaller do not have
      # --- much meaning since a value of 1.0 means that there is already
      # --- only one particle in that cell.
      grid1 = log10(where(less(grid,0.1),0.1,grid))
    else:
      # --- Before taking the log of the user supplied grid data, make sure
      # --- that there are no negative values. Zero is ok since they will
      # --- be replaced with a minimum value.
      dmax = maxnd(grid)
      dmin = minnd(where(equal(grid,0.),dmax,grid))
      if dmin <= 0.:
        raise "Can't take log since the grid has negative values"
      grid1 = log(where(less(grid,dmin/10.),dmin/10.,grid))
  else:
    grid1 = grid

  # --- Get grid min and max and generate contour levels if needed.
  if (hash or contours or color=='density' or chopped or surface or cellarray):
    if cmin is None: cmin = minnd(grid1)
    if cmax is None: cmax = maxnd(grid1)

  # --- Get grid mesh if it is needed
  if contours or hash or surface or cellarray:
    if xmesh is not None or ymesh is not None: usermesh = 1
    else:                                      usermesh = 0
    if xmesh is None: xmesh = xmin + dx*arange(nx+1)[:,NewAxis]*ones(ny+1,'d')
    if ymesh is None: ymesh = ymin + dy*arange(ny+1)*ones(nx+1,'d')[:,NewAxis]

  # --- Make filled contour plot of grid first since it covers everything
  # --- plotted before it.
  if contours and filled and nx > 1 and ny > 1:
    if cmax != cmin:
      plotc(transpose(grid1),transpose(ymesh),transpose(xmesh),
            color=ccolor,contours=contours,filled=filled,cmin=cmin,cmax=cmax)

  # --- Make cell-array plot. This also is done early since it covers anything
  # --- done before it. The min and max are adjusted so that the patch for
  # --- each grid cell has the correct centering.
  # --- If the user supplies a mesh, then use plf, a filled mesh plot, since
  # --- the meshes may not be Cartesian.
  if cellarray and nx > 1 and ny > 1:
    if not usermesh:
      if centering == 'node':
        xminc = xmin
        xmaxc = xmax + dx
        yminc = ymin
        ymaxc = ymax + dy
      elif centering == 'cell':
        xminc = xmin - dx/2.
        xmaxc = xmax + dx/2.
        yminc = ymin - dy/2.
        ymaxc = ymax + dy/2.
      elif centering == 'old':
        xminc = xmin
        xmaxc = xmax
        yminc = ymin
        ymaxc = ymax
      pli(transpose(grid1),xminc,yminc,xmaxc,ymaxc,top=ctop,cmin=cmin,cmax=cmax)
    else:
      plf(grid1,ymesh,xmesh)

  # --- Plot particles
  if particles:
    if color == 'density':
      z1 = zeros(len(x),'d')
      getgrid2d(len(x),x,yms,z1,nx,ny,grid1,xmin,xmax,ymin,ymax)
    if chopped or denmin or denmax:
      dd = zeros(len(x),'d')
      getgrid2d(len(x),x,yms,dd,nx,ny,grid,xmin,xmax,ymin,ymax)
      maxdensity = maxnd(grid)
      dd = dd/maxdensity
      ipick = ones(shape(x))
      if chopped:
        ipick[:] = ipick*less(RandomArray.random(shape(x)),chopped/dd)
      if denmin:
        ipick[:] = ipick*less(denmin,dd)
      if denmax:
        ipick[:] = ipick*less(dd,denmax)
      x = compress(ipick,x)
      yms = compress(ipick,yms)
      if color == 'density':
        z1 = compress(ipick,z1)
    if color == 'density':
      # --- Plot particles with color based on the density from the grid.
      ppco(yms,x,z1,uz=1.,marker=marker,msize=msize,lframe=0,
           xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=cmin,zmax=cmax,
           ncolor=ncolor,usepalette=usepalette)
    else:
      # --- Plot particles as a solid color.
      warpplp(yms,x,color=color,marker=marker,msize=msize)

  # --- Now plot unfilled contours, which are easier to see on top of the
  # --- particles
  if contours and not filled and nx > 1 and ny > 1:
    if cmax != cmin:
      plotc(transpose(grid1),transpose(ymesh),transpose(xmesh),
            color=ccolor,contours=contours,filled=filled,cmin=cmin,cmax=cmax)

  # --- Plot hash last since it easiest seen on top of everything else.
  if hash:
    # --- Set line length
    if nx != 0 and cmax != cmin:
      sss = line_scale*(xmax-xmin)/nx/(cmax - cmin)
    else:
      sss = 1.
    # --- Make plot of tick marks
    for ix in range(nx+1):
      for iy in range(ny+1):
        plg(ymesh[ix,iy]+zeros(2),xmesh[ix,iy]+array([0.,sss*grid1[ix,iy]]),
            color=hcolor,width=width)

  # --- Add colorbar if needed
  if lcolorbar and \
     ((contours and filled==1) or (color == 'density' and len(x) > 0) or \
      (cellarray)):
    if (contours and filled==1):
      try:
        nc = len(contours) + 1
        levs = contours
      except TypeError:
        nc = contours
        levs = None
    elif (color == 'density' and len(x) > 0):
      nc = ncolor
      levs = None
    elif (cellarray):
      nc = ctop
      levs = None
    if colbarunitless:
      dmin = 0.
      dmax = 1.0
    elif cellarray:
      dmin = cmin
      dmax = cmax
    else:
      dmin = cmin
      dmax = cmax
    colorbar(dmin,dmax,uselog=uselog,ncolor=nc,view=view,levs=levs,
             colbarlinear=colbarlinear,ctop=ctop)

  # --- Make surface plot
  if surface and me == 0 and nx > 1 and ny > 1:
    try:
      import VPythonobjects
      if type(color) != ListType: scolor = None
      else:                       scolor = color
      xrange = 1.5*max(abs(xmin),abs(xmax))
      yrange = 1.5*max(abs(ymin),abs(ymax))
      zrange = 1.5*maxnd(abs(grid1))
      vo = VPythonobjects.VisualMesh(zvalues=grid1,title=default_titlet,
                                     color=scolor,vrange=(xrange,yrange,zrange))
      vpythonscenelist.append(vo.scene)
    except ImportError:
      pl3d.orient3()
      pl3d.light3()
      plwf.plwf(grid1,xmesh,ymesh,fill=grid1,edges=0)
      [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
      #limits(xmin3,xmax3,ymin3,ymax3)

  # --- Finish off the plot, adding titles and setting the frame limits.
  if titles: ptitles(v=view)
  settitles() 
  if (lframe):
    ppp = list(pplimits)
    if ppp[0] != 'e': ppp[0] = ppp[0]*xscale
    if ppp[1] != 'e': ppp[1] = ppp[1]*xscale
    if ppp[2] != 'e': ppp[2] = ppp[2]*yscale
    if ppp[3] != 'e': ppp[3] = ppp[3]*yscale
    limits(ppp[0],ppp[1],ppp[2],ppp[3])
if sys.version[:5] != "1.5.1":
  ppgeneric.__doc__ = ppgeneric.__doc__ + ppgeneric_doc('x','y')


#############################################################################
#############################################################################
def ppvector(gridy=None,gridx=None,kwdict={},**kw):
  """
Generic vector plotting routine.
Note that both the x and y grids must be passed in.
  - gridy, gridx: x and y vector comnponents
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {'titles':1,'lframe':0,
                'xmin':None,'xmax':None,'ymin':None,'ymax':None,
                'pplimits':('e','e','e','e'),'scale':1.,
                'color':'fg',
                'xbound':dirichlet,'ybound':dirichlet,
                'xmesh':None,'ymesh':None,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  assert (not badargs or allowbadargs), \
         "bad argument: %s"%string.join(badargs.keys())

  # --- Do some error checking on the consistency of the input
  assert gridx is not None and gridy is not None,"both gridx and gridy must be specified"

  nx = shape(gridx)[0] - 1
  ny = shape(gridx)[1] - 1

  assert (shape(gridy)[0] - 1) == nx and (shape(gridy)[1] - 1) == ny,"gridx and gridy must be the same shape"

  if xmin is None: xmin = 0
  if xmax is None: xmax = nx
  if ymin is None: ymin = 0
  if ymax is None: ymax = ny

  # --- Get meshes
  dx = (xmax - xmin)/nx
  dy = (ymax - ymin)/ny
  xx,yy = getmesh2d(xmin,dx,nx,ymin,dy,ny)

  # --- Compute scale
  scale = scale*min(dx,dy)/dvnz(max(maxnd(abs(gridx)),maxnd(abs(gridy))))
  print scale

  # --- Make plot
  plv(gridy,gridx,yy,xx,scale=scale)

#############################################################################
#############################################################################
# ed williams' colorbar stuff / modified for Warp by J.-L. Vay on 01/22/2001
def nicelevels(z,n=8) :
  """nicelevels(z,n=8) finds approximately n "nice values"
between min(z) and max(z) for axis labels. n defaults to eight.
  """
  zmax = max(ravel(z))
  zmin = min(ravel(z))
  if zmin == zmax: return array([zmin,zmax])
  finest = abs(zmax - zmin)/float (n)
  # blows up on zmin=zmax
  unit = 10.**floor (log10 (finest))
  finest = finest/unit
  if   finest > 5.0: finest = 10.
  elif finest > 2.:  finest = 5.
  elif finest > 1.:  finest = 2.
  unit = unit*finest
  cmin = unit*ceil(zmin/unit)
  if (abs(cmin - zmin) < 0.01*unit) :
     cmin = cmin
  cmax = unit*floor(zmax/unit)
  if (abs(cmax - zmax) < 0.01*unit) :
     cmax = cmax
  n = int(((cmax - cmin)/unit + 0.5) + 1)
  levs = cmin + arange(n)*unit
  llist = nonzero(less(abs(levs),0.1*unit))
  if len(llist) > 0:
     array_set(levs,llist,0.0)
  return levs

#-----------------------------------------------------------------------
colorbar_placement = [[0.62,0.64,0.43,0.86],[0.62,0.64,0.43,0.86],
                      [0.354,0.364,0.692,0.859],[0.617,0.627,0.692,0.859],
                      [0.354,0.364,0.428,0.596],[0.617,0.627,0.428,0.596],
                      [0.354,0.364,0.43,0.86],[0.617,0.627,0.43,0.86],
                      [0.617,0.627,0.692,0.859],[0.617,0.627,0.428,0.596]]
colorbar_fontsize = [14.,14.,8.,8.,8.,8.,8.,8.,8.,8.]

def colorbar(zmin,zmax,uselog=0,ncolor=100,view=1,levs=None,colbarlinear=1,
             ctop=199):
  """
Plots a color bar to the right of the plot square labelled by the z
values from zmin to zmax.
  - zmin, zmax: lower and upper range for color bar
  - uselog=0: when true, labels are printed in the form 10^x
  - ncolor=100: default number of colors to include
  - view=1: specifies the view that is associated with the color bar
  - levs: an optional list of color levels
  - ctop=199: number of colors from palette to use
  """
  plsys(0)
  xmin,xmax,ymin,ymax = colorbar_placement[view-1]
  fontsize = colorbar_fontsize[view-1]
  # --- draw the bar
  if colbarlinear and levs is not None:
    # --- Use the contour plotting routine plfc for this case. The data
    # --- plotted is uniformly spaced between zmin and zmax. The contour
    # --- levels are those specified. The result is that the colorbar
    # --- shows the contours levels by their values relative to zmin and zmax.
    plotval = span(zmin,zmax,255)[:,NewAxis]*ones(2)
    xx = array([xmin,xmax])*ones(255)[:,NewAxis]
    yy = span(ymin,ymax,255)[:,NewAxis]*ones(2)
    ireg = ones((255,2))
    plfc(plotval,yy,xx,ireg,contours=array(levs))
  else:
    # --- Use cell array plotting for this case. All of the colors get a block
    # --- of the same size. If levs is not specified, the uniform spacing 
    # --- matches the uniform spacing of the contours. If levs is specified,
    # --- each equal sized block represents one contour level, independent of
    # --- the range of the level relative to other levels.
    if type(zmin) == type(zmax) == type(1) and \
       zmin >= 0 and zmax <=199:
       plotval = arange(zmin,zmax+1,typecode='b')[:,NewAxis]*ones(2)
    else:
       plotval = (arange(ncolor)/(ncolor-1.))[:,NewAxis]*ones(2)
    pli(plotval,xmin,ymin,xmax,ymax,top=ctop)
  # --- Draw a black box around it
  pldj([xmin,xmin,xmin,xmax],[ymin,ymax,ymin,ymin],
       [xmax,xmax,xmin,xmax],[ymin,ymax,ymax,ymax])
  # --- Generate nice levels for the labels and tick marks.
  if levs is None:
    # --- Use the nicelevels routine to get evenly spaced labels.
    nicelevs = nicelevels(array([zmin,zmax]))
  else:
    # --- If there are less than 15 specified contour levels, put a label
    # --- at each of the labels. If there are more, pick out roughly 15
    # --- evenly spaced values. Also, if the levels do not extend beyond
    # --- zmin and zmax, add labels at those points too.
    nicelevs = levs
    if zmin < levs[0]:  nicelevs = array([zmin] + list(nicelevs))
    if zmax > levs[-1]: nicelevs = array(list(nicelevs) + [zmax])
  llev = len(nicelevs)
  # --- Create the labels
  labels = []
  # --- Calculate the location of the labels.
  if not colbarlinear and levs is not None:
    # --- The ys are evenly spaced
    ys = ymin + arange(len(nicelevs))/(len(levs)+1.)*(ymax - ymin)
    # --- If the lowest level is less than zmin, then bump up the y's
    # --- by one block size.
    if levs[0] < zmin: ys = ys + 1./(len(levs)+1.)*(ymax - ymin)
  elif llev==2 and (nicelevs[0] == nicelevs[1]):
    ys = array([ymin,ymax])
  else:
    ys = ymin + (ymax - ymin)*(nicelevs - zmin)/(zmax - zmin)
  # --- Plot the labels, skipping ones that are too close together.
  if uselog: ss = " 10^%.5g"
  else:      ss = " %.5g"
  ylast = 0.
  for i in xrange(llev):
    if ys[i] - ylast > (ymax-ymin)/30:
      plt(ss%nicelevs[i],xmax+0.005,ys[i]-0.005,height=fontsize)
      ylast = ys[i]
  # --- Plot the tick marks
  pldj(llev*[xmax],ys,llev*[xmax+0.005],ys)
  # --- Return to plot system 1.
  plsys(view)

#############################################################################
#############################################################################
def changepalette(returnpalette=0,help=0,view=1):
  """
Dynamically change the color palette.
  - returnpalette=0: when true, returns tuple of (red, green, blue)
  - help=0: when true, prints this message
Mouse actions:
  Button 1: shifts a point, compressing and stretching the rest of the colors
  Button 2: reset palette to original
  Button 3: shifts a point, sliding the colors up and down
  Shift Button 1: reverse the palette
  Control Button 1: add black point
  Control Button 3: add white point
  """
  # --- Print out help if wanted
  if help: print changepalette.__doc__
  # --- min's and max's are the same as in the colorbar routine
  xmin,xmax,ymin,ymax = colorbar_placement[view-1]
  # --- Create storate arrays
  # --- rr, gg, bb hold the original palette
  rr = zeros(200,'b')
  gg = zeros(200,'b')
  bb = zeros(200,'b')
  palette(rr,gg,bb,query=1)
  # --- newrr, newgg, newbb hold the new palette
  newrr = zeros(200,'b')
  newgg = zeros(200,'b')
  newbb = zeros(200,'b')
  # --- position relative to the original palette
  cc = arange(0,200)*1.
  newcc = arange(0,200)*1.
  # --- List of black and white points
  blacklist = []
  whitelist = []
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    # --- Get mouse positions. Skip if outside the colorbar
    (x1, y1, x2, y2) = tuple(mm[:4])
    if x1 < xmin or x1 > xmax or x2 < xmin or x2 > xmax: continue
    if y1 < ymin or y1 > ymax or y2 < ymin or y2 > ymax: continue

    if mm[9] == 1 and mm[10] == 0:
      # --- Button 1, no keys
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      up = (ymax - y1)/(ymax - y2)
      down = (y1 - ymin)/(y2 - ymin)
      for i in xrange(1,i2):
        iold = int(i*down)
        wold =     i*down - iold
        newcc[i] = cc[iold]*(1.-wold) + cc[iold+1]*wold
      for i in xrange(i2,199):
        iold = 199 - int((199-i)*up)
        wold = iold - (199 -    ((199-i)*up))
        newcc[i] = cc[iold]*(1.-wold) + cc[iold-1]*wold

    if mm[9] == 2:
      # --- Button 2, no keys
      # --- Restore original palette
      newcc = arange(0,200)*1.
      blacklist = []
      whitelist = []

    if mm[9] == 3:
      # --- Button 3, no keys
      # --- slide whole palette
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      for i in xrange(0,200):
        iold = i - (i2 - i1)
        if iold < 0: newcc[i] = cc[0]
        elif iold > 199: newcc[i] = cc[-1]
        else: newcc[i] = cc[iold]

    if mm[9] == 1 and mm[10] == 1:
      # --- Button 1, shift
      # --- Reverse the palette
      newcc[:] = cc[::-1]

    if mm[9] == 1 and mm[10] == 4:
      # --- button 1, control
      # --- Add black point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      blacklist.append(i1)

    if mm[9] == 3 and mm[10] == 4:
      # --- button 3, control
      # --- Add white point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      whitelist.append(i1)

    # --- Calculate the new palette based on the position relative to the
    # --- original palette.
    for i in xrange(0,200):
      ii = int(newcc[i])
      ww =     newcc[i]  - ii
      iip1 = min(ii+1,199)
      newrr[i] = (nint(rr[ii]*(1.-ww) + rr[iip1]*ww))
      newgg[i] = (nint(gg[ii]*(1.-ww) + gg[iip1]*ww))
      newbb[i] = (nint(bb[ii]*(1.-ww) + bb[iip1]*ww))
    for ii in blacklist: (newrr[ii], newgg[ii], newbb[ii]) = 0,0,0
    for ii in whitelist: (newrr[ii], newgg[ii], newbb[ii]) = 255,255,255
    cc[:] = newcc
    palette(newrr,newgg,newbb)
  if returnpalette: return (newrr,newgg,newbb)

#############################################################################
#############################################################################
def viewsurface(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
Returns the final values of the parameters that can be passed to pl3d.rot3
to reproduce the same orientation.
  """
  pl3d.gnomon(gnomon)
  [xmin3min,xmax3max,ymin3min,ymax3max,sys] = limits()
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    (xa, ya, za) = (0.,0.,0.)
    if mm[9] == 1:
      ya = - (mm[3] - mm[1])*scale
      za = - (mm[2] - mm[0])*scale
    if mm[9] == 3:
      xa = (mm[2] - mm[0])*scale
    pl3d.rot3(xa,ya,za)
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)
  print xa,ya,za

def _viewsurfacetest(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
  """
  pl3d.gnomon(gnomon)
  pl3d.orient3(phi=0.,theta=0.)
  [xmin3min,xmax3max,ymin3min,ymax3max] = pl3d.draw3(1)
  phi = 0.
  theta = 0.
  (xa, ya, za) = (0.,0.,0.)
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    dphi   = (mm[3] - mm[1])*scale
    dtheta = (mm[2] - mm[0])*scale
    print theta,phi
    newxa = xa + dtheta*sin(phi)*cos(theta) + dphi*cos(phi)*cos(theta)
    newya = ya + dtheta*sin(phi)*sin(theta) + dphi*cos(phi)*sin(theta)
    newza = za + dtheta*cos(phi)*cos(theta) + dphi*sin(phi)*sin(theta)
    phi = xa*cos(za) + ya*sin(za)
    theta = za
    pl3d.rot3(newxa-xa,newya-ya,newza-za)
    xa = newxa
    ya = newya
    za = newza
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)

#############################################################################
def zoombox():
  """
When called, use the mouse left button (press and hold) to draw a
box around the area to be zoomed to.
  """
  m1 = mouse(1,1,'')
  xmin = min(m1[0],m1[2])
  xmax = max(m1[0],m1[2])
  ymin = min(m1[1],m1[3])
  ymax = max(m1[1],m1[3])
  limits(xmin,xmax,ymin,ymax)

#############################################################################
def hidesurfaces():
  for scene in vpythonscenelist:
    scene.visible = 0

#############################################################################
#############################################################################
def ppmultispecies(pp,args,kw):
  """checks if js defined and assign it to a list if plotting multispecies.
  Also assign colors accordingly
  """
  if kw.has_key('js'):
    js = kw['js']
    if js != -1 and type(js) != ListType:
      return false
    else:
      if js == -1: js = range(top.ns)
      color = kw.get('color',range(0,240,240/len(js)))
      for i in xrange(len(js)):
        kw['js'] = js[i]
        kw['color'] = color[i]
        apply(pp,args,kw)
      return true
  else:
    return false

########################################################################
########################################################################
########################################################################
########################################################################
def checkparticleplotarguments(kw):
  """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
  """
  badargs = selectparticles(checkargs=1,kwdict=kw)
  badargs = pptitleright(checkargs=1,kwdict=badargs)
  badargs = ppgeneric(checkargs=1,kwdict=badargs)
  kw['allowbadargs'] = 1
  if badargs: raise "bad arguments",string.join(badargs.keys())
########################################################################
def ppzxy(iw=0,**kw):
  "Plots Z-X and Z-Y in single page"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxy,(iw,),kw): return
  kw['view'] = 9
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,win=top.ywindows,z=top.yp,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.xp,ii),take(top.zp,ii),kwdict=kw)

  kw['view'] = 10
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,win=top.xwindows,z=top.xp,kwdict=kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.yp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzxy.__doc__ = ppzxy.__doc__ + ppgeneric_doc('z','x')

##########################################################################
def ppzx(iw=0,**kw):
  "Plots Z-X"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzx,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.xp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzx.__doc__ = ppzx.__doc__ + ppgeneric_doc('z','x')

##########################################################################
def ppzy(iw=0,**kw):
  "Plots Z-Y"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzy,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.yp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzy.__doc__ = ppzy.__doc__ + ppgeneric_doc('z','y')

##########################################################################
def ppzr(iw=0,**kw):
  "Plots Z-R"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzr,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("R vs Z","Z","R",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(sqrt(top.xp**2+top.yp**2),ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzr.__doc__ = ppzr.__doc__ + ppgeneric_doc('z','r')

##########################################################################
def ppzxp(iw=0,**kw):
  "Plots Z-X'"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxp,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin,top.xpplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("X' vs Z","Z","X'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uxp,ii)/take(top.uzp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzxp.__doc__ = ppzxp.__doc__ + ppgeneric_doc('z',"x'")

##########################################################################
def ppzvx(iw=0,**kw):
  "Plots Z-Vx"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvx,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vx vs Z","Z","Vx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uxp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvx.__doc__ = ppzvx.__doc__ + ppgeneric_doc('z',"vx")

##########################################################################
def ppzyp(iw=0,**kw):
  "Plots Z-Y'"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzyp,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin,top.ypplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Y' vs Z","Z","Y'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uyp,ii)/take(top.uzp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzyp.__doc__ = ppzyp.__doc__ + ppgeneric_doc('z',"y'")

##########################################################################
def ppzvy(iw=0,**kw):
  "Plots Z-Vy"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvy,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vy vs Z","Z","Vy",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uyp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvy.__doc__ = ppzvy.__doc__ + ppgeneric_doc('z',"vy")

##########################################################################
def ppzvz(iw=0,**kw):
  "Plots Z-Vz"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange()
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vz vs Z","Z","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvz.__doc__ = ppzvz.__doc__ + ppgeneric_doc('z',"vz")

##########################################################################
def ppxy(iw=0,**kw):
  "Plots X-Y"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxy,(iw,),kw): return
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Y vs X","X","Y",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.yp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxy.__doc__ = ppxy.__doc__ + ppgeneric_doc('x','y')

##########################################################################
def ppxxp(iw=0,**kw):
  "Plots X-X'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxxp,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz',None))
    kw['slope'] = slope
    kw['yoffset'] = xpoffset
    kw['xoffset'] = xoffset
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("X' vs X","X","X'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uxp,ii)/take(top.uzp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxxp.__doc__ = ppxxp.__doc__ + ppgeneric_doc("x","x'")

##########################################################################
def ppyyp(iw=0,**kw):
  "Plots Y-Y'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyyp,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz',None))
    kw['slope'] = slope
    kw['yoffset'] = ypoffset
    kw['xoffset'] = yoffset
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Y' vs Y","Y","Y'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uyp,ii)/take(top.uzp,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyyp.__doc__ = ppyyp.__doc__ + ppgeneric_doc("y","y'")

##########################################################################
def ppxpyp(iw=0,**kw):
  "Plots X'-Y'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxpyp,(iw,),kw): return
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    (xslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz',None))
    (yslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz',None))
    kw['slope'] = 0.
  else:
    (xslope,xoffset,xpoffset) = (slope,0.,0.)
    (yslope,yoffset,ypoffset) = (slope,0.,0.)
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
  settitles("Y' vs X'","X'","Y'",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  xp = (take(top.uxp,ii)/take(top.uzp,ii) - xslope*(take(top.xp,ii)-xoffset) -
       xpoffset)
  yp = (take(top.uyp,ii)/take(top.uzp,ii) - yslope*(take(top.yp,ii)-yoffset) -
       ypoffset)
  return ppgeneric(yp,xp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxpyp.__doc__ = ppxpyp.__doc__ + ppgeneric_doc("x'","y'")

##########################################################################
def ppxvx(iw=0,**kw):
  "Plots X-Vx. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvx,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz',None))
    kw['slope'] = slope*vz
    kw['yoffset'] = xpoffset*vz
    kw['xoffset'] = xoffset
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vx vs X","X","Vx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uxp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxvx.__doc__ = ppxvx.__doc__ + ppgeneric_doc("x","Vx")

##########################################################################
def ppyvy(iw=0,**kw):
  "Plots Y-Vy. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvy,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz',None))
    kw['slope'] = slope*vz
    kw['yoffset'] = ypoffset*vz
    kw['xoffset'] = yoffset
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vy vs Y","Y","Vy",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uyp,ii)*take(top.gaminv,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyvy.__doc__ = ppyvy.__doc__ + ppgeneric_doc("y","Vy")

##########################################################################
def ppxvz(iw=0,**kw):
  "Plots X-Vz."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange()
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vz vs X","X","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxvz.__doc__ = ppxvz.__doc__ + ppgeneric_doc("x","Vz")

##########################################################################
def ppyvz(iw=0,**kw):
  "Plots Y-Vz."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange()
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vz vs Y","Y","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyvz.__doc__ = ppyvz.__doc__ + ppgeneric_doc("y","Vz")

##########################################################################
def pprrp(iw=0,scale=0,**kw):
  """Plots R-R', If slope='auto', it is calculated from the moments.
  - scale=0: when true, scale particle by 2*rms"""
  checkparticleplotarguments(kw)
  if ppmultispecies(pprrp,(iw,scale),kw): return
  xscale = 1.
  yscale = 1.
  xpscale = 1.
  ypscale = 1.
  if scale:
    iiw = max(0,iw)
    xscale = 2.*top.xrms[iiw]
    yscale = 2.*top.yrms[iiw]
    xpscale = 2.*top.vxrms[iiw]/top.vzbar[iiw]
    ypscale = 2.*top.vyrms[iiw]/top.vzbar[iiw]
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = take(top.xp,ii)/xscale
  yy = take(top.yp,ii)/yscale
  xp = take(top.uxp,ii)/take(top.uzp,ii)/xpscale
  yp = take(top.uyp,ii)/take(top.uzp,ii)/ypscale
  rr = sqrt(xx**2 + yy**2)
  tt = arctan2(yy,xx)
  rp = xp*cos(tt) + yp*sin(tt)
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    aversq = globalave(rr**2)
    averrp = globalave(rr*rp)
    if aversq > 0.:
      slope = averrp/aversq
    else:
      slope = 0.
    kw['slope'] = slope
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                      top.xpplmin/xpscale,top.xpplmax/ypscale)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("R' vs R","R","R'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(rp,rr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pprrp.__doc__ = pprrp.__doc__ + ppgeneric_doc("r","r'")

##########################################################################
def pprvz(iw=0,**kw):
  "Plots R-Vz"
  checkparticleplotarguments(kw)
  if ppmultispecies(pprvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange()
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax,top.yplmax),vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  settitles("Vz vs R","R","Vz",pptitleright(iw=iw,kwdict=kw))
  rr = sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2)
  vz = take(top.uzp,ii)*take(top.gaminv,ii)
  return ppgeneric(vz,rr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pprvz.__doc__ = pprvz.__doc__ + ppgeneric_doc("r","vz")

##########################################################################
def pptrace(iw=0,normalize=0,**kw):
  """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pptrace,(iw,normalize),kw): return
  ii = selectparticles(iw=iw,kwdict=kw)
  x = take(top.xp,ii)
  y = take(top.yp,ii)
  xp = take(top.uxp,ii)/take(top.uzp,ii)
  yp = take(top.uyp,ii)/take(top.uzp,ii)
  if(top.wpid!=0):kw['weights']=take(top.pid[:,top.wpid-1],ii)
  slope = kw.get('slope',0.)
  if type(slope)==type(''):
    del kw['slope']
    iz = kw.get('iz',None)
    (xxpslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=iz)
    (yypslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=iz)
    xp = xp - xxpslope*(x - xoffset) - xpoffset
    yp = yp - yypslope*(y - yoffset) - ypoffset
  if kw.get('titles',1):
    titler=pptitleright(iw=iw,kwdict=kw)
    ptitles(titler=titler)
  defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                     (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                     (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                     (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
  pplimits = kw.get('pplimits',None)
  if pplimits is None:
    pplimits = defaultpplimits
  else:
    kw['lframe'] = 1
    if type(pplimits[0]) != type(()):
      pplimits = 4*[pplimits]
    else:
      for i in xrange(4):
        if i == len(pplimits): pplimits.append(defaultpplimits[i])
        if not pplimits[i]: pplimits[i] = defaultpplimits[i]

  # --- First make plot with returngrid set to true. If no grids are used,
  # --- then ppgeneric returns None (and makes the plot), otherwise it
  # --- returns the grid. Also, if normalize is false, then set so plots
  # --- are made (and grid is not returned).
  if normalize: rg = 1
  else:         rg = 0
  kw['view'] = 3
  kw['pplimits'] = pplimits[0]
  settitles("Y vs X","X","Y")
  gxy = ppgeneric(y,x,returngrid=rg,kwdict=kw)
 
  kw['view'] = 4
  kw['pplimits'] = pplimits[1]
  settitles("Y' vs Y","Y","Y'")
  gyyp = ppgeneric(yp,y,returngrid=rg,kwdict=kw)
 
  kw['view'] = 5
  kw['pplimits'] = pplimits[2]
  settitles("X' vs X","X","X'")
  gxxp = ppgeneric(xp,x,returngrid=rg,kwdict=kw)
 
  kw['view'] = 6
  kw['pplimits'] = pplimits[3]
  settitles("X' vs Y'","Y'","X'")
  gxpyp = ppgeneric(xp,yp,returngrid=rg,kwdict=kw)

  # --- If the return value is None, then return since plots have already been
  # --- made.
  if gxy is None: return

  # --- If the return value is not None, then call ppgeneric again to 
  # --- actually make the plots with the appropriate cmin and cmax
  cmin = kw.get('cmin',None)
  cmax = kw.get('cmax',None)
  if kw.get('cmin',None) is None:
    kw['cmin']=min(minnd(gxy[0]),minnd(gyyp[0]),minnd(gxxp[0]),minnd(gxpyp[0]))
  if kw.get('cmax',None) is None:
    kw['cmax']=max(maxnd(gxy[0]),maxnd(gyyp[0]),maxnd(gxxp[0]),maxnd(gxpyp[0]))

  kw['view'] = 3
  kw['pplimits'] = pplimits[0]
  settitles("Y vs X","X","Y")
  ppgeneric(y,x,kwdict=kw)
 
  kw['view'] = 4
  kw['pplimits'] = pplimits[1]
  settitles("Y' vs Y","Y","Y'")
  ppgeneric(yp,y,kwdict=kw)
 
  kw['view'] = 5
  kw['pplimits'] = pplimits[2]
  settitles("X' vs X","X","X'")
  ppgeneric(xp,x,kwdict=kw)
 
  kw['view'] = 6
  kw['pplimits'] = pplimits[3]
  settitles("X' vs Y'","Y'","X'")
  ppgeneric(xp,yp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pptrace.__doc__ = pptrace.__doc__ + ppgeneric_doc("x","x'")

##########################################################################
##########################################################################
##########################################################################
def ppzxco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1,
           ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-X with color based in particle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","X"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xplmin,top.xplmax)

##########################################################################
def ppzyco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1,
           ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-Y with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
           color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","Y"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.yplmin,top.yplmax)

##########################################################################
def ppzxyco(js=0,marker="\1",msize=1.0,lframe=0,titles=1,
            ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-X and Z-Y in single frame with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  plsys(9)
  ppzxco(js,sys=9,marker=marker,msize=msize,lframe=lframe,titles=titles,
         ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol)
  plsys(10)
  ppzyco(js,sys=10,marker=marker,msize=msize,lframe=lframe,titles=titles,
         ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol)

##########################################################################
def ppzvzco(js=0,marker="\1",msize=1.0,lframe=0,titles=1,
            ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-Vz with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  (vzmin,vzmax) = getvzrange()
  istart = top.ins[js]-1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles("Vz vs Z","Z","Vz"," ")
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)

##########################################################################
def ppco(y,x,z,uz=1.,xmin=None,xmax=None,ymin=None,ymax=None,
         zmin=None,zmax=None,ncolor=None,usepalette=1,
         marker="\1",msize=1.0,lframe=0):
  """Plots y versus x with color based in z
     - y is y coordinate
     - x is x coordinate
     - z is used to calculate the color
     - xmin, xmax, ymin, ymax, zmin, zmax optional bounds
     - ncolor is number of colors to use, defaults to top.ncolor
     - usepalette=1 when true, uses palette, otherwise uses colors in array
                    color
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles
  """

  # --- If there are not particles to plot, just return
  np = globalsum(len(x))
  if np == 0: return

  # --- Make sure the lengths of the input are the same
  assert (len(y) == len(x) == len(z)),"x, y, and z must all be the same length"

  # --- This routine can be expensive in parallel when there are many
  # --- colors since synchronization is needed for each color.
  # --- So, if there arn't too many particles, transfer everything to PE0
  # --- and let it do the work.
  if np < 1000000:
    alllocal = 1
    y = gatherarray(y)
    x = gatherarray(x)
    z = gatherarray(z)
    if type(uz) == ArrayType: uz = gatherarray(uz)
  else:
    alllocal = 0

  # --- Make sure arrays are 1-D
  rx = ravel(x)
  ry = ravel(y)
  rz = ravel(z)

  # --- Find extrema
  if xmin is None: xmin = globalmin(rx)
  if xmax is None: xmax = globalmax(rx)
  if ymin is None: ymin = globalmin(ry)
  if ymax is None: ymax = globalmax(ry)
  if zmin is None: zmin = globalmin(rz)
  if zmax is None: zmax = globalmax(rz)

  if ncolor is None: ncolor = top.ncolor
  dd = (zmax - zmin)/ncolor
  for ic in xrange(ncolor):
    ii = compress(logical_and(logical_and(not_equal(uz,0.),
           less(zmin+ic*dd,rz)),less(rz,zmin+(ic+1)*dd)), iota(0,len(rx)))
    if usepalette:
      c = nint(199*ic/(ncolor-1.))
    else:
      c = color[ic%len(color)]
    if alllocal:
      plp(take(y,ii),take(x,ii),color=c,marker=marker,msize=msize)
    else:
      warpplp(take(y,ii),take(x,ii),
              color=c,linetype="none",marker=marker,msize=msize)
  if (lframe): limits(xmin,xmax,ymin,ymax)


##########################################################################
# To be implemented
#ppzx4
#ppzy4
#ppzxp4
#ppzyp4
#ppzvz4
#ppxy4
#ppxxp4
#ppyyp4
#ppxpyp4
#ppxxpco
#ppyypco


##########################################################################
# This routine allows plotting of multi-dimensioned arrays.
def pla(y,x=None,color="fg",linetype="solid",marks=0,marker=None,msize=1.0,
        width=1.,decomposed=0):
  "Same as plg but can plot multidimensional array"
  if x:
    if shape(x)!=shape(y) and (len(shape(x))==1 and shape(x)[0]!=shape(y)[0]):
      raise TypeError,"pla: x must either be the same shape as y, or it must be 1-D with the same len as the 1st dimension of y."
  else:
    x = arange(0,y.shape[0],1,'d')
  if len(shape(x)) > 2:
    # This next 'if' is needed since non-contiguous arrays cannot be reshaped.
    if x.iscontiguous():
      xx = x
    else:
      xx = x + 0.
    # Change xx into a 2-D array, with all of the upper dims lumped into
    # the second dimension.
    xx.shape = (xx.shape[0],product(array(xx.shape[1:])))
  elif len(shape(x)) == 2:
    # The input x is usable as is.
    xx = x
  else:
    # Extend xx into a 2-D array, with a second dimension of length 1.
    xx=x[:,NewAxis]
  if len(shape(y)) > 2:
    # This next 'if' is needed since non-contiguous arrays cannot be reshaped.
    if y.iscontiguous():
      yy = y
    else:
      yy = y + 0.
    # Change yy into a 2-D array, with all of the upper dims lumped into
    # the second dimension.
    yy.shape = (yy.shape[0],product(array(yy.shape[1:])))
  elif len(shape(y)) == 2:
    # The input y is usable as is.
    yy = y
  else:
    # Extend yy into a 2-D array, with a second dimension of length 1.
    yy=y[:,NewAxis]
  n = shape(xx)[1]
  for i in xrange(yy.shape[1]):
    if decomposed:
      warpplg(yy[:,i],xx[:,i%n],color=color,linetype=linetype,
              marks=marks,marker=marker,msize=msize,width=width)
    else:
      if len(yy[:,i]) > 0:
        plg(yy[:,i],xx[:,i%n],color=color,type=linetype,
            marks=marks,marker=marker,msize=msize,width=width)

##########################################################################
##########################################################################
# history plotting routines have been replaced by those in histplots.py
##########################################################################
##########################################################################

##########################################################################
def penv(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a and b envelope
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles"""
  if not me==0: return
  plg(env.aenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(env.benv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Envelope","Z")
##########################################################################
def penvp(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a' and b' of envelope
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles"""
  if not me==0: return
  plg(env.apenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(env.bpenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Envelope slope","Z")
##########################################################################
def penvaedge(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a envelope +/- x centroid
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles"""
  if not me==0: return
  plg(+env.aenv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(-env.aenv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("X Envelope edges","Z")
##########################################################################
def penvbedge(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots b envelope +/- x centroid
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles"""
  if not me==0: return
  plg(+env.benv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(-env.benv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Y Envelope edges","Z")
##########################################################################
##########################################################################
# --- These functions returns or sets slices of phi and rho.
##########################################################################
def getrho(ix=None,iy=None,iz=None,bcast=0):
  """Returns slices of rho, the charge density array. The shape of the object
returned depends on the number of arguments specified, which can be from none
to all three.
  - ix = None
  - iy = None
  - iz = None
  - bcast=0: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      return w3d.rho
    if ix is not None and iy is None     and iz is None    :
      return w3d.rho[ix,:,:]
    if ix is None     and iy is not None and iz is None    :
      return w3d.rho[:,iy,:]
    if ix is None     and iy is None     and iz is not None:
      return w3d.rho[:,:,iz]
    if ix is not None and iy is not None and iz is None    :
      return w3d.rho[ix,iy,:]
    if ix is not None and iy is None     and iz is not None:
      return w3d.rho[ix,:,iz]
    if ix is None     and iy is not None and iz is not None:
      return w3d.rho[:,iy,iz]
    if ix is not None and iy is not None and iz is not None:
      return w3d.rho[ix,iy,iz]
  else:
    iz1 = top.izfsslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    ppp = w3d.rho[:,:,iz1:iz2]
    if ix is not None and iy is None:
      ppp = ppp[ix,:,:]
    elif ix is None and iy is not None:
      ppp = ppp[:,iy,:]
    elif ix is not None and iy is not None:
      ppp = ppp[ix,iy,:]
    if iz is None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertizfstope(iz)
      if pe is None: return None
      if me == pe: ppp = ppp[...,iz-top.izfsslave[me]]
      else:        ppp = zeros(shape(ppp[...,0]),'d')
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = broadcast(ppp)
    return ppp
# --------------------------------------------------------------------------
def setrho(val,ix=None,iy=None,iz=None):
  """Sets slices of rho, the charge density array. The shape of the input
object depends on the number of arguments specified, which can be from none
to all three.
  - val input array (must be supplied)
  - ix = None
  - iy = None
  - iz = None
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      w3d.rho[:,:,:] = val
    if ix is not None and iy is None     and iz is None    :
      w3d.rho[ix,:,:] = val
    if ix is None     and iy is not None and iz is None    :
      w3d.rho[:,iy,:] = val
    if ix is None     and iy is None     and iz is not None:
      w3d.rho[:,:,iz] = val
    if ix is not None and iy is not None and iz is None    :
      w3d.rho[ix,iy,:] = val
    if ix is not None and iy is None     and iz is not None:
      w3d.rho[ix,:,iz] = val
    if ix is None     and iy is not None and iz is not None:
      w3d.rho[:,iy,iz] = val
    if ix is not None and iy is not None and iz is not None:
      w3d.rho[ix,iy,iz] = val
  else:
    print "Warning, setrho this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.rho[:,:,:-top.grid_overlap]
   #else:
   #  ppp = w3d.rho[:,:,:]
   #if ix is not None and iy is None    :
   #  ppp = ppp[ix,:,:]
   #elif ix is None and iy is not None:
   #  ppp = ppp[:,iy,:]
   #elif ix is not None and iy is not None:
   #  ppp = ppp[ix,iy,:]
   #if iz is None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe is None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = broadcast(ppp)
   #return ppp
# --------------------------------------------------------------------------
def getphi(ix=None,iy=None,iz=None,bcast=0):
  """Returns slices of phi, the electrostatic potential array. The shape of
the object returned depends on the number of arguments specified, which can
be from none to all three.
  - ix = None
  - iy = None
  - iz = None Value is relative to the fortran indexing, so iz ranges
              from -1 to nz+1
  - bcast=0: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      return w3d.phi[:,:,1:-1]
    if ix is not None and iy is None     and iz is None    :
      return w3d.phi[ix,:,1:-1]
    if ix is None     and iy is not None and iz is None    :
      return w3d.phi[:,iy,1:-1]
    if ix is None     and iy is None     and iz is not None:
      return w3d.phi[:,:,iz+1]
    if ix is not None and iy is not None and iz is None    :
      return w3d.phi[ix,iy,1:-1]
    if ix is not None and iy is None     and iz is not None:
      return w3d.phi[ix,:,iz+1]
    if ix is None     and iy is not None and iz is not None:
      return w3d.phi[:,iy,iz+1]
    if ix is not None and iy is not None and iz is not None:
      return w3d.phi[ix,iy,iz+1]
  else:
    iz1 = top.izfsslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    ppp = w3d.phi[:,:,iz1+1:iz2+1]
    if ix is not None and iy is None:
      ppp = ppp[ix,:,:]
    elif ix is None and iy is not None:
      ppp = ppp[:,iy,:]
    elif ix is not None and iy is not None:
      ppp = ppp[ix,iy,:]
    if iz is None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertizfstope(iz)
      if pe is None: return None
      if me == pe: ppp = ppp[...,iz-top.izfsslave[me]]
      else:        ppp = zeros(shape(ppp[...,0]),'d')
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = broadcast(ppp)
    return ppp
# --------------------------------------------------------------------------
def setphi(val,ix=None,iy=None,iz=None):
  """Sets slices of phi, the electrostatic potential array. The shape of
the input object depends on the number of arguments specified, which can
be from none to all three.
  - val input array (must be supplied)
  - ix = None
  - iy = None
  - iz = None Value is relative to the fortran indexing, so iz ranges
              from -1 to nz+1
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      w3d.phi[:,:,1:-1] = val
    if ix is not None and iy is None     and iz is None    :
      w3d.phi[ix,:,1:-1] = val
    if ix is None     and iy is not None and iz is None    :
      w3d.phi[:,iy,1:-1] = val
    if ix is None     and iy is None     and iz is not None:
      w3d.phi[:,:,iz+1] = val
    if ix is not None and iy is not None and iz is None    :
      w3d.phi[ix,iy,1:-1] = val
    if ix is not None and iy is None     and iz is not None:
      w3d.phi[ix,:,iz+1] = val
    if ix is None     and iy is not None and iz is not None:
      w3d.phi[:,iy,iz+1] = val
    if ix is not None and iy is not None and iz is not None:
      w3d.phi[ix,iy,iz+1] = val
  else:
    print "Warning, setphi this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.phi[:,:,1:w3d.nz-top.grid_overlap+2]
   #else:
   #  ppp = w3d.phi[:,:,1:-1]
   #if ix is not None and iy is None:
   #  ppp = ppp[ix,:,:]
   #elif ix is None and iy is not None:
   #  ppp = ppp[:,iy,:]
   #elif ix is not None and iy is not None:
   #  ppp = ppp[ix,iy,:]
   #if iz is None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe is None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = broadcast(ppp)
   #return ppp
# --------------------------------------------------------------------------
def getselfe(comp=None,ix=None,iy=None,iz=None,bcast=0,fullplane=0):
  """Returns slices of selfe, the electrostatic field array. The shape of
the object returned depends on the number of arguments specified, which can
be from none to all three.
  - comp: field component to get, either 'x', 'y', or 'z'
  - ix = None
  - iy = None
  - iz = None
  - bcast=0: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  """
  assert comp in ['x','y','z'],"comp must be one of 'x', 'y', or 'z'"
  if top.efetch != 3 or w3d.nx_selfe == 0:
    # --- If not already using selfe, then allocate it and set it.
    # --- Note that this could be an unexpected expense for a user.
    w3d.nx_selfe = w3d.nxp
    w3d.ny_selfe = w3d.nyp
    w3d.nz_selfe = w3d.nzp
    if w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.XZgeom:
      w3d.ny_selfe = 0
    if w3d.solvergeom==w3d.Zgeom: w3d.nx_selfe = 0
    gchange("Efields3d")
    getselfe3d(w3d.phip,w3d.nxp,w3d.nyp,w3d.nzp,w3d.selfe,
               w3d.nx_selfe,w3d.ny_selfe,w3d.nz_selfe,w3d.dx,w3d.dy,w3d.dz,
               top.pboundxy)
  ic = ['x','y','z'].index(comp)
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      eee = w3d.selfe[ic,:,:,:]
    elif ix is not None and iy is None     and iz is None    :
      eee = w3d.selfe[ic,ix,:,:]
    elif ix is None     and iy is not None and iz is None    :
      eee = w3d.selfe[ic,:,iy,:]
    elif ix is None     and iy is None     and iz is not None:
      eee = w3d.selfe[ic,:,:,iz]
    elif ix is not None and iy is not None and iz is None    :
      eee = w3d.selfe[ic,ix,iy,:]
    elif ix is not None and iy is None     and iz is not None:
      eee = w3d.selfe[ic,ix,:,iz]
    elif ix is None     and iy is not None and iz is not None:
      eee = w3d.selfe[ic,:,iy,iz]
    elif ix is not None and iy is not None and iz is not None:
      eee = w3d.selfe[ic,ix,iy,iz]
  else:
    iz1 = top.izpslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izpslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzpslave[me] + 1
    eee = w3d.selfe[ic,:,:,iz1:iz2]
    if ix is not None and iy is None:
      eee = eee[ix,:,:]
    elif ix is None and iy is not None:
      eee = eee[:,iy,:]
    elif ix is not None and iy is not None:
      eee = eee[ix,iy,:]
    if iz is None:
      eee = transpose(gatherarray(transpose(eee)))
    else:
      pe = convertizptope(iz)
      if pe is None: return None
      if me == pe: eee = eee[...,iz-top.izpslave[me]]
      else:        eee = zeros(shape(eee[...,0]),'d')
      if (me == pe or me == 0) and (pe != 0): eee = getarray(pe,eee,0)
    if bcast: eee = broadcast(eee)

  if not fullplane:
    return eee
  else:
    ii = 0
    if ix is None and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
      ss = array(shape(eee))
      nn = ss[ii] - 1
      ss[ii] = 2*nn + 1
      ee1 = zeros(tuple(ss),'d')
      if comp == 'x': fsign = -1
      else:           fsign = +1
      ee1[nn:,...] = eee
      ee1[nn::-1,...] = fsign*eee
      eee = ee1
    if ix is None: ii = ii + 1
    if iy is None and (w3d.l2symtry or w3d.l4symtry):
      ss = array(shape(eee))
      nn = ss[ii] - 1
      ss[ii] = 2*nn + 1
      ee1 = zeros(tuple(ss),'d')
      if comp == 'y': fsign = -1
      else:           fsign = +1
      if ii == 0:
        ee1[nn:,...] = eee
        ee1[nn::-1,...] = fsign*eee
      else:
        ee1[:,nn:,...] = eee
        ee1[:,nn::-1,...] = fsign*eee
      eee = ee1
    return eee

##########################################################################
##########################################################################
def pcrhozy(ix=None,fullplane=1,lbeamframe=1,**kw):
  """Plots contours of charge density in the Z-Y plane
  - ix=w3d.ix_axis X index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  """
  if ix is None: ix = w3d.ix_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.ymmin,w3d.ymmax)
  settitles("Charge density in z-y plane","Z","Y","ix = "+repr(ix))
  rrr = getrho(ix=ix)
  if me > 0: rrr = zeros((w3d.ny+1,w3d.nzfull+1),'d')
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    rr1 = zeros((2*w3d.ny+1,w3d.nzfull+1),'d')
    rr1[w3d.ny:,:] = rrr
    rr1[w3d.ny::-1,:] = rrr
    rrr = rr1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=transpose(rrr),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhozy.__doc__ = pcrhozy.__doc__ + ppgeneric_doc("z","y")
##########################################################################
def pcrhozx(iy=None,fullplane=1,lbeamframe=1,**kw):
  """Plots contours of charge density in the Z-X plane
  - iy=w3d.iy_axis Y index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  """
  if iy is None: iy = w3d.iy_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.xmmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.xmmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.xmmin,w3d.xmmax)
  settitles("Charge density in z-x plane","Z","X","iy = "+repr(iy))
  rrr = getrho(iy=iy)
  if me > 0: rrr = zeros((w3d.nx+1,w3d.nzfull+1),'d')
  if fullplane and w3d.l4symtry:
    rr1 = zeros((2*w3d.nx+1,w3d.nzfull+1),'d')
    rr1[w3d.nx:,:] = rrr
    rr1[w3d.nx::-1,:] = rrr
    rrr = rr1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=transpose(rrr),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhozx.__doc__ = pcrhozx.__doc__ + ppgeneric_doc("z","x")
##########################################################################
def pcrhoxy(iz=None,fullplane=1,**kw):
  """Plots contours of charge density in the X-Y plane
  - iz=w3d.iz_axis: Z index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
  """
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if not kw.has_key('xmin'): kw['xmin'] = w3d.xmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.xmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
  settitles("Charge density in x-y plane","X","Y","iz = "+repr(iz))
  rrr = getrho(iz=iz)
  if me > 0: rrr = zeros((w3d.nx+1,w3d.ny+1),'d')
  if fullplane and w3d.l4symtry:
    rr1 = zeros((2*w3d.nx+1,2*w3d.ny+1),'d')
    rr1[w3d.nx:,w3d.ny:] = rrr
    rr1[w3d.nx::-1,w3d.ny:] = rrr
    rr1[w3d.nx:,w3d.ny::-1] = rrr
    rr1[w3d.nx::-1,w3d.ny::-1] = rrr
    rrr = rr1
    kw['ymin'] = - kw['ymax']
    kw['xmin'] = - kw['xmax']
  elif fullplane and w3d.l2symtry:
    rr1 = zeros((w3d.nx+1,2*w3d.ny+1),'d')
    rr1[:,w3d.ny:] = rrr
    rr1[:,w3d.ny::-1] = rrr
    rrr = rr1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=rrr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhoxy.__doc__ = pcrhoxy.__doc__ + ppgeneric_doc("x","y")
##########################################################################
def pcphizy(ix=None,fullplane=1,lbeamframe=1,**kw):
  """Plots contours of electrostatic potential in the Z-Y plane
  - ix=w3d.ix_axis X index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  """
  if ix is None: ix = w3d.ix_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.ymmin,w3d.ymmax)
  settitles("Electrostatic potential in z-y plane","Z","Y","ix = "+repr(ix))
  ppp = getphi(ix=ix)
  if me > 0: ppp = zeros((w3d.ny+1,w3d.nzfull+1),'d')
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    pp1 = zeros((2*w3d.ny+1,w3d.nzfull+1),'d')
    pp1[w3d.ny:,:] = ppp
    pp1[w3d.ny::-1,:] = ppp
    ppp = pp1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=transpose(ppp),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphizy.__doc__ = pcphizy.__doc__ + ppgeneric_doc("z","y")
##########################################################################
def pcphizx(iy=None,fullplane=1,lbeamframe=1,**kw):
  """Plots contours of electrostatic potential in the Z-X plane
  - iy=w3d.iy_axis Y index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  """
  if iy is None: iy = w3d.iy_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.xmmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.xmmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.xmmin,w3d.xmmax)
  settitles("Electrostatic potential in z-x plane","Z","X","iy = "+repr(iy))
  ppp = getphi(iy=iy)
  if me > 0: ppp = zeros((w3d.nx+1,w3d.nzfull+1),'d')
  if fullplane and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
    pp1 = zeros((2*w3d.nx+1,w3d.nzfull+1),'d')
    pp1[w3d.nx:,:] = ppp
    pp1[w3d.nx::-1,:] = ppp
    ppp = pp1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=transpose(ppp),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphizx.__doc__ = pcphizx.__doc__ + ppgeneric_doc("z","x")
##########################################################################
def pcphixy(iz=None,fullplane=1,**kw):
  """Plots contours of electrostatic potential in the X-Y plane
  - iz=w3d.iz_axis Z index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
  """
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if not kw.has_key('xmin'): kw['xmin'] = w3d.xmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.xmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('cellarray') or not kw['cellarray']:
    if not kw.has_key('contours'): kw['contours'] = 20
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
  settitles("Electrostatic potential in x-y plane","X","Y","iz = "+repr(iz))
  ppp = getphi(iz=iz)
  if me > 0: ppp = zeros((w3d.nx+1,w3d.ny+1),'d')
  if fullplane and w3d.l4symtry:
    pp1 = zeros((2*w3d.nx+1,2*w3d.ny+1),'d')
    pp1[w3d.nx:,w3d.ny:] = ppp
    pp1[w3d.nx::-1,w3d.ny:] = ppp
    pp1[w3d.nx:,w3d.ny::-1] = ppp
    pp1[w3d.nx::-1,w3d.ny::-1] = ppp
    ppp = pp1
    kw['ymin'] = - kw['ymax']
    kw['xmin'] = - kw['xmax']
  elif fullplane and w3d.l2symtry:
    pp1 = zeros((w3d.nx+1,2*w3d.ny+1),'d')
    pp1[:,w3d.ny:] = ppp
    pp1[:,w3d.ny::-1] = ppp
    ppp = pp1
    kw['ymin'] = - kw['ymax']
  ppgeneric(grid=ppp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphixy.__doc__ = pcphixy.__doc__ + ppgeneric_doc("x","y")
##########################################################################
def pcselfezy(comp='',ix=None,fullplane=1,lbeamframe=1,vec=0,sz=1,sy=1,**kw):
  """Plots contours of electrostatic field in the Z-Y plane
  - comp: field component to plot, either 'x', 'y', or 'z'
  - ix=w3d.ix_axis X index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - vec=0: when true, plots E field vectors
  - sz,sy=1: step size in grid for plotting fewer points
  """
  if ix is None: ix = w3d.ix_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.ymmin,w3d.ymmax)
  settitles("Electrostatic E%s in z-y plane"%comp,"Z","Y","ix = "+repr(ix))
  if fullplane and (w3d.l2symtry or w3d.l4symtry): kw['ymin'] = - kw['ymax']
  if not vec:
    if not kw.has_key('cellarray') or not kw['cellarray']:
      if not kw.has_key('contours'): kw['contours'] = 20
    eee = getselfe(comp=comp,ix=ix,fullplane=fullplane)
    ppgeneric(grid=transpose(eee),kwdict=kw)
  else:
    ey = getselfe(comp='y',ix=ix,fullplane=fullplane)
    ez = getselfe(comp='z',ix=ix,fullplane=fullplane)
    ppvector(transpose(ey[::sy,::sz]),transpose(ez[::sy,::sz]),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcselfezy.__doc__ = pcselfezy.__doc__ + ppgeneric_doc("z","y")
##########################################################################
def pcselfezx(comp=None,iy=None,fullplane=1,lbeamframe=1,vec=0,sz=1,sx=1,**kw):
  """Plots contours of electrostatic potential in the Z-X plane
  - comp: field component to plot, either 'x', 'y', or 'z'
  - iy=w3d.iy_axis Y index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - vec=0: when true, plots E field vectors
  - sz,sx=1: step size in grid for plotting fewer points
  """
  if iy is None: iy = w3d.iy_axis
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if not kw.has_key('xmin'): kw['xmin'] = top.zplmin + zbeam
  if not kw.has_key('xmax'): kw['xmax'] = top.zplmax + zbeam
  if not kw.has_key('ymin'): kw['ymin'] = w3d.xmmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.xmmax
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+zbeam,top.zplmax+zbeam,w3d.xmmin,w3d.xmmax)
  settitles("Electrostatic E%s in z-x plane"%comp,"Z","X","iy = "+repr(iy))
  if fullplane and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
    kw['ymin'] = - kw['ymax']
  if not vec:
    if not kw.has_key('cellarray') or not kw['cellarray']:
      if not kw.has_key('contours'): kw['contours'] = 20
    eee = getselfe(comp=comp,iy=iy,fullplane=fullplane)
    ppgeneric(grid=transpose(eee),kwdict=kw)
  else:
    ex = getselfe(comp='x',iy=iy,fullplane=fullplane)
    ez = getselfe(comp='z',iy=iy,fullplane=fullplane)
    ppvector(transpose(ex[::sx,::sz]),transpose(ez[::sx,::sz]),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcselfezx.__doc__ = pcselfezx.__doc__ + ppgeneric_doc("z","x")
##########################################################################
def pcselfexy(comp=None,iz=None,fullplane=1,vec=0,sx=1,sy=1,**kw):
  """Plots contours of electrostatic potential in the X-Y plane
  - comp: field component to plot, either 'x', 'y', or 'z'
  - iz=w3d.iz_axis Z index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - vec=0: when true, plots E field vectors
  - sx,sy=1: step size in grid for plotting fewer points
  """
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if not kw.has_key('xmin'): kw['xmin'] = w3d.xmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.xmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if kw.has_key('pplimits'):
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
  settitles("Electrostatic E%s in x-y plane"%comp,"X","Y","iz = "+repr(iz))
  if fullplane and w3d.l4symtry:
    kw['ymin'] = - kw['ymax']
    kw['xmin'] = - kw['xmax']
  elif fullplane and w3d.l2symtry:
    kw['ymin'] = - kw['ymax']
  if not vec:
    if not kw.has_key('cellarray') or not kw['cellarray']:
      if not kw.has_key('contours'): kw['contours'] = 20
    eee = getselfe(comp=comp,iz=iz,fullplane=fullplane)
    ppgeneric(grid=eee,kwdict=kw)
  else:
    ex = getselfe(comp='x',iz=iz,fullplane=fullplane)
    ey = getselfe(comp='y',iz=iz,fullplane=fullplane)
    ppvector(ey[::sx,::sy],ex[::sx,::sy],kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcselfexy.__doc__ = pcselfexy.__doc__ + ppgeneric_doc("x","y")
##########################################################################
##########################################################################
def ppdecomposition(scale=1.,minscale=0.,gap=0.2):
  """Shows the domain decomposition in a graphical way. For each
processor, the total mesh extent is plotted as a filled rectangle
covering the z-length and with height determined by 'scale' and the
number of processors. Another filled rectangle is plotted in the top
half showing the particle domains, and one on the lower half shows the
field domain.
  - scale=1.: the maximum vertical extent of the graph
  - minscale=0.: the minimum vertical extent of the graph
  - gap=0.2: fractional vertical gap between rectangles
  """
  z = []
  x = []
  y = []
  dd = 1.*scale/top.maxslaves
  mm = 1. - gap
  for i in xrange(top.maxslaves):
    z = z + [1.]
    zmin = top.zmslmin[i]
    zmax = top.zmslmax[i]
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,mm,mm,-mm]))
  for i in xrange(top.maxslaves):
    z = z + [2.]
    zmin = top.zpslmin[i]
    zmax = top.zpslmax[i]
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([0,0,mm,mm,0]))
  for i in xrange(top.maxslaves):
    z = z + [3.]
    zmin = top.izfsslave[i]*w3d.dz
    zmax = top.izfsslave[i]*w3d.dz + top.nzfsslave[i]*w3d.dz
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,0,0,-mm]))
  plfp(array(z),y,x,5*ones(len(z)),cmin=0,cmax=4)
  for i in xrange(len(z)):
    pldj(x[i*5:i*5+4],y[i*5:i*5+4],x[i*5+1:i*5+5],y[i*5+1:i*5+5])
      

##########################################################################
##########################################################################
plps = []
plfreq = []
plseldom = []
plalways = []
#--------------------------------------------------------------------------
def installplseldom(f):
  "Adds a function to the list of functions called when seldom plots are made"
  plseldom.append(f)
#--------------------------------------------------------------------------
def uninstallplseldom(f):
  """Removes the function from the list of functions called when seldom plots
     are made"""
  if f in plseldom:
    plseldom.remove(f)
  else:
    raise 'Warning: uninstallplseldom: no such function had been installed'
#--------------------------------------------------------------------------
def installplalways(f):
  "Adds a function to the list of functions called when always plots are made"
  plalways.append(f)
#--------------------------------------------------------------------------
def uninstallplalways(f):
  """Removes the function from the list of functions called when always plots
     are made"""
  if f in plalways:
    plalways.remove(f)
  else:
    raise 'Warning: uninstallplalways: no such function had been installed'

##########################################################################
def pltfld3d(fld='phi',freqflag=always):
  """Makes fields plots which have been turned on
     - fld='phi' quantity to plot, either 'phi' or 'rho'
     - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = current_window()
  window(0)
  nwindows = 9
  for i in xrange(nwindows):
    if (top.icrhoxy[i] == freqflag and fld == "rho"): pcrhoxy[i]
    if (top.icrhozx[i] == freqflag and fld == "rho"): pcrhozx[i]
    if (top.icrhozy[i] == freqflag and fld == "rho"): pcrhozy[i]
    if (top.icphixy[i] == freqflag and fld == "phi"): pcphixy[i]
    if (top.icphizx[i] == freqflag and fld == "phi"): pcphizx[i]
    if (top.icphizy[i] == freqflag and fld == "phi"): pcphizy[i]
  #if (top.icrhoxy4 == freqflag and fld == "rho"): pcrhoxy4
  #if (top.icrhozx4 == freqflag and fld == "rho"): pcrhozx4
  #if (top.icrhozy4 == freqflag and fld == "rho"): pcrhozy4
  #if (top.icphixy4 == freqflag and fld == "phi"): pcphixy4
  #if (top.icphizx4 == freqflag and fld == "phi"): pcphizx4
  #if (top.icphizy4 == freqflag and fld == "phi"): pcphizy4
  oldlimits = limits()
  window(currentwindow)

##########################################################################
def onedplts(freqflag=always):
  """Makes 1-D plots which have been turned on
     - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = current_window()
  window(0)
  if freqflag == top.ipcurr: pzcurr()
  if freqflag == top.ipegap: pzegap()
  if freqflag == top.iplchg: pzlchg()
  if freqflag == top.ipvzofz: pzvzofz()
  if freqflag == top.iprhoax: pzrhoax()
  if freqflag == top.ipphiax: pzphiax()
  if freqflag == top.ipezax: pzezax()
  oldlimits = limits()
  window(currentwindow)

# --- Thses are defined for the fortran interface. If WARP is not imported
# --- main, then the functions and the always and seldom parameters will
# --- not be accessible from the fortran call. This way avoids that by
# --- declaring parameterless functions and explicitly adding them to main.
def onedpltsalways():
  onedplts(always)
def onedpltsseldom():
  onedplts(seldom)
__main__.__dict__['onedpltsalways'] = onedpltsalways
__main__.__dict__['onedpltsseldom'] = onedpltsseldom

##########################################################################
def psplots(freqflag=always,js=0):
  """Makes particle phase space plots which have been turned on
     - freqflag=always frequency flag, either always, seldom, or never
     - js=0 specifies the species of particles to plot"""
# --- Phase space plots, both "frequent" ones and others
# --- Do z-x,y 2-to-a-page subset and all-particle plots
  # --- Save current device and set active device to window(0). This
  # --- ensures that plots created by this routine will be dumped to
  # --- the appropriate plot file.
  currentwindow = current_window()
  window(0)

  nsubsets = 3
  nwindows = 9

  for i in xrange(-nsubsets,1):
    if (top.ipzxy[i] == freqflag):
      ppzxy(i,lframe=true)
      fma()

# --- Do z-x,y 2-to-a-page in color, skipping NSKIPCOL particles
  if (top.ipzxyco == freqflag):
    ppzxyco(js,lframe=true)
    fma()

# --- Do z-vz in color, skipping NSKIPCOL particles
  if (top.ipzvzco == freqflag):
    ppzvzco(js,lframe=true)
    fma()

# --- Do x-xp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipxxpco[i] == freqflag):
     ppxxpco(i,lframe=true)
     fma()

# --- Do y-yp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipyypco[i] == freqflag):
     ppyypco(i,lframe=true)
     fma()

# --- Do z-x and z-xp subset and y-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzx[i] == freqflag):
      ppzx(i,lframe=true)
      fma()
  #if (top.ipzx4 == freqflag):
    #ppzx4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzxp[i] == freqflag):
      ppzxp(i,lframe=true)
      fma()
  #if (top.ipzxp4 == freqflag):
    #ppzxp4

# --- Do z-y and z-yp subset and x-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzy[i] == freqflag):
      ppzy(i,lframe=true)
      fma()
  #if (top.ipzy4 == freqflag):
    #ppzy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzyp[i] == freqflag):
      ppzyp(i,lframe=true)
      fma()
  #if (top.ipzyp4 == freqflag):
    #ppzyp4
    #fma()

# --- Do z-vz subset and r-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzvz[i] == freqflag):
      ppzvz(i,lframe=true)
      fma()
  #if (top.ipzvz4 == freqflag):
    #ppzvz4
    #fma()

# --- Do transverse phase-space subset and z-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxy[i] == freqflag):
      ppxy(i,lframe=true)
      fma()
  #if (top.ipxy4 == freqflag):
    #ppxy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxxp[i] == freqflag):
      ppxxp(i,lframe=true)
      fma()
  #if (top.ipxxp4 == freqflag):
    #ppxxp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipyyp[i] == freqflag):
      ppyyp(i,lframe=true)
      fma()
  #if (top.ipyyp4 == freqflag):
    #ppyyp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxpyp[i] == freqflag):
      ppxpyp(i,lframe=true)
      fma()
  #if (top.ipxpyp4 == freqflag):
    #ppxpyp4
    #fma()

# --- Do trace-space z-window plots
  for i in xrange(nwindows+1):
    if (top.iptrace[i] == freqflag and i >= 0):
      pptrace(i,lframe=true)
      fma()

# --- Do the user defined plots
  oldlimits = limits()
  if freqflag == always:
    for p in plalways:
      p()
      fma()
      oldlimits = limits()
    for p in plfreq:
      p()
      fma()
      oldlimits = limits()
  if freqflag == seldom:
    for p in plseldom:
      p()
      fma()
      oldlimits = limits()
    for p in plps:
      p()
      fma()
      oldlimits = limits()

# --- Reset the current window to it previous value.
  oldlimits = limits()
  window(currentwindow)

# --- Thses are defined for the fortran interface. If WARP is not imported
# --- main, then the functions and the always and seldom parameters will
# --- not be accessible from the fortran call. This way avoids that by
# --- declaring parameterless functions and explicitly adding them to main.
def psplotsalways():
  psplots(always)
def psplotsseldom():
  psplots(seldom)
__main__.__dict__['psplotsalways'] = psplotsalways
__main__.__dict__['psplotsseldom'] = psplotsseldom

