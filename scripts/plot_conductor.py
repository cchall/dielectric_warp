from warp import *
plot_conductor_version = "$Id: plot_conductor.py,v 1.25 2002/01/29 22:44:41 dave Exp $"

def plot_conductordoc():
  print """
The following functions plot contours of the potential in various planes
along with the conductors in that plane. The first three plot with the axis
having units of meters. The suffix 'g' means that it plots with the axis
having units of number of grid cells. The 'box' suffix means
that it plots a box around each grid point inside of a conductor.

pfxy, pfzx, pfzy
pfxyg, pfzxg, pfzyg
pfxybox, pfzxbox, pfzybox, pfzxboxi, pfzyboxi

plotgrid: plots the x-z mesh in the lab frame (including any bends)
pfzxlab: makes the pfzx plot in the lab frame (including any bends)
plotsrfrv: handy command to plot r versus z for a suface of revolution, giving
           the function describing it

plotquadoutline: plots outline of quadrupole structure

cleanconductors: not a plot routine, buts removes conductor points not
		 within the the range of the field solve
  """

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotsubgrid(yy,xx,zz,pp,iz,numb,ymin,xmin,dy,dx,color,subgridlen,
                signy,signx):
  assert (pp == 'e' or pp == 'o'),"pp has invalid data"
  nn = eval('f3d.n'+pp+'cndbdy')
  if nn > 0:
    ixc = eval('f3d.i'+pp+'cnd'+xx)*signx
    iyc = eval('f3d.i'+pp+'cnd'+yy)*signy
    izc = eval('f3d.i'+pp+'cnd'+zz)
    delmx = eval('f3d.'+pp+'cdelm'+xx)*signx
    delpx = eval('f3d.'+pp+'cdelp'+xx)*signx
    delmy = eval('f3d.'+pp+'cdelm'+yy)*signy
    delpy = eval('f3d.'+pp+'cdelp'+yy)*signy
    try:
      numbmx = eval('f3d.'+pp+'cnumbm'+xx)
      numbpx = eval('f3d.'+pp+'cnumbp'+xx)
      numbmy = eval('f3d.'+pp+'cnumbm'+yy)
      numbpy = eval('f3d.'+pp+'cnumbp'+yy)
    except:
      numbmx = eval('f3d.'+pp+'cnumb')
      numbpx = eval('f3d.'+pp+'cnumb')
      numbmy = eval('f3d.'+pp+'cnumb')
      numbpy = eval('f3d.'+pp+'cnumb')
  else:
    ixc = array([])
    iyc = array([])
    izc = array([])
    delmx = array([])
    delpx = array([])
    delmy = array([])
    delpy = array([])
    numbmx = array([])
    numbpx = array([])
    numbmy = array([])
    numbpy = array([])
  ii = compress(equal(izc[:nn],iz),arange(nn))
  xx = take(ixc,ii)*dx+xmin
  yy = take(iyc,ii)*dy+ymin
  delmx = take(delmx,ii)*dx
  delpx = take(delpx,ii)*dx
  delmy = take(delmy,ii)*dy
  delpy = take(delpy,ii)*dy
  if numb is not None:
    numbmx = take(numbmx,ii)
    numbpx = take(numbpx,ii)
    numbmy = take(numbmy,ii)
    numbpy = take(numbpy,ii)
  if lparallel:
    xx = gatherarray(xx)
    yy = gatherarray(yy)
    delmx = gatherarray(delmx)
    delpx = gatherarray(delpx)
    delmy = gatherarray(delmy)
    delpy = gatherarray(delpy)
  # --- This code combines all of the individual lines into the list pp.
  # --- This vectorized code avoids slower explicit loops.
  niota = arange(len(xx))
  ii = compress(less(abs(delmy),dy*subgridlen),niota)
  pp = map(lambda x,y,d:[y,y-d,x,x],take(xx,ii),take(yy,ii),take(delmy,ii))
  ii = compress(less(abs(delmx),dx*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y,x,x-d],take(xx,ii),take(yy,ii),take(delmx,ii))
  ii = compress(less(abs(delpy),dy*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y+d,x,x],take(xx,ii),take(yy,ii),take(delpy,ii))
  ii = compress(less(abs(delpx),dx*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y,x,x+d],take(xx,ii),take(yy,ii),take(delpx,ii))
  # --- Convert the list to an array and plot it
  if len(pp) > 0:
    pp = array(pp)
    pldj(pp[:,2],pp[:,0],pp[:,3],pp[:,1],color=color)

# x-y plane
def pfxy(iz=None,izf=None,fullplane=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  """
  kw.update(kwdict)
  # --- This logic is needed since in the parallel version, iz_axis already
  # --- has izslave subtracted from it. If the user passes in a value,
  # --- it must be checked for consistency, otherwise coding below could lead
  # --- to a deadlock in the parallel version
  if izf is not None: iz = izf
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if iz < 0 or w3d.nzfull < iz: return
  izlocal = iz - top.izslave[me]
  if scale:
    dx = w3d.dx
    dy = w3d.dy
    xmmin = w3d.xmmin
    ymmin = w3d.ymmin
    xmmax = w3d.xmmax
    ymmax = w3d.ymmax
  else:
    dx = 1.
    dy = 1.
    xmmin = 0.
    ymmin = 0.
    xmmax = w3d.nx
    ymmax = w3d.ny
  if plotphi:
    kw['xmin'] = xmmin
    kw['xmax'] = xmmax
    kw['ymin'] = ymmin
    kw['ymax'] = ymmax
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphixy,(iz,fullplane),kw)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.izcond[0:f3d.ncond],izlocal),arange(f3d.ncond))
    yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
  else:
    yy = array([])
    xx = array([])
  warpplp(yy,xx,color=condcolor)
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    warpplp(-yy,xx,color=condcolor)
  if fullplane and w3d.l4symtry:
    warpplp( yy,-xx,color=condcolor)
    warpplp(-yy,-xx,color=condcolor)
  if (plotsg):
    plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                subgridlen,1,1)
    plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                subgridlen,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,1,-1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,1,-1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,-1,1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,-1,1)
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,-1,-1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,-1,-1)

# z-x plane
def pfzx(iy=None,iyf=None,fullplane=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  """
  kw.update(kwdict)
  if iyf is not None: iy = iyf
  if iy is None: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if scale:
    dx = w3d.dx
    dz = w3d.dz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin
    xmmax = w3d.xmmax
    zmmax = w3d.zmmax
  else:
    dx = 1.
    dz = 1.
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
    xmmax = w3d.nx
    zmmax = w3d.nz
  if plotphi:
    kw['xmin'] = zmmin
    kw['xmax'] = zmmax
    kw['ymin'] = xmmin
    kw['ymax'] = xmmax
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizx,(iy,fullplane),kw)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iy),arange(f3d.ncond))
    xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    xx = array([])
    zz = array([])
  warpplp(xx,zz,color=condcolor)
  if fullplane and w3d.l4symtry:
    warpplp(-xx,zz,color=condcolor)
  if (plotsg):
    plotsubgrid('x','z','y','e',iy,numb,xmmin,zmmin,dx,dz,evencolor,
                subgridlen,1,1)
    plotsubgrid('x','z','y','o',iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','e',iy,numb,xmmin,zmmin,dx,dz,evencolor,
                  subgridlen,-1,1)
      plotsubgrid('x','z','y','o',iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                  subgridlen,-1,1)

# z-y plane
def pfzy(ix=None,ixf=None,fullplane=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  """
  kw.update(kwdict)
  if ixf is not None: ix = ixf
  if ix is None: ix = w3d.ix_axis
  if ix < 0 or w3d.nx < ix: return
  if scale:
    dy = w3d.dy
    dz = w3d.dz
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
    ymmax = w3d.ymmax
    zmmax = w3d.zmmax
  else:
    dy = 1.
    dz = 1.
    ymmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
    ymmax = w3d.ny
    zmmax = w3d.nz
  if plotphi:
    kw['xmin'] = zmmin
    kw['xmax'] = zmmax
    kw['ymin'] = ymmin
    kw['ymax'] = ymmax
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizy,(ix,fullplane),kw)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ix),arange(f3d.ncond))
    yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    yy = array([])
    zz = array([])
  warpplp(yy,zz,color=condcolor)
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    warpplp(-yy,zz,color=condcolor)
  if (plotsg):
    plotsubgrid('y','z','x','e',ix,numb,ymmin,zmmin,dy,dz,evencolor,
                subgridlen,1,1)
    plotsubgrid('y','z','x','o',ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','z','x','e',ix,numb,ymmin,zmmin,dy,dz,evencolor,
                  subgridlen,-1,1)
      plotsubgrid('y','z','x','o',ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                  subgridlen,-1,1)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(iz=None,izf=None,fullplane=1,plotsg=1,plotphi=1,
          phicolor=blue,subgridlen=1.,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane in grid
frame
Same arguments as pfxy
  """
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-x plane
def pfzxg(iy=None,iyf=None,fullplane=1,plotsg=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane in grid
frame
Same arguments as pfzx
  """
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-y plane
def pfzyg(ix=None,ixf=None,fullplane=1,plotsg=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane in grid
frame
Same arguments as pfzy
  """
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in real units, with inverted x.  Used to make complete plots when  #
# 2-fold symmetry is used.                                           #
######################################################################
 
# x-y plane
def pfxyi(iz=None,izf=None,fullplane=1,plotsg=1,scale=1,plotphi=1,
          phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full X-Y plane,
Same arguments as pfxy
  """
  print "Notice: pfxyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfxy"
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,plotsg=plotsg,scale=scale,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-x plane
def pfzxi(iy=None,iyf=None,fullplane=1,plotsg=1,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-X plane
Same arguments as pfzx
  """
  print "Notice: pfzxi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzx"
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,plotsg=plotsg,scale=scale,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-y plane
def pfzyi(ix=None,ixf=None,fullplane=1,plotsg=1,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-Y plane
Same arguments as pfzy
  """
  print "Notice: pfzyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzy"
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,plotsg=plotsg,scale=scale,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)



############################################################################
# These plot a box at each conductor point
############################################################################

# x-y plane
def pfxybox(iz=None,izf=None,contours=8,plotsg=1,scale=1,signx=1,signy=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - subgridlen=1 maximum length of subgrid line which are plotted
  """
  kw.update(kwdict)
  if izf is not None: iz = izf
  if not iz: iz = w3d.iz_axis
  if iz < 0 or w3d.nzfull < iz: return
  izlocal = iz - top.izslave[me]
  if scale:
    dy = w3d.dy*signy
    dx = w3d.dx*signx
    ymmin = w3d.ymmin
    xmmin = w3d.xmmin
  else:
    dy = 1.*signy
    dx = 1.*signx
    ymmin = 0.
    xmmin = 0.
  if plotphi:
    ppp = getphi(iz=iz)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=xmmin,xmax=xmmin+w3d.nx*dx,
                ymin=ymmin,ymax=ymmin+w3d.ny*dy,kwdict=kw)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.izcond[0:f3d.ncond],izlocal),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
  else:
    x = array([])
    y = array([])
  if lparallel:
    x = gatherarray(x)
    y = gatherarray(y)
  if len(x) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([x-dx/2,x+dx/2,x+dx/2,x-dx/2,x-dx/2]),
        color=condcolor)

# z-x plane
def pfzxbox(iy=None,iyf=None,contours=8,plotsg=1,scale=1,signz=1,signx=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  kw.update(kwdict)
  if iyf is not None: iy = iyf
  if not iy: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  if plotphi:
    ppp = getphi(iy=iy)
    ppp = transpose(ppp)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
                ymin=xmmin,ymax=xmmin+w3d.nx*dx,kwdict=kw)
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iy),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    x = array([])
    z = array([])
  if lparallel:
    x = gatherarray(x)
    z = gatherarray(z)
  if len(x) > 0:
    pla(array([x-dx/2,x-dx/2,x+dx/2,x+dx/2,x-dx/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=condcolor)

# z-y plane
def pfzybox(ix=None,ixf=None,contours=8,plotsg=1,scale=1,signz=1,signy=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  kw.update(kwdict)
  if ixf is not None: ix = ixf
  if not ix: ix = w3d.ix_axis
  if ix < 0 or w3d.nx < ix: return
  if scale:
    dy = w3d.dy*signy
    dz = w3d.dz*signz
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
  else:
    dy = 1.*signy
    dz = 1.*signz
    ymmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  if plotphi:
    ppp = getphi(ix=ix)
    ppp = transpose(ppp)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
                ymin=ymmin,ymax=ymmin+w3d.ny*dy,kwdict=kw)
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ix),arange(f3d.ncond))
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    y = array([])
    z = array([])
  if lparallel:
    y = gatherarray(y)
    z = gatherarray(z)
  if len(y) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=condcolor)

# z-x plane
def pfzxboxi(iy=None,iyf=None,contours=8,plotsg=1,scale=1,signz=1,
             plotphi=1,filled=0,phicolor=blue,condcolor=cyan,**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-X) plane
  - iy=w3d.iy_axis y index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  if iyf is not None: iy = iyf
  pfzxbox(iy=iy,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signx=-1,plotphi=plotphi,filled=filled,
          phicolor=phicolor,condcolor=condcolor,kwdict=kw)

# z-y plane
def pfzyboxi(ix=None,ixf=None,contours=8,plotsg=1,scale=1,signz=1,signy=-1,
             plotphi=1,filled=0,phicolor=blue,condcolor=cyan,**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-Y) plane
  - ix=w3d.ix_axis x index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  if ixf is not None: ix = ixf
  pfzybox(ix=ix,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signy=-1,plotphi=plotphi,filled=filled,
          phicolor=phicolor,condcolor=condcolor,kwdict=kw)




############################################################################
# These plots plot the conductor points colored based on the conductor
# number. The list of colors is input by the user.

# --- convenience function
def findunique(i):
  ii = sort(i)
  result = list(compress(ii[:-1]!=ii[1:],ii[:-1])) + [ii[-1]]
  return result

def plotcondn(iz,nc,cx,cy,cz,cn,dx,dy,xmmin,ymmin,marker,color):
  ncolor = len(color)
  if f3d.ncond > 0:
    ii = compress(equal(cz[0:nc],iz),arange(nc))
    xx = take(cx[0:nc],ii)*dx+xmmin
    yy = take(cy[0:nc],ii)*dy+ymmin
    nn = take(cn[0:nc],ii)
  else:
    xx = array([])
    yy = array([])
    nn = array([])
  nlist = gatherarray(nn)
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    x = compress(equal(nn,nlist[i]),xx)
    y = compress(equal(nn,nlist[i]),yy)
    warpplp(y,x,color=color[i%ncolor],marker=marker)

def pfzxn(iy=None,numbs=None,colors=None,cmarker=point,smarker=circle,
          scale=1,signz=1,signx=1,subgridlen=1.,fullplane=1):
  if iy is None: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if colors is None: colors = color
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  plotcondn(iy,f3d.ncond,f3d.izcond,f3d.ixcond,f3d.iycond,f3d.condnumb,
            dx,dz,xmmin,zmmin,cmarker,color)
  ncolor = len(colors)
  nlist = gatherarray(f3d.ecnumb[:f3d.necndbdy])
  nlist = findunique(nlist)
  nlist.remove(0)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y','e',iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','e',iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,-1,1)
  nlist = gatherarray(f3d.ocnumb[:f3d.nocndbdy])
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y','o',iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','o',iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,-1,1)


############################################################################
# These plot the conductors in laboratory frame, using the tolabfrm routine
# to convert from code frame to lab frame.  There is also a routine to plot
# the computational x-z grid in lab frame.
############################################################################

# --- plot grid in lab frame (including bends)
def plotgrid(zz=None,ii=2,plotcond=1):
  """
Plots Z-X grid in the lab frame (including bends)
  - zz=top.zbeam is the center position
  - ii=2 is the step size in the grid points plotted
    2 means that every other grid line is plotted
  - plotcond=1 when true, plots conductors
  """
  if not zz: zz=top.zbeam
  # --- declare temporary data space, 2 2-D arrays to hold grid coordinates
  xxx = zeros((w3d.nx/ii+1,w3d.nz/ii+1),'d')
  zzz = zeros((w3d.nx/ii+1,w3d.nz/ii+1),'d')
  for iz in xrange(w3d.nz/ii+1):
    xxx[:,iz] = w3d.xmesh[::ii]
  for ix in xrange(w3d.nx/ii+1):
    zzz[ix,:] = w3d.zmmin + iota(0,w3d.nz,ii)*w3d.dz + w3d.zz

  # --- If in a bend, convert the grid data to the lab frame
  if top.linbend:
    # --- reshape arrays to make 1-D arrays to pass to tolabfrm
    nn = int((w3d.nx/ii+1)*(w3d.nz/ii+1))
    xxx.shape = (nn)
    zzz.shape = (nn)

    # --- Convert data to lab frame
    tolabfrm(zz,nn,xxx,zzz)

    # --- Reshape back into 2-D arrays
    xxx.shape = (w3d.nx/ii+1,w3d.nz/ii+1)
    zzz.shape = (w3d.nx/ii+1,w3d.nz/ii+1)

  # --- Make plots
  pla(xxx,zzz,marks=0)
  pla(transpose(xxx),transpose(zzz),marks=0)

  if plotcond: pfzxlab(zz)

# --- Make pfzx plot in lab frame
def pfzxlab(zz=None,iy=None,condcolor=cyan):
  """Plots conductors in Z-X lab frame (including bends)
  - zz=top.zbeam is the center position
  - condcolor=cyan color of conductor points inside conductors
  """
  if not zz: zz=top.zbeam
  if iy is None: iy = w3d.iy_axis
  # --- if zz is not equal to zbeam, then calculate conductors for new location
  if (zz != top.zbeam):
    z = top.zbeam
    g = top.zgrid
    top.zbeam = zz
    top.zgrid = zz
    setlatt()
    fieldsol(1)
  # --- gather conductor data
  if f3d.ncond > 0:
    xxxx=compress(equal(f3d.iycond[0:f3d.ncond],iy),
                        f3d.ixcond[0:f3d.ncond])*w3d.dx+w3d.xmmin
    zzzz=compress(equal(f3d.iycond[0:f3d.ncond],iy),
                        f3d.izcond[0:f3d.ncond])*w3d.dz+w3d.zmmin+zz
    # --- convert to lab frame
    tolabfrm(zz,len(xxxx),xxxx,zzzz)   
    # --- make plot
    plg(xxxx,zzzz,marker='\2',color=condcolor)
  # --- restore original conductor data at zbeam
  if (zz != top.zbeam):
    top.zbeam = z
    top.zgrid = g
    setlatt()
    fieldsol(1)


#####################################################################
def plotsrfrv(srfrv,zmin,zmax,n=1000,color='fg',gridframe=0,rscale=1,zscale=1,
              roff=0,zoff=0,rmin=0.,rmax=top.largepos):
  """Handy function for plotting the r versus z for a surface of revolution
 - srfrv: surface of revolution function to plot
 - zmin,zmax: z range to plot
 - n=1000: number of points to plot
 - color='fg': color of line
 - gridframe=0: when true, plots in grid frame
 - rscale=1: scaling for radius
 - zscale=1: scaling for z
 - roff=0: offset for radius
 - zoff=0: offset for z
 - rmin=0: minimum value of r plotted (before applying rscale and roff)
 - rmax=0: maximum value of r plotted (before applying rscale and roff)
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rr = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrv()
    rr[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rr = (rr)/w3d.dx
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  plg(rscale*rr+roff,zscale*zz+zoff,color=color)


#####################################################################
#####################################################################
def plotelementoutline(color,gridframe,axis,ie,ne,
                       ezs,eze,eap,err,erl,egl,egp,eox,eoy,epa,epr,epw,
                       dpal,dpar):
  """Plots the outline of electrostatic elements
  - color: line color
  - gridframe: when true, make plot in grid coordinates
  - axis: selects axis to plot, either 'x' or 'y'
  """
  if axis == 'x': gpsign = 1
  else:           gpsign = -1
  for i in range(ie,ie+ne):
    # --- plot rods
    # --- If aperture is zero, then this quad is skipped
    rodap = eap[i]
    if erl[i] > 0.:
      rodlen = erl[i]
      gp = egp[i]*gpsign
      gaplen = egl[i]
    else:
      rodlen = (eze[i] - ezs[i])
      gp = 1*gpsign
      gaplen = 0.
    if err[i] > 0.: rodrr = err[i]
    else:           rodrr = 8./7.*eap[i]
    if axis == 'x': offset = eox[i]
    else:           offset = eoy[i]
    if rodap > 0. and rodlen > 0.:
      rr = rodap + rodrr + rodrr*array([1.,1.,-1.,-1.,1.])
      zz = gp*(-0.5*(rodlen+gaplen) + rodlen*array([0.,1.,1.,0.,0.]))
      rr1 = offset + rr
      rr2 = offset - rr
      zz = 0.5*(eze[i] + ezs[i]) + top.zlatstrt + zz
      if gridframe:
        rr1 = rr1/w3d.dx
        rr2 = rr2/w3d.dx
        zz = (zz - w3d.zmmin)/w3d.dz
      plg(rr1,zz,color=color)
      plg(rr2,zz,color=color)
    # --- Plot end plates
    pw = epw[i]
    if pw > 0.:
      if epa[i] > 0.: pa = epa[i]
      else:           pa = eap[i]
      if epr[i] > 0.: pr = epr[i]
      else:           pr = rodap + 2.*rodrr
      pal = pa + dpal[i]
      par = pa + dpar[i]
      rrl = array([pr,pr,pal,pal,pr])
      rrr = array([pr,pr,par,par,pr])
      zz = pw*array([0.,1.,1.,0.,0.])
      rrl1 = offset + rrl
      rrl2 = offset - rrl
      rrr1 = offset + rrr
      rrr2 = offset - rrr
      zzl = 0.5*(eze[i] + ezs[i]) - 0.5*(rodlen+gaplen) - zz + \
            top.zlatstrt
      zzr = 0.5*(eze[i] + ezs[i]) + 0.5*(rodlen+gaplen) + zz + \
            top.zlatstrt
      if gridframe:
        rrl1 = rrl1/w3d.dx
        rrl2 = rrl2/w3d.dx
        rrr1 = rrr1/w3d.dx
        rrr2 = rrr2/w3d.dx
        zzl = (zzl - w3d.zmmin)/w3d.dz
        zzr = (zzr - w3d.zmmin)/w3d.dz
      plg(rrl1,zzl,color=color)
      plg(rrl2,zzl,color=color)
      plg(rrr1,zzr,color=color)
      plg(rrr2,zzr,color=color)


#---------------------------------------------------------------------------
def plotquadoutline(iq=0,nq=None,color='fg',gridframe=0,axis='x'):
  """Plots the outline of quadrupole elements
  - iq=0: starting quad to plot
  - nq=top.nquad+1: number of quads to plot
  - color='fg': line color
  - gridframe=0: when true, make plot in grid coordinates
  - axis='x': selects axis to plot, either 'x' or 'y'
  """
  if nq is None: nq = top.nquad + 1
  plotelementoutline(color,gridframe,axis,iq,nq,
                     top.quadzs,top.quadze,top.quadap,top.quadrr,top.quadrl,
                     top.quadgl,top.quadgp,top.qoffx,top.qoffy,
                     top.quadpa,top.quadpr,top.quadpw,
                     top.qdelpal,top.qdelpar)

#---------------------------------------------------------------------------
def plotemltoutline(ie=0,ne=None,color='fg',gridframe=0,axis='x'):
  """Plots the outline of emlt elements
  - ie=0: starting emlt to plot
  - ne=top.nemlt+1: number of emlts to plot
  - color='fg': line color
  - gridframe=0: when true, make plot in grid coordinates
  - axis='x': selects axis to plot, either 'x' or 'y'
  """
  if ne is None: ne = top.nemlt + 1
  plotelementoutline(color,gridframe,axis,ie,ne,
                     top.emltzs,top.emltze,top.emltap,top.emltrr,top.emltrl,
                     top.emltgl,top.emltgp,top.emltox,top.emltoy,
                     top.emltpa,zeros(top.nemlt+1,'d'),top.emltpw,
                     zeros(top.nemlt+1,'d'),zeros(top.nemlt+1,'d'))



#########################################################################
#########################################################################
def cleanconductors():
  """This routine clears out conductor points which are not within the
range of the field solution, w3d.izfsmin and w3d.izfsmax. This is done
for optimization so that time is not wasted on those points.
  """
  if f3d.ncond > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.izcond[:f3d.ncond]), \
                              less(f3d.izcond[:f3d.ncond],w3d.izfsmax+1)), \
                  arange(f3d.ncond))
    xx = take(f3d.ixcond,ii)
    yy = take(f3d.iycond,ii)
    zz = take(f3d.izcond,ii)
    vv = take(f3d.condvolt,ii)
    f3d.ncond = len(ii)
    f3d.ncondmax = f3d.ncond
    gchange("PSOR3d")
    if f3d.ncond > 0:
      f3d.ixcond[:] = xx
      f3d.iycond[:] = yy
      f3d.izcond[:] = zz
      f3d.condvolt[:] = vv
  if f3d.necndbdy > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.iecndz[:f3d.necndbdy]), \
                              less(f3d.iecndz[:f3d.necndbdy],w3d.izfsmax+1)), \
                  arange(f3d.necndbdy))
    xx = take(f3d.iecndx,ii)
    yy = take(f3d.iecndy,ii)
    zz = take(f3d.iecndz,ii)
    mx = take(f3d.ecdelmx,ii)
    my = take(f3d.ecdelmy,ii)
    mz = take(f3d.ecdelmz,ii)
    px = take(f3d.ecdelpx,ii)
    py = take(f3d.ecdelpy,ii)
    pz = take(f3d.ecdelpz,ii)
    vv = take(f3d.ecvolt,ii)
    f3d.necndbdy = len(ii)
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
    gchange("PSOR3d")
    if f3d.necndbdy > 0:
      f3d.iecndx[:f3d.necndbdy] = xx
      f3d.iecndy[:f3d.necndbdy] = yy
      f3d.iecndz[:f3d.necndbdy] = zz
      f3d.ecdelmx[:f3d.necndbdy] = mx
      f3d.ecdelmy[:f3d.necndbdy] = my
      f3d.ecdelmz[:f3d.necndbdy] = mz
      f3d.ecdelpx[:f3d.necndbdy] = px
      f3d.ecdelpy[:f3d.necndbdy] = py
      f3d.ecdelpz[:f3d.necndbdy] = pz
      f3d.ecvolt[:f3d.necndbdy] = vv
  if f3d.nocndbdy > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.iocndz[:f3d.nocndbdy]), \
                              less(f3d.iocndz[:f3d.nocndbdy],w3d.izfsmax+1)), \
                  arange(f3d.nocndbdy))
    xx = take(f3d.iocndx,ii)
    yy = take(f3d.iocndy,ii)
    zz = take(f3d.iocndz,ii)
    mx = take(f3d.ocdelmx,ii)
    my = take(f3d.ocdelmy,ii)
    mz = take(f3d.ocdelmz,ii)
    px = take(f3d.ocdelpx,ii)
    py = take(f3d.ocdelpy,ii)
    pz = take(f3d.ocdelpz,ii)
    vv = take(f3d.ocvolt,ii)
    f3d.nocndbdy = len(ii)
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
    gchange("PSOR3d")
    if f3d.nocndbdy > 0:
      f3d.iocndx[:f3d.nocndbdy] = xx
      f3d.iocndy[:f3d.nocndbdy] = yy
      f3d.iocndz[:f3d.nocndbdy] = zz
      f3d.ocdelmx[:f3d.nocndbdy] = mx
      f3d.ocdelmy[:f3d.nocndbdy] = my
      f3d.ocdelmz[:f3d.nocndbdy] = mz
      f3d.ocdelpx[:f3d.nocndbdy] = px
      f3d.ocdelpy[:f3d.nocndbdy] = py
      f3d.ocdelpz[:f3d.nocndbdy] = pz
      f3d.ocvolt[:f3d.nocndbdy] = vv

