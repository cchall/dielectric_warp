from warp import *
import mpi
import __main__
warpparallel_version = "$Id: warpparallel.py,v 1.30 2002/07/10 17:42:27 dave Exp $"

top.my_index = me
top.nslaves = npes

#-------------------------------------------------------------------------
# --- This reorganizes the particles for the parallel version after the mesh
# --- has been changed. The compiled routine reorg_particles should probably
# --- be called instead (but it is not callable from python yet).
def reorgparticles():
  lout = 1
  while lout:
    zpartbnd(w3d.zmmax,w3d.zmmin,w3d.dz,top.zgrid)
    lout = 0
    for js in range(top.ns):
      i1 = top.ins[js]-1
      i2 = top.ins[js]+top.nps[js]-1
      ii=compress(not_equal(top.uzp[i1:i2],0.),iota(i1,i2-1))
      if len(ii) > 0:
        if (min(take(top.zp,ii))-top.zgrid <  top.zpslmin[me] or
            max(take(top.zp,ii))-top.zgrid >= top.zpslmax[me]):
          lout = 1
 
     #ip1 = top.ins[js]-1
     #ip2 = top.ins[js]+top.nps[js]-1
     #if ip2 > ip1:
     #  if logical_and(not_equal(top.uzp[ip1:ip1],0.),
     #       logical_or(less(top.zp[ip1:ip1],top.zpslmin[me]+top.zgrid),
     #         greater_equal(top.zp[ip1:ip1],top.zpslmax[me]+top.zgrid))):
     #   lout = 1

    lout = globalmax(lout)

# ---------------------------------------------------------------------------
def gatherallzarray(a,zaxis=0):
  """Gathers and broadcasts the data in a z-array. Each processor contributes
the data from within the particle decomposition region it owns. This works
with any array from the groups Z_arrays and Z_Moments.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Get start and end of particle decomposition region
  iz1 = 0
  if me < npes-1: iz2 = top.izpslave[me+1] - 1 - top.izpslave[me]
  else:           iz2 = w3d.nzfull - top.izpslave[me]
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Gather and broadcast it
  result = gatherarray(a[iz1:iz2+1,...])
  result = broadcast(result)
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result
 
# ---------------------------------------------------------------------------
def scatterallzarray(a,zaxis=0):
  if not lparallel: return a
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Get the appropriate subsection
  result = a[top.izpslave[me]:top.izpslave[me]+top.nzpslave[me] + 1,...]
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result

# ---------------------------------------------------------------------------
def gatherallzfsarray(a,zaxis=0):
  """Gathers and broadcasts the data in a field-solve array. Each
processor contributes the data from within the field-solve decomposition
region it owns.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Get start and end of field-solve decomposition region
  iz1 = 0
  if me < npes-1: iz2 = top.izfsslave[me+1] - 1 - top.izfsslave[me]
  else:           iz2 = w3d.nzfull - top.izfsslave[me]
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Gather and broadcast it
  result = gatherarray(a[iz1:iz2+1,...])
  result = broadcast(result)
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result
 
# ---------------------------------------------------------------------------
def scatterallzfsarray(a,zaxis=0):
  if not lparallel: return a
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Get the appropriate subsection
  result = a[top.izfsslave[me]:top.izfsslave[me]+top.nzfsslave[me] + 1,...]
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result

# ---------------------------------------------------------------------------
def gatherallzmarray(a,zaxis=0):
  """Gathers and broadcasts the data in a z-array. Each processor contributes
the data from within the field decomposition region it owns.
 - first argument is the z-array
 - zaxis: axis which is decomposed in z
  """
  if not lparallel: return a
  # --- Get start and end of particle decomposition region
  iz1 = 0
  if me < npes-1: iz2 = top.izslave[me+1] - 1 - top.izslave[me]
  else:           iz2 = w3d.nzfull - top.izslave[me]
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Gather and broadcast it
  result = gatherarray(a[iz1:iz2+1,...])
  result = broadcast(result)
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result
 
# ---------------------------------------------------------------------------
def scatterallzmarray(a,zaxis=0):
  if not lparallel: return a
  # --- Rearrange array to put the decomposed axis first
  if zaxis != 0: a = swapaxes(a,0,zaxis)
  # --- Get the appropriate subsection
  result = a[top.izslave[me]:top.izslave[me]+top.nzslave[me] + 1,...]
  # --- Rearrange array to put the decomposed axis back where it started
  if zaxis != 0: result = swapaxes(result,0,zaxis)
  return result

#-------------------------------------------------------------------------
def convertiztope(iz):
  """Given an iz value, returns the processor number whose particle region
contains that value."""
  if 0 <= iz <= w3d.nzfull:
    # --- This finds all of the processors for which have iz within their
    # --- domain. The last one is selected since in the regions which
    # --- overlap, the standard is that the processor to the right has
    # --- priority for that region, i.e. the processor which has the
    # --- overlapping on it left hand edge.
    pe = compress(logical_and(less_equal(top.izpslave,iz),
                    less_equal(iz,top.izpslave+top.nzpslave)),arange(npes))[-1]
  else:
    pe = None
  return pe
convertizptope = convertiztope
def convertizfstope(iz):
  """Given an iz value, returns the processor number whose field solve region
contains that value."""
  if 0 <= iz <= w3d.nzfull:
    # --- This finds all of the processors for which have iz within their
    # --- domain. The last one is selected since in the regions which
    # --- overlap, the standard is that the processor to the right has
    # --- priority for that region, i.e. the processor which has the
    # --- overlapping on it left hand edge.
    pe = compress(logical_and(less_equal(top.izfsslave,iz),
                  less_equal(iz,top.izfsslave+top.nzfsslave)),arange(npes))[-1]
  else:
    pe = None
  return pe

def convertiwtope(iw):
  """Given an z window number, returns the processor number whose region
contains that range (actually the center of the window)."""
  if 0 <= iw <= top.nzwind:
    zz = 0.5*(top.zwindows[0,iw]+top.zwindows[1,iw])
    pe = compress(logical_and(less_equal(top.zpslmin,zz),
                                    less(zz,top.zpslmax)),arange(npes))
    if len(pe) > 0:
      return pe[-1]
    else:
      # --- The zwindow is outside of the grid boundary.
      return None
  else:
    # --- The input iw is not a proper value.
    return None

# --- Gathers windows data onto PE0.
def getwin_moments():
  """Broadcasts the window moments data to all processors."""
  # --- First, get a list of all arrays in the group Win_Moments
  vlistall = top.varlist('Win_Moments')
  vlist = []
  for v in vlistall:
    if eval('type(top.'+v+')==type(array([]))'): vlist.append(v)
  # --- Loop over the number of windows and get the pe that owns it.
  # --- The processor is then the root for the broadcast call.
  # --- All of the moment data for each window is sent as a single array.
  for iw in range(1,top.nzwind+1):
    pe = convertiwtope(iw)
    if pe is not None:
      vdata = []
      for v in vlist:
        h = eval('top.'+v)
        vdata.append(h[iw])
      vdata = mpi.bcast(array(vdata),pe)
      i = 0
      for v in vlist:
        h = eval('top.'+v)
        h[iw] = vdata[i]
        i = i + 1

# --- Gathers windows history data onto PE0.
# --- Still need to deal with linechg and hvzofz and other zmoments histories.
def gethist():
  # --- Do nothing if there is no history data
  if top.jhist < 0: return
  # --- All zwindow histories have the attribute winhist.
  vlist = top.varlist('winhist')
  # --- Loop over the number of windows and get the pe that owns it.
  # --- The processor is then the root for the broadcast call.
  # --- All of the history data for each window is sent as a single array.
  for iw in range(1,top.nzwind+1):
    pe = convertiwtope(iw)
    if pe is not None:
      vdata = []
      for v in vlist:
        h = eval('top.'+v)
        vdata.append(h[iw,:top.jhist+1])
      vdata = mpi.bcast(array(vdata),pe)
      for i in range(len(vlist)):
        h = eval('top.'+vlist[i])
        h[iw,:top.jhist+1] = vdata[i]

# --- Gather lab moments onto PE0. Note that after the moments are gathered,
# --- more data can still be calculated if the run continues.
def getlabmoments():
  if not top.iflabwn or top.nlabwn == 0: return
  # --- First, get total number of data points saved.
  i = parallelsum(top.ilabwn)
  # --- Make sure there is space on processor 0
  if me == 0:
    if top.ntlabwn < max(i):
      top.ntlabwn = max(i)
      gchange('Lab_Moments')
  # --- Now gather the data
  varlist = top.varlist('Lab_Moments')
  for v in varlist:
    h = eval('top.'+v)
    if type(h) == type(array([])) and len(shape(h)) == 2:
      for il in range(top.nlabwn):
        gatherh = gatherarray(h[:top.ilabwn[il],il])
        if me == 0:
          h[:i[il],il] = gatherh
  if me == 0:
    top.ilabwn[:] = i
  else:
    top.ilabwn = 0

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Make a restart dump for a parallel simulation. There are three basic steps.
#  1 PE0 writes out all non-parallel data (data which is the same on all
#    processors, such as the lattice)
#  2 For parallel data, PE0 creates entries in the dump file big enough to
#    hold the data from all of the processors
#  3 All of the processors, in parallel, write out the parallel data
# An attempt was made to make the parallel restart file the same format as
# a serial restart dump. This allows a serial job to use a parallel restart
# dump without modification. It also would make is easier for a restart
# dump to be use by jobs with differing numbers of processors (not implemented
# yet). To do that, global values of scalars are written out and arrays are
# written out in the same format as serial dump.
def paralleldump(fname,attr='dump',vars=[],serial=0,histz=2,varsuffix=None,
                 verbose=false):
  getwin_moments()
  gethist()
  getlabmoments()

  # --- Convert attr into a list if needed
  if not (type(attr) == type([])): attr = [attr]

  # --- Gather nps from all processors
  nps_p = gatherarray(top.nps)
  nps_p = mpi.bcast(nps_p,0)
  nps_p.shape = (top.nslaves,top.ns)
  nps_p0 = zeros((top.nslaves+1,top.ns+1))
  nps_p0[1:,1:] = nps_p

  # --- Need boundnz from the right most processor.
  boundnz_p = gatherarray(f3d.boundnz)

  # --- Gather conductor information from all processors
  # --- Also setup temp arrays which have gaurd cells at the lower end.
  # --- These are used to calculate partial sums.
  ncond_p = gatherarray(f3d.ncond)
  ncond_p = mpi.bcast(ncond_p,0)
  ncond_p0 = array([0] + list(ncond_p))
  necndbdy_p = gatherarray(f3d.necndbdy)
  necndbdy_p = mpi.bcast(necndbdy_p,0)
  necndbdy_p0 = array([0] + list(necndbdy_p))
  nocndbdy_p = gatherarray(f3d.nocndbdy)
  nocndbdy_p = mpi.bcast(nocndbdy_p,0)
  nocndbdy_p0 = array([0] + list(nocndbdy_p))

  # --- PE0 first writes out all of the non-parallel stuff and sets up
  # --- space for parallel stuff.
  if me==0:
    # --- PE0 creates the file
    ff = PW.PW(fname)

    # --- Might as well write out *_p data right now while I'm thinking of it.
    # --- But only if they have the attribute attr.
    for a in attr:
      if re.search(a,top.getvarattr("nps")):
        ff.write('nps_p@parallel',nps_p)
      if re.search(a,f3d.getvarattr("ncond")):
        ff.write('ncond_p@parallel',ncond_p)
      if re.search(a,f3d.getvarattr("necndbdy")):
        ff.write('necndbdy_p@parallel',necndbdy_p)
      if re.search(a,f3d.getvarattr("nocndbdy")):
        ff.write('nocndbdy_p@parallel',nocndbdy_p)

    # --- Write out the list of python variables to the file.
    # --- This assumes that all interpreted variables are the same on
    # --- all processors.
    pydump(ff=ff,attr=None,vars=vars,varsuffix=varsuffix,verbose=verbose)

    # --- Loop through all variables, getting the ones with attribute attr
    packagelist = package()
    packagelist.reverse()
    for p in packagelist:
      pkg = eval(p,__main__.__dict__)
      # --- Get list of variables in the package p with attribute attr
      vlist = []
      for a in attr: vlist = vlist + pkg.varlist(a)
      for vname in vlist:
        # --- Check if can get a python object for the variable.
        # --- This check is for dynamic arrays - if array is unallocated,
        # --- the getpyobject routine returns None.
        v = pkg.getpyobject(vname)
        if type(v) == ArrayType and v.typecode() == Complex: continue
        if v is not None:
          # --- Get attributes of the variable
          a = pkg.getvarattr(vname)
          # --- Set name of variable in pdb file
          if varsuffix is None:
            pdbname = vname+'@'+p
          else:
            pdbname = vname+varsuffix
          if verbose: print "writing "+p+"."+vname+" as "+pdbname+" -first pass"
          # --- Check if variable has the attribute 'parallel'
          parallelvar = re.search('parallel',a)
          # --- Method of writing variables differs between parallel and
          # --- non-parallel variables
          if not parallelvar:
            # --- Only PE0 writes out non-parallel variables.
            # --- The value is the same on all processors anyway.
            ff.write(pdbname,v)
          else:
            # --- First, deal with scalars. The values written for the scalars
            # --- are the global values. For many of these, the value on PE0
            # --- is the correct value, for example, zmmin. Exceptions are
            # --- those that differ, for example zmmax and nz. Also, space
            # --- is created in the pdb file so that all of the processors
            # --- can write there own values of each of the paralle scalars.
            # --- That is done with the 'defent' call below.
            # --- Note that parallel scalars are written out during a serial
            # --- dump.
            if type(v) != type(array([])):
              # --- First, deal with exceptions
              if (p == 'top' and vname in ['zlmax','zzmax','zmmntmax']) or \
                 (p == 'w3d' and vname in ['zmmax']):
                ff.write(pdbname,top.zmslmax[-1])
              elif (p == 'top' and vname in ['nzl','nzlmax']):
                if v == 0:
                  ff.write(pdbname,0)
                else:
                  ff.write(pdbname,w3d.nzfull)
              elif (p == 'top' and vname in ['nzzarr','nzmmnt']) or \
                   (p == 'w3d' and vname in ['nz','izfsmax','nz_selfe']):
                ff.write(pdbname,w3d.nzfull)
              elif (p=='top' and vname in ['np','nplive','npmax','npmaxb']) or\
                   (p=='wxy' and vname in ['npmaxxy']):
                ff.write(pdbname,sum(nps_p)[0])
              elif (p=='f3d' and vname in ['boundnz']):
                ff.write(pdbname,boundnz_p[-1])
              elif (p=='f3d' and vname in ['ncond','ncondmax']):
                ff.write(pdbname,sum(ncond_p))
              elif (p=='f3d' and vname in ['necndbdy']):
                ff.write(pdbname,sum(necndbdy_p))
              elif (p=='f3d' and vname in ['nocndbdy']):
                ff.write(pdbname,sum(nocndbdy_p))
              elif (p=='f3d' and vname in ['ncndmax']):
                ff.write(pdbname,max(sum(necndbdy_p),sum(nocndbdy_p)))
              else:
                # --- Otherwise, write out value as is on PE0
                ff.write(pdbname,v)
              # --- For all parallel scalars, create entry to write data too
              if not serial:
                ff.defent(vname+'@parallel',array([v]),(top.nslaves,))

            elif not serial:
              # --- Now arrays...
              # --- No parallel arrays are written out during a serial dump.
              # --- The majority of the arrays are 1-D arrays which are
              # --- domain decomposed in z. For these, an entry is made in
              # --- the file big enough to fit the data from all of the
              # --- processors. Other arrays have special requirements.
              # --- First deal with the exceptions
              if p == 'top' and vname in ['sm','sw','sq']:
                # --- These arrays are not actually parallel
                ff.write(pdbname,v)
              elif (p == 'top' and vname == 'zwindows'):
                # --- Window 0 needs to set to cover the full system
                zwin = top.zwindows + 0.
                zwin[1,0] = top.zmslmax[-1]
                ff.write(pdbname,zwin)
                ff.defent(vname+'@parallel',v,(top.nslaves,2,10))
              elif (p == 'top' and vname in ['zmmnts0','zmmnts']) or \
                   (p == 'f3d' and vname in ['sorerrar','boundarr']):
                # --- These arrays don't need to be saved.
                pass
              elif vname == 'npmax_s' and p == 'top':
                # --- This is set to be correct globally
                ff.write(pdbname,array([0]+list(cumsum(sum(nps_p[:,:])))))
                ff.defent(vname+'@parallel',v,(top.nslaves,top.ns+1))
              elif vname == 'ins' and p == 'top':
                # --- This is set to be correct globally
                iii = array([1])
                if top.ns > 1:
                  iii = array([1]+list(cumsum(sum(nps_p[:,:-1]))+array([1])))
                ff.write(pdbname,iii)
                ff.defent(vname+'@parallel',v,(top.nslaves,top.ns))
              elif vname == 'nps' and p == 'top':
                # --- This is set to be correct globally
                ff.write(pdbname,sum(nps_p))
                ff.defent(vname+'@parallel',v,(top.nslaves,top.ns))
              elif p == 'top' and vname in ['xp','yp','zp','uxp','uyp','uzp', \
                                            'gaminv']:
                # --- For the particle data, a space big enough to hold
                # --- all of the data is created.
                if sum(sum(nps_p)) > 0:
                  ff.defent(pdbname,v,(sum(sum(nps_p)),))
              elif p == 'wxy' and vname in ['dtp']:
                # --- A WARPxy particle array
                if wxy.npmaxxy > 0 and sum(sum(nps_p)) > 0:
                  ff.defent(pdbname,v,(sum(sum(nps_p)),))
              elif re.search('zhist',a):
                # --- z moments histories
                ff.write('histz@parallel',histz)
                if not histz:
                  # --- This is a temporary solution.
                  # --- This method can lead to erroneaously large
                  # --- dump files if the decomposition is not nearly
                  # --- uniform.
                  ff.defent(vname+'@parallel',v,
                            (top.nslaves,max(1+top.nzpslave),top.lenhist+1))
                elif histz == 1:
                # --- This would be the proper way, but writing out the
                # --- would be expensive since it would have to be done
                # --- in small chunks.
                  ff.defent(pdbname,v,(w3d.nzfull+1,top.lenhist+1))
                elif histz == 2:
                # --- The tranpose of the array is written out. It is written
                # --- out as one chunk but is in the wrong order.
                  ff.defent(vname+'@parallel',v,(top.lenhist+1,w3d.nzfull+1))
              elif p == 'f3d' and vname in ['ixcond','iycond','izcond', \
                                            'condvolt','icondlxy','icondlz']:
                # --- The conductor data is gathered into one place
                if sum(ncond_p) > 0:
                  ff.defent(pdbname,v,(sum(ncond_p),))
              elif p == 'f3d' and \
                   vname in ['ecndpvph','iecndx','iecndy','iecndz',
                             'ecdelmx','ecdelmy','ecdelmz','ecdelpx',
                             'ecdelpy','ecdelpz','ecvolt','ecvoltmx',
                             'ecvoltpx','ecvoltmy','ecvoltpy','ecvoltmz',
                             'ecvoltpz','iecndlxy','iecndlz']:
                # --- The conductor data is gathered into one place
                if sum(necndbdy_p) > 0:
                  ff.defent(pdbname,v,(sum(necndbdy_p),))
              elif p == 'f3d' and \
                   vname in ['ocndpvph','iocndx','iocndy','iocndz',
                             'ocdelmx','ocdelmy','ocdelmz','ocdelpx',
                             'ocdelpy','ocdelpz','ocvolt','ocvoltmx',
                             'ocvoltpx','ocvoltmy','ocvoltpy','ocvoltmz',
                             'ocvoltpz','iocndlxy','iocndlz']:
                # --- The conductor data is gathered into one place
                if sum(nocndbdy_p) > 0:
                  ff.defent(pdbname,v,(sum(nocndbdy_p),))
              elif p == 'w3d' and vname in ['rho']:
                # --- Be prepared to dump out rho in case it is needed.
                # --- For example the egun script wants rho saved.
                ff.defent(pdbname,v,(w3d.nx+1,w3d.ny+1,w3d.nzfull+1))
              elif p == 'w3d' and vname in ['phi']:
                # --- Be prepared to dump out phi in case it is needed.
                ff.defent(pdbname,v,
                                 (w3d.nx+1,w3d.ny+1,w3d.nzfull+2+w3d.izextra))
              else:
                # --- The rest are domain decomposed Z arrays
                ff.defent(pdbname,v,(w3d.nzfull+1,))

    # --- PE0 closes the file at this point
    ff.close()

  # --- None of the processors can procede past this point until PE0 has
  # --- completed the above.
  mpi.barrier()

  # --- If only writing non-parallel data, then return here
  if serial: return

  # --- Now that that is all done, the parallel data can actually be written
  # --- out now.

  # --- All of the processors open the file for appending
  ff = PW.PW(fname,'a')

  # --- Now we gotta go through all of this mess again!
  # --- Loop through all variables, getting the ones with attribute attr
  # --- See comments above for more details.
  for p in package():
    pkg = eval(p,__main__.__dict__)
    vlist = []
    for a in attr: vlist = vlist + pkg.varlist(a)
    for vname in vlist:
      v = pkg.getpyobject(vname)
      if v is not None:
        if verbose: print "writing "+p+"."+vname+" as "+pdbname+" -second pass"
        a = pkg.getvarattr(vname)
        if varsuffix is None:
          pdbname = vname+'@'+p
        else:
          pdbname = vname+varsuffix
        parallelvar = re.search('parallel',a)
        if parallelvar:
          # --- First, deal with scalars. Each processor writes out it's
          # --- own value of the scalars into a special array.
          if type(v) != type(array([])):
            ff.write(vname+'@parallel',array([v]),indx=(me,))
            # --- That was easy.

          else:
            # --- Now arrays...
            # --- The data is now written into the space which was set aside
            # --- above by PE0.
            # --- First the exceptions
            if (p == 'top' and vname in ['sm','sw','sq','zmmnts0','zmmnts']) or\
               (p == 'f3d' and vname in ['sorerrar','boundarr']):
              # --- Nothing need be done for these.
              pass
            elif (p == 'top' and vname == 'zwindows'):
              ff.write(vname+'@parallel',array([v]),indx=(me,0,0))
            elif p == 'top' and vname in ['npmax_s','ins','nps']:
              # --- Write out to parallel space
              ff.write(vname+'@parallel',array([v]),indx=(me,0))
            elif (p == 'top' and vname in ['xp','yp','zp','uxp','uyp','uzp', \
                                          'gaminv']) or \
                 (p == 'wxy' and vname in ['dtp']):
              # --- Write out each species seperately.
              for js in xrange(top.ns):
                if top.nps[js] > 0:
                  ipmin = sum(sum(nps_p0[:,0:js+1])) + sum(nps_p0[:me+1,js+1])
                  ff.write(pdbname,v[top.ins[js]-1:top.ins[js]+top.nps[js]-1],
                           indx=(ipmin,))
            elif re.search('zhist',a):
              # --- z moments histories
              if not histz:
                # --- Write out the data as one chunk.
                ff.write(vname+'@parallel',array([v]),indx=(me,0,0))
              elif histz == 1:
                # --- Write out the data in proper order. The loop
                # --- over the second dimension of the arrays
                # --- writes out the data for each time step seperately.
                # --- This is VERY slow.
                for ih in xrange(top.jhist+1):
                  ff.write(pdbname,transpose(array([v[:,ih]])),
                           indx=(top.izpslave[me],ih))
              elif histz == 2:
                # --- Write out the transpose of the array.
                ff.write(vname+'@parallel',transpose(v),
                         indx=(0,top.izpslave[me]))
            elif p == 'f3d' and vname in ['ixcond','iycond','izcond', \
                                          'condvolt','icondlxy','icondlz']:
              # --- Write out conductor data.
              offset = 0
              if vname=='izcond': offset = top.izslave[me]
              if f3d.ncond > 0:
                ff.write(pdbname,v[:f3d.ncond]+offset,
                         indx=(sum(ncond_p0[:me+1]),))
            elif p == 'f3d' and \
                 vname in ['ecndpvph','iecndx','iecndy','iecndz',
                           'ecdelmx','ecdelmy','ecdelmz','ecdelpx',
                           'ecdelpy','ecdelpz','ecvolt','ecvoltmx',
                           'ecvoltpx','ecvoltmy','ecvoltpy','ecvoltmz',
                           'ecvoltpz','iecndlxy','iecndlz']:
              # --- Write out conductor data.
              offset = 0
              if vname=='iecndz': offset = top.izslave[me]
              if f3d.necndbdy > 0:
                ff.write(pdbname,v[:f3d.necndbdy]+offset,
                         indx=(sum(necndbdy_p0[:me+1]),))
            elif p == 'f3d' and \
                 vname in ['ocndpvph','iocndx','iocndy','iocndz',
                           'ocdelmx','ocdelmy','ocdelmz','ocdelpx',
                           'ocdelpy','ocdelpz','ocvolt','ocvoltmx',
                           'ocvoltpx','ocvoltmy','ocvoltpy','ocvoltmz',
                           'ocvoltpz','iocndlxy','iocndlz']:
              # --- Write out conductor data.
              offset = 0
              if vname=='iocndz': offset = top.izslave[me]
              if f3d.nocndbdy > 0:
                ff.write(pdbname,v[:f3d.nocndbdy]+offset,
                         indx=(sum(nocndbdy_p0[:me+1]),))
            elif p == 'w3d' and vname in ['rho']:
              iz1 = top.izfsslave[me] - top.izslave[me]
              if me < npes-1: iz2 = top.izfsslave[me+1] - top.izslave[me]
              else:           iz2 = iz1 + top.nzfsslave[me] + 1
              ppp = w3d.rho[:,:,iz1:iz2]
              ff.write(pdbname,ppp,indx=(0,0,top.izfsslave[me]))
            elif p == 'w3d' and vname in ['phi']:
              iz1 = top.izfsslave[me] - top.izslave[me]
              if me == 0: iz1 = iz1 - 1
              if me < npes-1: iz2 = top.izfsslave[me+1] - top.izslave[me]
              else:           iz2 = iz1 + top.nzfsslave[me] + 1
              ppp = w3d.phi[:,:,iz1+1:iz2+1]
              if me == 0: izmin = 0
              else:       izmin = top.izslave[me]+1
              ff.write(pdbname,ppp,indx=(0,0,izmin))
            else:
              # --- The rest are domain decomposed Z arrays
              # --- Assuming they have the same decomposition as the particles.
              iz1 = 0
              iz2 = top.nzpslave[me] + 1
              ff.write(pdbname,v[iz1:iz2],indx=(top.izpslave[me],))

  # --- Everybody closes the file at this point
  ff.close()

  # --- The parallel dump is done.

############################################################################
# This assumes that the dump was made using the paralleldump routine
# with the default attribute, 'dump'.
# There are two main steps
#  1 Read in some initial data which is needed to setup array sizes or
#    needed in the second step
#  2 Read the rest of the data in
#
def parallelrestore(fname,verbose=false,skip=[]):
  # --- All PE's open the file for reading.
  ff = PR.PR(fname)

  # --- Get list of all of the variables in the restart file
  vlist = ff.inquire_ls()

  # --- Remove skipped variables from vlist
  vlistcopy = 1*vlist
  for v in vlistcopy:
    if v in skip:
      vlist.remove(v)
      continue
    if len(v) > 4 and v[-4] == '@':
      if v[:-4] in skip or v[-3:]+'.'+v[:-4] in skip: vlist.remove(v)
  del vlistcopy

  # --- First, setup some arrays that need special handling.
  # --- The particles arrays are returned to there original size as in
  # --- the simulation which made the restart dump. The following
  # --- variables are needed to get that setup correctly.
  # --- In all cases, check if the variable was written out first.
  if 'ns@top' in vlist:
    top.ns = ff.read('ns@top')
  if 'nps_p@parallel' in vlist:
    nps_p = ff.read('nps_p@parallel')
    nps_p0 = zeros((top.nslaves+1,top.ns+1))
    nps_p0[1:,1:] = nps_p
  itriple = array([me,me,1])
  if 'npmax@parallel' in vlist:
    top.npmax = ff.read_part('npmax@parallel',itriple)[0]
  if 'npmaxb@parallel' in vlist:
    top.npmaxb = ff.read_part('npmaxb@parallel',itriple)[0]
  if 'np@parallel' in vlist:
    top.np = ff.read_part('np@parallel',itriple)[0]
  if 'npmaxxy@parallel' in vlist:
    wxy.npmaxxy = ff.read_part('npmaxxy@parallel',itriple)[0]
  gchange("Particles")
  gchange("Particlesxy")
  itriple = array([me,me,1,0,top.ns-1,1])
  if 'ins@parallel' in vlist:
    top.ins[:] = ff.read_part('ins@parallel',itriple)[0,...]
  if 'nps@parallel' in vlist:
    top.nps[:] = ff.read_part('nps@parallel',itriple)[0,...]

  # --- These arrays need to be read in to get the indices for the
  # --- correct data for each processor.
  if 'ncond_p@parallel' in vlist:
    ncond_p = ff.read('ncond_p@parallel')
    ncond_p0 = array([0] + list(ncond_p))
    f3d.ncond = ncond_p[me]
  if 'necndbdy_p@parallel' in vlist:
    necndbdy_p = ff.read('necndbdy_p@parallel')
    necndbdy_p0 = array([0] + list(necndbdy_p))
    f3d.necndbdy = necndbdy_p[me]
  if 'nocndbdy_p@parallel' in vlist:
    nocndbdy_p = ff.read('nocndbdy_p@parallel')
    nocndbdy_p0 = array([0] + list(nocndbdy_p))
    f3d.nocndbdy = nocndbdy_p[me]

  # --- These are needed below and, since there is no gaurantee about the
  # --- order in which variables are read in, they must be set ahead of time.
  # --- Note that the non parallel ones here are actually reread in again
  # --- in the loop below.
  if 'nx@w3d' in vlist:
    w3d.nx = ff.read('nx@w3d')
  if 'ny@w3d' in vlist:
    w3d.ny = ff.read('ny@w3d')
  if 'nz@parallel' in vlist:
    w3d.nz = ff.read_part('nz@parallel',array([me,me,1]))[0]
  if 'izslave@top' in vlist:
    top.forceassign('izslave',ff.__getattr__('izslave@top'))
  if 'lenhist@top' in vlist:
    top.lenhist = ff.read('lenhist@top')
  if 'jhist@top' in vlist:
    top.jhist = ff.read('jhist@top')

  # --- Loop over the list of all of the variables in the restart file.
  # --- Read in all of the scalars first - this ensures that all of the
  # --- integers which describe the size of dynamics arrays are read in
  # --- before the arrays, a requirement of the f90 version.
  for v in vlist:
    if len(v) > 4 and v[-4]=='@':
      # --- Variable is a fortran variable
      vname = v[:-4]
      p = v[-3:]
      pkg = eval(p,__main__.__dict__)
      pname = p+'.'+vname
      # --- The shape is used determine whether the variable is an array
      # --- or not. When the length of the shape is zero, then the
      # --- variable is a scalar.
      if len(ff.inquire_shape(v)) != 0: continue
      if verbose: "reading "+p+"."+vname
      # --- Make sure that the variable is still valid. If not
      # --- (e.g. it has been deleted) then don't try to restore it.
      try:
        a = pkg.getvarattr(vname)
      except pybasisC.error:
        print "Warning: There was a problem %s - it can't be found."%(pname)
      parallelvar = re.search('parallel',a)
      if not parallelvar:
        # --- Simply read variable directly in.
        s = pname+'=ff.__getattr__(v)'
      else:
        # --- Scalars: get data saved with parallel suffix.
        itriple = array([me,me,1])
        s = pname+'= ff.read_part(vname+"@parallel",itriple)[0]'
      try:
        exec(s,__main__.__dict__,locals())
      except:
        print "Warning: There was a problem restoring %s"%(pname)
    elif v[-7:] == '@global':
      try:
        if verbose: print "reading python variable "+v[:-7]
        exec('%s=ff.__getattr__("%s");__main__.__dict__["%s"]=%s'%
             (v[:-7],v,v[:-7],v[:-7]))
      except:
        pass
    elif v[-9:] == '@parallel':
      # --- Do not directly read in variables with the @parallel suffix
      pass
    else:
      try:
        if verbose: print "reading python variable "+v
        exec('%s=ff.%s;__main__.__dict__["%s"]=%s'%(v,v,v,v))
      except:
        pass
  # --- Now read in the arrays.
  for v in vlist:
    if len(v) > 4 and v[-4]=='@':
      # --- Variable is a fortran variable
      vname = v[:-4]
      p = v[-3:]
      pkg = eval(p,__main__.__dict__)
      pname = p+'.'+vname
      # --- The shape is used determine whether the variable is an array
      # --- or not. When the length of the shape is zero, then the
      # --- variable is a scalar.
      if len(ff.inquire_shape(v)) == 0: continue
      if verbose: "reading "+p+"."+vname
      # --- Make sure that the variable is still valid. If not
      # --- (e.g. it has been deleted) then don't try to restore it.
      try:
        a = pkg.getvarattr(vname)
      except pybasisC.error:
        print "Warning: There was a problem %s - it can't be found."%(pname)
      parallelvar = re.search('parallel',a)
      if not parallelvar:
        # --- Simply read variable directly in.
        s = p+'.forceassign(vname,ff.__getattr__(v))'
      else:
        # --- Many arrays need special handling. These are dealt with first.
        if p == 'top' and vname in ['sm','sw','sq']:
          # --- These arrays are not actually parallel and so can just
          # --- be read in.
          s = p+'.forceassign(vname,ff.__getattr__(v))'
        elif (p == 'top' and vname in ['zmmnts0','zmmnts']) or \
             (p == 'f3d' and vname in ['sorerrar','boundarr']):
          # --- These arrays don't need to be restored.
          s = 'pass'
        elif vname == 'zwindows' and p == 'top':
          # --- This doesn't want to work
          #itriple = array([me,me,1,0,1,1,0,9,1])
          #s = p+'.forceassign(vname,\
          #       ff.read_part(vname+"@parallel",itriple)[0,...])'
          zwin = ff.read(vname+"@parallel")
          top.zwindows[:,:] = zwin[me,...]
          s = 'pass'
        elif vname == 'npmax_s' and p == 'top':
          itriple = array([me,me,1,0,top.ns,1])
          s = p+'.forceassign(vname,\
                 ff.read_part(vname+"@parallel",itriple)[0,...])'
        elif p == 'top' and vname in ['ins','nps']:
          # --- These have already been restored above since they are
          # --- needed to read in the particles.
          s = 'pass'
        elif (p == 'top' and vname in ['xp','yp','zp','uxp','uyp','uzp', \
                                      'gaminv']) or \
             (p == 'wxy' and vname in ['dtp']):
          # --- Read in each species seperately.
          # --- The assumption is made that if wxy.dtp was written out,
          # --- it has the same shape as the other particle arrays.
          # --- The command is exec'ed here since a different command
          # --- is needed for each species.  Errors are not caught.
          s = 'pass'
          for js in xrange(top.ns):
            if top.nps[js] > 0:
              ipmin = sum(sum(nps_p0[:,0:js+1])) + sum(nps_p0[:me+1,js+1])
              itriple = array([ipmin,ipmin+top.nps[js]-1,1])
              ip = '[top.ins[js]-1:top.ins[js]+top.nps[js]-1]'
              exec(pname+ip+' = ff.read_part(v,itriple)',
                   __main__.__dict__,locals())
        elif re.search('zhist',a):
          # --- z moments histories
          try:
            histz = ff.read("histz@parallel")
          except:
            histz = 2
          if histz == 0:
            # --- This is a temporary but fast solution.
            itriple = array([me,me,1,0,top.nzpslave[me],1,0,top.lenhist,1])
            s = p+'.forceassign(vname,\
                  ff.read_part(vname+"@parallel",itriple)[0,...])'
          elif histz == 1:
            # --- The proper solution would have to read the data in
            # --- in chunks for each time step.
            itriple = array([top.izpslave[me],
                             top.izpslave[me]+top.nzpslave[me],1,0,0,0])
            tmp = zeros((1+top.nzpslave[me],1+top.lenhist),'d')
            for ih in xrange(top.jhist+1):
              itriple[3:] = [ih,ih,1]
              tmp[:,ih] = ff.read_part(vname+"@parallel",itriple)
            s = p+'.forceassign(vname,tmp)'
          elif histz == 2:
            # --- Read in data and untranspose it.
            itriple = array([0,top.lenhist,1,
                        top.izpslave[me],top.izpslave[me]+top.nzpslave[me],1])
            tmp = ff.read_part(vname+"@parallel",itriple)
            s = p+'.forceassign(vname,transpose(tmp))'
        elif p == 'f3d' and vname in ['ixcond','iycond','izcond', \
                                      'condvolt','icondlxy','icondlz']:
          # --- The conductor data was gathered into one place
          if ncond_p[me] > 0:
            imin = sum(ncond_p0[:me+1])
            itriple = array([imin,imin+ncond_p[me]-1,1])
            s = p+'.forceassign(vname,ff.read_part(v,itriple))'
            if vname == 'izcond':
              s = s + ';f3d.izcond[:]=f3d.izcond[:]-top.izslave[me]'
          else:
            s = 'pass'
        elif p == 'f3d' and \
             vname in ['ecndpvph','iecndx','iecndy','iecndz',
                       'ecdelmx','ecdelmy','ecdelmz','ecdelpx',
                       'ecdelpy','ecdelpz','ecvolt','ecvoltmx',
                       'ecvoltpx','ecvoltmy','ecvoltpy','ecvoltmz',
                       'ecvoltpz','iecndlxy','iecndlz']:
          # --- The conductor data was gathered into one place
          if necndbdy_p[me] > 0:
            imin = sum(necndbdy_p0[:me+1])
            itriple = array([imin,imin+necndbdy_p[me]-1,1])
            s = p+'.forceassign(vname,ff.read_part(v,itriple))'
            if vname == 'iecndz':
              s = s + ';f3d.iecndz[:]=f3d.iecndz[:]-top.izslave[me]'
          else:
            s = 'pass'
        elif p == 'f3d' and \
             vname in ['ocndpvph','iocndx','iocndy','iocndz',
                       'ocdelmx','ocdelmy','ocdelmz','ocdelpx',
                       'ocdelpy','ocdelpz','ocvolt','ocvoltmx',
                       'ocvoltpx','ocvoltmy','ocvoltpy','ocvoltmz',
                       'ocvoltpz','iocndlxy','iocndlz']:
          # --- The conductor data was gathered into one place
          if nocndbdy_p[me] > 0:
            imin = sum(nocndbdy_p0[:me+1])
            itriple = array([imin,imin+nocndbdy_p[me]-1,1])
            s = p+'.forceassign(vname,ff.read_part(v,itriple))'
            if vname == 'iocndz':
              s = s + ';f3d.iocndz[:]=f3d.iocndz[:]-top.izslave[me]'
          else:
            s = 'pass'
        elif p == 'w3d' and vname in ['rho']:
          itriple = array([0,w3d.nx,1,0,w3d.ny,1,
                      top.izpslave[me],top.izpslave[me]+top.nzpslave[me],1])
          s = p+'.forceassign(vname,ff.read_part(v,itriple))'
        elif p == 'w3d' and vname in ['phi']:
          itriple = array([0,w3d.nx,1,0,w3d.ny,1,
              top.izfsslave[me]-1+1,top.izfsslave[me]+top.nzfsslave[me]+2,1])
          s = p+'.forceassign(vname,ff.read_part(v,itriple))'
        else:
          # --- The rest are domain decomposed Z arrays
          itriple = array([top.izpslave[me],
                           top.izpslave[me]+top.nzpslave[me],1])
          s = p+'.forceassign(vname,ff.read_part(v,itriple))'

      try:
        exec(s,__main__.__dict__,locals())
      except:
        print "Warning: There was a problem restoring %s"%(pname)

  ff.close()


##############################################################################
# --- The following dump routines create a file for each processor. This
# --- is faster than the above version but create many files.
# --- Make dump.
def makedump():
  fname = arraytostr(top.runid) + "%06d_%03d.dump" % (top.it, me)
  print "Dumping to file "+fname
  basisdump("dump",fname)
  print "Dump successful"

# --- Read dump back in.
def readdump(fnameprefix):
  fname = fnameprefix + "_%03d.dump" % me
  print "Restore from file "+fname
  basisrestore("dump",fname)
  print "Restore successful"

