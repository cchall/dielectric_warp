#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys,os,os.path,string
from Forthon.compilers import FCompiler
import getopt

try:
    import distutils
    from distutils.core import setup, Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
except:
    raise SystemExit, "Distutils problem"

optlist,args = getopt.getopt(sys.argv[1:],'gt:F:',
                             ['parallel','with-numpy',
                              'fargs=','cargs=','fcompexec='])
machine = sys.platform
debug   = 0
fcomp   = None
parallel = 0
with_numpy = 1
for o in optlist:
  if   o[0] == '-g': debug = 1
  elif o[0] == '-t': machine = o[1]
  elif o[0] == '-F': fcomp = o[1]
  elif o[0] == '--parallel': parallel = 1
  elif o[0] == '--with-numpy': with_numpy = 1

sys.argv = ['setup.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)

dummydist = Distribution()
dummybuild = build(dummydist)
dummybuild.finalize_options()
builddir = dummybuild.build_temp

warppkgs = ['top','env','w3d','f3d','wxy','fxy','wrz','frz','her','cir','cho','em2d','em3d']

def makeobjects(pkg):
  return [pkg+'.o',pkg+'_p.o',pkg+'pymodule.o']

warpobjects = []
for pkg in warppkgs:
  warpobjects = warpobjects + makeobjects(pkg)

warpobjects = warpobjects + ['top_lattice.o','dtop.o',
                             'dw3d.o','w3d_injection.o','w3d_interp.o',
                             'w3d_collisions.o','w3d_utilities.o','w3d_load.o',
                             'f3d_mgrid.o','f3d_ImplicitES.o','f3d_mgrid_be.o',
                             'f3d_conductors.o','f3d_bfield.o',
                             'fft.o','util.o',
                             'fxy_mgrid.o',
                             'dwrz.o',
                             'frz_mgrid.o','frz_mgrid_be.o','frz_ImplicitES.o',
                             'em2d_apml.o','em2d_apml_cummer.o',
                             'em2d_maxwell.o','em3d_maxwell.o']
if parallel:
  warpobjects = warpobjects + ['f3dslave.o','frzslave.o','topslave.o',
                               'w3dslave.o']

warpobjects = map(lambda p:os.path.join(builddir,p),warpobjects)

library_dirs = fcompiler.libdirs
libraries = fcompiler.libs
if parallel:
  library_dirs = library_dirs + ['/usr/lpp/ppe.poe/lib']
  libraries = fcompiler.libs + ['mpi']
  #warpobjects = warpobjects + ['/usr/local/mpi/ifc_farg.o']

include_dirs = [builddir]
if with_numpy:
  import numpy
  define_macros = []
  include_dirs.append(numpy.get_include())
else:
  define_macros = [('WITH_NUMERIC','1')]

# --- The behavior of distutils changed from 2.2 to 2.3. In 2.3, the object
# --- files are always put in a build/temp directory relative to where the
# --- source file is, rather than relative to the main build directory.
# --- This tells distutils to put the objects in the same directory
# --- as the source files.
if sys.hexversion >= 0x020300f0:
  sys.argv += ['--build-temp','']

if machine == 'darwin': os.environ['ARCHFLAGS'] = '-arch i386'
if machine == 'ppc': os.environ['ARCHFLAGS'] = '-arch ppc'

setup (name = "warpC",
       version = '3.0',
       author = 'David P. Grote',
       author_email = "DPGrote@lbl.gov",
       description = "Combines warp's packages into one",
       platforms = "Unix, Windows (cygwin), Mac OSX",
       ext_modules = [Extension('warpC',
                                ['warpC_Forthon.c',
                                 os.path.join(builddir,'Forthon.c'),
                                 'pmath_rng.c','ranf.c','ranffortran.c'],
                                include_dirs=include_dirs,
                                library_dirs=library_dirs,
                                libraries=libraries,
                                define_macros=define_macros,
                                extra_objects=warpobjects,
                                extra_link_args=['-g']+
                                             fcompiler.extra_link_args,
                                extra_compile_args=fcompiler.extra_compile_args
                               )]

       )
                                #'superlu_dgssv.c'],
                                #include_dirs=[builddir,'SuperLU'],
                                #library_dirs=library_dirs+['SuperLU'],
                                #libraries=libraries+['superlu_3.0','blas'],
