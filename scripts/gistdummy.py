# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
gistdummy_version = "$Id: gistdummy.py,v 1.6 2001/10/23 17:08:37 dave Exp $"
from Numeric import *
import sys, os	# To be sure expand_path has posixpath and we have sys.path
#from gistC import *
#import gistC
from help import help
from shapetest import *
from arrayfns import *
from types import *

def hcp_file(*x,**xx):
  pass
def hcpon(*x,**xx):
  pass
def hcp(*x,**xx):
  pass
def window(*x,**xx):
  pass
def plg(*x,**xx):
  pass
def fma():
  pass
def plt(*x,**xx):
  pass
def plsys(i=0):
  pass
def winkill(*N):
  pass
def pltitle(title):
  pass
def limits(*x,**xx):
  pass
def ylimits(ymin='u',ymax='u'):
  pass
def moush(*arg):
  pass
def eps(name):
  import os
  pass
def xytitles(xtitle = "", ytitle = "", delta = (0.,0.)):
  pass
def _spanz(lb,ub,n):
  pass
def plmk(y,x=None,marker=None,width=None,color=None,msize=None):
  pass
def spann (zmin, zmax, n = 8, fudge = 0, force = 0) :
  pass
def plfc (z, y, x, ireg, contours = 20, colors = None, region = 0,
          triangle = None, scale = "lin") :
  pass
def plc(*x,**xx):
  pass
def redraw():
  pass
def current_window():
  pass
def pldefault(*x,**xx):
  pass
def palette(*x,**xx):
  pass
def pli(*x,**xx):
  pass
def pldj(*x,**xx):
  pass
def plfp(*x,**xx):
  pass
