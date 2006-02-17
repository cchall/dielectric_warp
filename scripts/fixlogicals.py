from warp import *
fixlogicals_version = "$Id: fixlogicals.py,v 1.2 2006/02/17 23:52:36 dave Exp $"

def fixlogicals():
  if cir.lcirout: cir.lcirout = true
  if cir.ltdependent: cir.ltdependent = true
  if cir.lezbeam: cir.lezbeam = true
  if cir.lperveance: cir.lperveance = true
  if cir.lemittance: cir.lemittance = true
  if cir.lallez: cir.lallez = true
  if cir.llinear: cir.llinear = true
  if cir.limage: cir.limage = true
  if cir.lendzero: cir.lendzero = true
  if cir.lsavehist: cir.lsavehist = true
  if env.lenvout: env.lenvout = true
  if env.lefofz: env.lefofz = true
  if env.libeame_z: env.libeame_z = true
  if env.lemitne_z: env.lemitne_z = true
  if f3d.loadquad: f3d.loadquad = true
  if f3d.lchebshv: f3d.lchebshv = true
  if f3d.lzerophiedge: f3d.lzerophiedge = true
  if f3d.lplates: f3d.lplates = true
  if f3d.lcndbndy: f3d.lcndbndy = true
  if f3d.lsrlinr: f3d.lsrlinr = true
  if f3d.lsrminlinr: f3d.lsrminlinr = true
  if f3d.lsrmaxlinr: f3d.lsrmaxlinr = true
  if top.lrelativ: top.lrelativ = true
  if top.lacclzl: top.lacclzl = true
  if top.diposet: top.diposet = true
  if top.laccumulate_zmoments: top.laccumulate_zmoments = true
  if top.lprntpara: top.lprntpara = true
  if top.loneliner: top.loneliner = true
  if top.lpsplots: top.lpsplots = true
  if top.lonedplts: top.lonedplts = true
  if top.ifgap: top.ifgap = true
  if top.allspecl: top.allspecl = true
  if top.laccumulate_rho: top.laccumulate_rho = true
  if top.periinz: top.periinz = true
  if top.stickyz: top.stickyz = true
  if top.stickyxy: top.stickyxy = true
  if top.lgridqnt: top.lgridqnt = true
  if top.lbeamcom: top.lbeamcom = true
  if top.lvinject: top.lvinject = true
  if top.lhlinechg: top.lhlinechg = true
  if top.lhvzofz: top.lhvzofz = true
  if top.lhepsxz: top.lhepsxz = true
  if top.lhepsyz: top.lhepsyz = true
  if top.lhepsnxz: top.lhepsnxz = true
  if top.lhepsnyz: top.lhepsnyz = true
  if top.lhepsgz: top.lhepsgz = true
  if top.lhepshz: top.lhepshz = true
  if top.lhepsngz: top.lhepsngz = true
  if top.lhepsnhz: top.lhepsnhz = true
  if top.lhxbarz: top.lhxbarz = true
  if top.lhybarz: top.lhybarz = true
  if top.lhxybarz: top.lhxybarz = true
  if top.lhxrmsz: top.lhxrmsz = true
  if top.lhyrmsz: top.lhyrmsz = true
  if top.lhxprmsz: top.lhxprmsz = true
  if top.lhyprmsz: top.lhyprmsz = true
  if top.lhxsqbarz: top.lhxsqbarz = true
  if top.lhysqbarz: top.lhysqbarz = true
  if top.lhvxbarz: top.lhvxbarz = true
  if top.lhvybarz: top.lhvybarz = true
  if top.lhvzbarz: top.lhvzbarz = true
  if top.lhxpbarz: top.lhxpbarz = true
  if top.lhypbarz: top.lhypbarz = true
  if top.lhvxrmsz: top.lhvxrmsz = true
  if top.lhvyrmsz: top.lhvyrmsz = true
  if top.lhvzrmsz: top.lhvzrmsz = true
  if top.lhxpsqbarz: top.lhxpsqbarz = true
  if top.lhypsqbarz: top.lhypsqbarz = true
  if top.lhxxpbarz: top.lhxxpbarz = true
  if top.lhyypbarz: top.lhyypbarz = true
  if top.lhxypbarz: top.lhxypbarz = true
  if top.lhyxpbarz: top.lhyxpbarz = true
  if top.lhxpypbarz: top.lhxpypbarz = true
  if top.lhxvzbarz: top.lhxvzbarz = true
  if top.lhyvzbarz: top.lhyvzbarz = true
  if top.lhvxvzbarz: top.lhvxvzbarz = true
  if top.lhvyvzbarz: top.lhvyvzbarz = true
  if w3d.bndfflag: w3d.bndfflag = true
  if w3d.bnprflag: w3d.bnprflag = true
  if w3d.bnezflag: w3d.bnezflag = true
  if w3d.bnjtflag: w3d.bnjtflag = true
  if w3d.l2symtry: w3d.l2symtry = true
  if w3d.l4symtry: w3d.l4symtry = true
  if w3d.lbeforefs: w3d.lbeforefs = true
  if w3d.lafterfs: w3d.lafterfs = true
  if w3d.lgetese3d: w3d.lgetese3d = true
  if w3d.lrhodia3d: w3d.lrhodia3d = true
  if w3d.lpltfld3d: w3d.lpltfld3d = true
  if w3d.cigarld: w3d.cigarld = true
  if w3d.cylinder: w3d.cylinder = true
  try:
    if wrz.fftdiag: wrz.fftdiag = true
    if wrz.vecrho: wrz.vecrho = true
    if wrz.lbeforefs: wrz.lbeforefs = true
    if wrz.lafterfs: wrz.lafterfs = true
    if wrz.cigarld: wrz.cigarld = true
    if wrz.cylinder: wrz.cylinder = true
    if wrz.shell: wrz.shell = true
    if wrz.linbeam: wrz.linbeam = true
  except:
    pass
  if wxy.lvzchang: wxy.lvzchang = true
  if wxy.lthick: wxy.lthick = true
  if wxy.lvp3d: wxy.lvp3d = true
  if wxy.ldiag: wxy.ldiag = true
  if wxy.lexbend: wxy.lexbend = true
