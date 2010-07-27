em3d
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for the 3-D EM solver of code WARP
# Jean-Luc Vay,   LBNL, (510)486-4934
# David P. Grote, LLNL, (510)423-7194
# Alex Friedman,  LLNL, (510)422-0827

*********** EM3D_APML:
pml              integer /1/
pml_sadjusted    integer /2/
apml_exponential integer /3/
apml_hybrid      integer /4/
apml_ssa         integer /5/
apml_lwa         integer /6/

*********** EM3D_bnd dump:
l_pml_cummer  logical    /.false./
s_max_init       real    /4./
s_max_x          real
s_max_y          real
s_delta          real    /5./
sb_coef          real    /0./
nn               real    /2./
bnd_cond         integer /2/

*********** EM3D_kyee dump:
alphax real /0.58333333333333337/  # 7./12.
betax  real /0.083333333333333329/ # 1./12.
gammax real /0.020833333333333332/ # 1./48.
alphay real /0.58333333333333337/  # 7./12.
betay  real /0.083333333333333329/ # 1./12.
gammay real /0.020833333333333332/ # 1./48.
alphaz real /0.58333333333333337/  # 7./12.
betaz  real /0.083333333333333329/ # 1./12.
gammaz real /0.020833333333333332/ # 1./48.

*********** EM3D_FIELDobjects dump:
l_onegrid                    logical /.true./
l_elaser_out_plane           logical /.false./
l_moving_window              logical /.false./
l_noinputfield               logical /.false./
l_copyfields                 logical /.false./
l_smoothdensity              logical /.false./
ntamp_apr                    integer /4/
rap                          integer /1/
ndelta_t                     integer /1/
nxpatch                      integer /1/
nypatch                      integer /1/
ixpatch                      integer /0/
iypatch                      integer /0/
ntamp_scatter                integer /2/
ntamp_gather                 integer /4/
transition_zone              real /0./ # length of zone for linear transition from coarse to fine force (in coarse cell units)
tmin_moving_main_window      real /0./
otherproc                    integer /10/
otherblock                   integer /11/
push_em3d_e(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_b(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_ef(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_f(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_phi(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_a(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_block(f:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_eef(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical,l_pushpot:logical) subroutine
push_em3d_bf(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical,l_pushpot:logical) subroutine
init_splitfield(sf:EM3D_SPLITYEEFIELDtype, 
                nx:integer,ny:integer,nz:integer, 
                nxguard:integer,nyguard:integer,nzguard:integer, 
                dt:real,dx:real,dy:real,dz:real,
                xmin:real,ymin:real,zmin:real,clight:real,
                lsx:integer,lsy:integer,lsz:integer, 
                nnx:integer, smaxx:real, sdeltax:real, 
                nny:integer, smaxy:real, sdeltay:real, 
                nnz:integer, smaxz:real, sdeltaz:real, 
                l_2dxz:logical, l_2drz:logical) subroutine
depose_jxjyjz_esirkepov_linear_serial(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical)
                           subroutine
depose_jxjyjz_esirkepov_n(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           vx(n):real,vy(n):real,vz(n):real,gaminv(n):real,
                           w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           nox:integer,noy:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_jxjyjz_pxpypz_esirkepov_linear_serial(cj:real,mp:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,m:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical,
                           l_relativ:logical)
                           subroutine
depose_rho_linear_serial(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w(n):real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical)
                           subroutine
depose_rho_n(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           nox:integer,noy:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_rho_n_2dxz(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w:real,q:real,
                           xmin:real,zmin:real,
                           dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical,l_2drz:logical)
                           subroutine
depose_j_n_2dxz(cj:real,
                           n:integer,x(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,
                           dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
getf3d_linear(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real)
                           subroutine
getf3d_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,ymin:real,zmin:real,
         dx:real,dy:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noy:integer,noz:integer,
         exg:real,eyg:real,ezg:real,l4symtry:logical)
                           subroutine
getf2dxz_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,zmin:real,
         dx:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noz:integer,
         exg:real,eyg:real,ezg:real,l4symtry:logical,l_2drz:logical)
                           subroutine
gete3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real)
                           subroutine
getb3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               bxg:real,byg:real,bzg:real)
                           subroutine
geteb3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real,
               bxg:real,byg:real,bzg:real)
                           subroutine
gete3d_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               nox:integer,noy:integer,noz:integer,
               exg:real,eyg:real,ezg:real,
                           l4symtry:logical)
                           subroutine
getb3d_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               nox:integer,noy:integer,noz:integer,
               bxg:real,byg:real,bzg:real,
                           l4symtry:logical)
                           subroutine
gete2dxz_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,zmin:real,
               dx:real,dz:real,
               nx:integer,nz:integer,
               nxguard:integer,nzguard:integer,
               nox:integer,noz:integer,
               exg:real,eyg:real,ezg:real,
                           l4symtry:logical,l_2drz:logical)
                           subroutine
getb2dxz_n_energy_conserving(n:integer,xp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,zmin:real,
               dx:real,dz:real,
               nx:integer,nz:integer,
               nxguard:integer,nzguard:integer,
               nox:integer,noz:integer,
               bxg:real,byg:real,bzg:real,
                           l4symtry:logical)
                           subroutine
yee2node3d(f:EM3D_YEEFIELDtype) subroutine
node2yee3d(f:EM3D_YEEFIELDtype) subroutine
em3d_exchange_e(b:EM3D_BLOCKtype) subroutine
em3d_exchange_b(b:EM3D_BLOCKtype) subroutine
em3d_exchange_f(b:EM3D_BLOCKtype) subroutine
em3d_exchange_j(b:EM3D_BLOCKtype) subroutine
em3d_exchange_rho(b:EM3D_BLOCKtype) subroutine
add_current_slice_3d(f:EM3D_YEEFIELDtype,i:integer) subroutine
add_rho_slice_3d(f:EM3D_YEEFIELDtype,i:integer) subroutine
set_incond(f:EM3D_YEEFIELDtype,n:integer,indx(3,n):integer) subroutine
em3d_applybc_rho(f:EM3D_YEEFIELDtype,xlbnd:integer,xrbnd:integer,
                                     ylbnd:integer,yrbnd:integer,
                                     zlbnd:integer,zrbnd:integer) subroutine
em3d_applybc_j(f:EM3D_YEEFIELDtype,xlbnd:integer,xrbnd:integer,
                                   ylbnd:integer,yrbnd:integer,
                                   zlbnd:integer,zrbnd:integer) subroutine
project_jxjyjz(jfine:real,jcoarse:real,jcoarse_mother:real,
               nxf:integer,nyf:integer,nzf:integer,
               nxc:integer,nyc,nzc:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               rapx:integer,rapy:integer,rapz:integer,
               ixc:integer,iyc:integer,izc:integer,l_2dxz:logical,
               icycle:integer,novercycle:integer) subroutine
project_rho(rhofine:real,rhocoarse:real,rhocoarse_mother:real,
               nxf:integer,nyf:integer,nzf:integer,
               nxc:integer,nyc,nzc:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               rapx:integer,rapy:integer,rapz:integer,
               ixc:integer,iyc:integer,izc:integer,l_2dxz:logical) subroutine
apply_dmask(rho:real,jc:real,dmaskx:real,dmasky:real,dmaskz:real,
            bounds(6):integer,nguarddepos(3):integer,ntrans(3):integer,
            nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,
            l_pushf:logical,l_2dxz:logical) subroutine
addsubstractfields(child:EM3D_BLOCKtype,child_coarse:EM3D_BLOCKtype,
                   parent:EM3D_BLOCKtype,lc(3):integer,ref(3):integer,l_2dxz:logical) subroutine
addsubstractfields_nodal(child:EM3D_BLOCKtype,child_coarse:EM3D_BLOCKtype,
                   parent:EM3D_BLOCKtype,lc(3):integer,ref(3):integer,l_2dxz:logical) subroutine
shift_em3dblock_ncells_x(b:EM3D_BLOCKtype,n:integer) subroutine
shift_em3dblock_ncells_z(b:EM3D_BLOCKtype,n:integer) subroutine
depose_jxjy_esirkepov_linear_serial_2d(j:real,
                           n:integer,x(n):real,y(n):real,
                           xold(n):real,yold(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,ymin:real,dt:real,dx:real,dy:real,
                           nx:integer,ny:integer,l_particles_weight:logical)
                           subroutine
depose_jxjyjz_esirkepov_n_2d(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical,l_2drz:logical)
                           subroutine
depose_jxjyjz_villasenor_n_2d(j:real,
                           n:integer,x(n):real,y(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
setebp(emblock:EM3D_YEEFIELDtype,icycle:integer,novercycle:integer) subroutine

%%%%%%%% EM3D_SPLITYEEFIELDtype:
fieldtype integer /-2/
stencil integer /0/ # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F
nx integer
ny integer
nz integer
nxguard integer /1/
nyguard integer /1/
nzguard integer /1/
ixmin integer /0/ # position of first node of grid interior in the x direction (FORTRAN indexing)
iymin integer /0/ # position of first node of grid interior in the y direction (FORTRAN indexing)
izmin integer /0/ # position of first node of grid interior in the z direction (FORTRAN indexing)
ixmax integer /obj__%nx/ # position of last node of grid interior in the x direction (FORTRAN indexing)
iymax integer /obj__%ny/ # position of last node of grid interior in the y direction (FORTRAN indexing)
izmax integer /obj__%nz/ # position of last node of grid interior in the z direction (FORTRAN indexing)
ixming integer /-obj__%nxguard/ # position of first node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iyming integer /-obj__%nyguard/ # position of first node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izming integer /-obj__%nzguard/ # position of first node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
ixmaxg integer /obj__%ixmax+obj__%nxguard/ # position of last node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iymaxg integer /obj__%iymax+obj__%nyguard/ # position of last node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izmaxg integer /obj__%izmax+obj__%nzguard/ # position of last node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
jxmin integer /0/ # position of first node of grid interior in the x direction (Python indexing)
jymin integer /0/ # position of first node of grid interior in the y direction (Python indexing)
jzmin integer /0/ # position of first node of grid interior in the z direction (Python indexing)
jxmax integer /0/ # position of last node of grid interior in the x direction (Python indexing)
jymax integer /0/ # position of last node of grid interior in the y direction (Python indexing)
jzmax integer /0/ # position of last node of grid interior in the z direction (Python indexing)
jxming integer /0/ # position of first node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jyming integer /0/ # position of first node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzming integer /0/ # position of first node of entire grid (interior+guard nodes) in the z direction (Python indexing)
jxmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jymaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the z direction (Python indexing)
nxpo integer /0/
nypo integer /0/
nzpo integer /0/
dx real
dy real
dz real
dxi real
dyi real
dzi real
dt real
xmin real
xmax real
ymin real
ymax real
zmin real
zmax real
clight real
lsx integer
nnx integer
smaxx real
sdeltax real
lsy integer
nny integer
smaxy real
sdeltay real
lsz integer
nnz integer
smaxz real
sdeltaz real
l_2dxz logical /.False./
l_2drz logical /.False./
exx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ax(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
ay(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
az(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
phi(1:3,-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
afx(-nxguard:nx+nxguard) _real
bpfx(-nxguard:nx+nxguard) _real
bmfx(-nxguard:nx+nxguard) _real
agx(-nxguard:nx+nxguard) _real
bpgx(-nxguard:nx+nxguard) _real
bmgx(-nxguard:nx+nxguard) _real
afy(-nyguard:ny+nyguard) _real
bpfy(-nyguard:ny+nyguard) _real
bmfy(-nyguard:ny+nyguard) _real
agy(-nyguard:ny+nyguard) _real
bpgy(-nyguard:ny+nyguard) _real
bmgy(-nyguard:ny+nyguard) _real
afz(-nzguard:nz+nzguard) _real
bpfz(-nzguard:nz+nzguard) _real
bmfz(-nzguard:nz+nzguard) _real
agz(-nzguard:nz+nzguard) _real
bpgz(-nzguard:nz+nzguard) _real
bmgz(-nzguard:nz+nzguard) _real

%%%%%%%% EM3D_YEEFIELDtype:
fieldtype integer /-1/
stencil integer /0/ # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F
nx integer /0/ # nb of mesh cells of grid interior in the x direction
ny integer /0/ # nb of mesh cells of grid interior in the y direction
nz integer /0/ # nb of mesh cells of grid interior in the z direction
nxguard integer /1/ # nb of guard cells in the x direction
nyguard integer /1/ # nb of guard cells in the y direction
nzguard integer /1/ # nb of guard cells in the z direction
ixmin integer /0/ # position of first node of grid interior in the x direction (FORTRAN indexing)
iymin integer /0/ # position of first node of grid interior in the y direction (FORTRAN indexing)
izmin integer /0/ # position of first node of grid interior in the z direction (FORTRAN indexing)
ixmax integer /obj__%nx/ # position of last node of grid interior in the x direction (FORTRAN indexing)
iymax integer /obj__%ny/ # position of last node of grid interior in the y direction (FORTRAN indexing)
izmax integer /obj__%nz/ # position of last node of grid interior in the z direction (FORTRAN indexing)
ixming integer /-obj__%nxguard/ # position of first node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iyming integer /-obj__%nyguard/ # position of first node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izming integer /-obj__%nzguard/ # position of first node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
ixmaxg integer /obj__%ixmax+obj__%nxguard/ # position of last node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iymaxg integer /obj__%iymax+obj__%nyguard/ # position of last node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izmaxg integer /obj__%izmax+obj__%nzguard/ # position of last node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
jxmin integer /0/ # position of first node of grid interior in the x direction (Python indexing)
jymin integer /0/ # position of first node of grid interior in the y direction (Python indexing)
jzmin integer /0/ # position of first node of grid interior in the z direction (Python indexing)
jxmax integer /0/ # position of last node of grid interior in the x direction (Python indexing)
jymax integer /0/ # position of last node of grid interior in the y direction (Python indexing)
jzmax integer /0/ # position of last node of grid interior in the z direction (Python indexing)
jxming integer /0/ # position of first node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jyming integer /0/ # position of first node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzming integer /0/ # position of first node of entire grid (interior+guard nodes) in the z direction (Python indexing)
jxmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jymaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the z direction (Python indexing)
nxp integer /0/
nyp integer /0/
nzp integer /0/
nxext integer /0/
nyext integer /0/
nzext integer /0/
nxdamp integer /0/
nydamp integer /0/
nzdamp integer /0/
nxpnext integer /0/
nypnext integer /0/
nzpnext integer /0/
nxf integer /0/
nyf integer /0/
nzf integer /0/
nxpo integer /0/
nypo integer /0/
nzpo integer /0/
nxmp integer /0/
nymp integer /0/
nzmp integer /0/
ntimes integer /1/
nconds integer /0/
nxcond integer /0/
nycond integer /0/
nzcond integer /0/
xmin real
ymin real
zmin real
xmax real
ymax real
zmax real
dx real
dy real
dz real
dxi real
dyi real
dzi real
clight real
mu0    real
theta_damp real /0./
l_2dxz logical /.False./
l_2drz logical /.False./
sigmae real /0./ # coefficient for extended solver
sigmab real /0./ # coefficient for extended solver
Ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
By(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Exp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Eyp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Ezp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Bxp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Byp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Bzp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Expnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Eypnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Ezpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bxpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bypnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bzpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
F(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard) _real
Rho(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard) _real
Rhoold(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard) _real
Rhoarray(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard,ntimes) _real
J(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) _real
Jarray(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3,ntimes) _real
incond(-nxguard:nxcond+nxguard,-nyguard:nycond+nyguard,-nzguard:nzcond+nzguard) _logical
Ax(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Ay(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Az(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Phi(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Mp(-nxguard:nxmp+nxguard,-nyguard:nymp+nyguard,-nzguard:nzmp+nzguard,3) _real
Exold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eyold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Exbar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eybar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezbar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Excp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eycp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezcp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
DEXY(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DEXZ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DEYX(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DEYZ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DEZX(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DEZY(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBXY(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBXZ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBYX(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBYZ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBZX(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
DBZY(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BXYCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BYXCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BXZCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BZXCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BYZCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BZYCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EXYCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EYXCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EXZCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EZXCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EYZCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EZYCJ(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BXYCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BYXCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BXZCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BZXCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BYZCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
BZYCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EXYCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EYXCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EXZCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EZXCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EYZCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
EZYCJT(-nxguard:nxext+nxguard,-nyguard:nyext+nyguard,-nzguard:nzext+nzguard) _real
E_inx_pos integer /-1/
E_inx_angle real  /0./
E_inx(-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
E_iny_pos integer /-1/
E_iny_angle real  /0./
E_iny(-nxguard:nx+nxguard,-nzguard:nz+nzguard) _real
E_inz_pos real /0./
E_inz_angle real  /0./
Ex_inz(-nxguard:nx+nxguard,-nyguard:ny+nyguard) _real
Ey_inz(-nxguard:nx+nxguard,-nyguard:ny+nyguard) _real
dmaskx(-nxguard:nx+nxguard) _real
dmasky(-nyguard:ny+nyguard) _real
dmaskz(-nzguard:nz+nzguard) _real

%%%%%%%% EM3D_FIELDtype:
fieldtype integer /0/
yf _EM3D_YEEFIELDtype
syf _EM3D_SPLITYEEFIELDtype
proc integer /0/
xl _EM3D_FIELDtype
xr _EM3D_FIELDtype
yl _EM3D_FIELDtype
yr _EM3D_FIELDtype
zl _EM3D_FIELDtype
zr _EM3D_FIELDtype

%%%%%%%% EM3D_BLOCKtype:
nx integer
ny integer
nz integer
nxguard integer 
nyguard integer 
nzguard integer 
nbndx integer
nbndy integer
nbndz integer
xmin real
ymin real
zmin real
xmax real
ymax real
zmax real
dx real
dy real
dz real
dxi real
dyi real
dzi real
xlbnd                      integer /0/
xrbnd                      integer /0/
ylbnd                      integer /0/
yrbnd                      integer /0/
zlbnd                      integer /0/
zrbnd                      integer /0/
core _EM3D_FIELDtype
sidexl _EM3D_FIELDtype 
sidexr _EM3D_FIELDtype
sideyl _EM3D_FIELDtype
sideyr _EM3D_FIELDtype
sidezl _EM3D_FIELDtype
sidezr _EM3D_FIELDtype
edgexlyl _EM3D_FIELDtype
edgexryl _EM3D_FIELDtype
edgexlyr _EM3D_FIELDtype
edgexryr _EM3D_FIELDtype
edgexlzl _EM3D_FIELDtype
edgexlzr _EM3D_FIELDtype
edgexrzl _EM3D_FIELDtype
edgexrzr _EM3D_FIELDtype
edgeylzl _EM3D_FIELDtype
edgeyrzl _EM3D_FIELDtype
edgeylzr _EM3D_FIELDtype
edgeyrzr _EM3D_FIELDtype
cornerxlylzl  _EM3D_FIELDtype
cornerxrylzl  _EM3D_FIELDtype
cornerxlyrzl  _EM3D_FIELDtype
cornerxryrzl  _EM3D_FIELDtype
cornerxlylzr  _EM3D_FIELDtype
cornerxrylzr  _EM3D_FIELDtype
cornerxlyrzr  _EM3D_FIELDtype
cornerxryrzr  _EM3D_FIELDtype
