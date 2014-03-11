def warpfortrandoc():
    print """
Definges warpfortran() which prints description of most useful fortran
routines accesible from python
"""


def warpfortran():
    print """
derivqty() Calculates global derived qtys.
resetlat() Resizes lattice arrays to their true lengths
setlatt() Sets lattice pointers for the current beam location
setlattzt() Sets lattice pointers for the specified location
getelemid() Find index for element at specified location
species() Sets species related arrays.

getpsgrd() Bins particles on a 2-D, optionally slanted, grid
emitthresh() Calculates emittance thresholding
emitellipse() Calculates emittance versus percent of current in nested ellipses
stepid() Set top line at the bottom of plots
savehist() Saves history data to the arrays
prntpara() Prints global parameters
drawlatt() Draws the lattice (Basis only)
tolabfrm() Converts data from lattice frame to lab frame (including bends)
chckpart() Makes sure that there is space in particles arrays
load2d() Loads particles based on a 2-D distribution
r2rev() Returns a bit reversed number
rnrev() Returns a digit reversed number
rnorm() Returns a Gaussian random number
rnormdig() Returns a Gaussian random number using digit reversed random number
rm() Returns Pseudo-Gaussian random numbers (from 6 uniform random numbers)
rma() Returns list of numbers generated by rm
dvdfit() For singular-value decomposition
dvdcmp() For singular-value decomposition
dvbksb() For singular-value decomposition
writarry() Write an array to file
wtime() Returns current absolute CPU time
wtimeon()  Turns timer on
wtimeoff() Returns time since last call to wtimeon or wtimeoff
divxy() Calculates RMS vx and vy versus x and y
multpole() Calculate the multipole moments of the potential
"""
