# Installing Warp on Debian

The instructions below are valid for Linux distributions that derive
from Debian (including Ubuntu for instance).

### Preparing your Linux environnement

Please make sure that you have `gfortran`, `git` and `make` installed.
On Ubuntu/Debian, this can be done as follows:
```
sudo apt-get update && sudo apt-get install gfortran make git wget libx11-dev
```

## Installing Anaconda

If you do not have Anaconda installed, install it by typing:
```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b
```
Open a new terminal and make sure that the default `python` is the Anaconda version. (Running `python` should print a message which starts with a line of the form `Python 2.7.10 |Anaconda 2.3.0`)

# Installing required packages

Run `conda install numpy scipy h5py mpi4py`

Then make sure that the executable that you call when typing `mpif90` and `mpirun` are those provided by Anaconda. This can be checked by typing `which mpif90` and `which mpirun`.

## Installing Forthon

Run `pip install Forthon`

## Installing pygist

Use the following commands 
```
git clone https://bitbucket.org/dpgrote/pygist.git
cd pygist
ptyhon setup.py config
python setup.py install
```

## Installing Warp

Run ```git clone https://bitbucket.org/berkeleylab/warp.git```

then `cd` into the repository `warp/pywarp90` and create two files:

- A file named `Makefile.local.pympi` which contains the following text:

```FCOMP= -F gfortran --fcompex mpif90```

- A file named `setup.local.py` which contains the following text:

```python
if parallel:
	library_dirs += ['~/anaconda2/lib/']
	libraries = fcompiler.libs + ['mpichf90', 'mpich', 'opa', 'mpl']
```

Then install Warp by running:
```
make install
make pinstall
```

## Running a simulation

In order to run a simulation, create a new directory,
copy your Warp input script to this directory, and rename this script
to `warp_script.py`. (The folder `scripts/examples/` of the
[Warp repository](https://bitbucket.org/berkeleylab/warp/src) contains
several examples of input scripts.)

Then launch the simulation by typing either (for serial simulations)
```
python warp_script.py
```
or (for e.g. a parallel simulation with a 2x3x2 domain decomposition in 3D)
```
mpirun -np 12 python warp_script.py -p 2 3 2
```