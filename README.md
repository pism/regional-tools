# PISM Drainage basin delineation tool

## Overview

Regional modeling of outlet glaciers using PISM requires isolating a drainage
basin.

Our goal was to develop a tool using surface topography to divide an ice sheet
into smaller basins without using measured ice surface speed (which may not be
available in some cases).

Please see `doc/method.tex` for details (unfinished).

## Requirements

The C code that does the heavy lifting uses an ODE solver from
[GSL](http://www.gnu.org/software/gsl/), although hand-coding one of standard
time-stepping methods would probably result in code that performs about as
well.

This code can be used separately (i.e. in a C-only program).

In addition to using GSL, some functions are parallelized using OpenMP, so a
compiler supporting OpenMP 2.5 or later is needed (GCC 4.2 and later is OK.)

[Cython](http://cython.org/) is used to make the C code mentioned above
available from Python.

The wrapper script `pism_regional.py` uses several Python modules, notably

- [netcdf4-python](http://code.google.com/p/netcdf4-python/) for NetCDF I/O
- [NumPy](http://numpy.scipy.org/) for 2D arrays, etc
- [matplotlib](http://matplotlib.sourceforge.net/) for plotting
- [Tkinter](http://wiki.python.org/moin/TkInter) for the GUI

All these libraries and packages are available via a package manager on
Linux systems.

## Installation

To install for the current user, run

    python setup.py install --user

To build and use from the current directory, run

    python setup.py build_ext --inplace

To disable OpenMP, run

    NO_OPENMP=1 python setup.py ...

If you have GSL in a non-standard location, run

    GSL_PREFIX=/path/to/gsl python setup.py ...

On systems where the default compiler does not support OpenMP you can specify
the compiler to use like this:

    CC=gcc-4.2 CXX=g++-4.2 python setup.py ...

## Usage

- Run `pism_regional.py`. Select a NetCDF file containing variables `x`,
   `y`, `usurf`, and `thk`. 2D arrays have to be stored in the `(y,x)`
   order.
- Select the terminus region using the mouse.
- Click "Compute the drainage basin mask"
- Repeat the terminus region selection and re-generate the mask if necessary
- Save the mask to a file.
- Quit the script by closing all windows.
