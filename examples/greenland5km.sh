#!/bin/bash

# Copyright (C) 2009-2012 Ed Bueler and Andy Aschwanden and Constantine Khroulev

# downloads SeaRISE "Present Day Greenland" master dataset NetCDF file,
# adjusts metadata, saves under new name

# depends on wget and NCO (ncrename, ncap, ncwa, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland

DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc
DATASIZE=11Mb

echo "fetching $DATASIZE master file ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo

WORKING=gris_5km.nc
echo "creating simplified file $WORKING with vars x,y,usurf,thk from master ..."

# copies over preserving history and global attrs:
ncks -O -v topg,thk,usrf,x1,y1 $DATANAME $WORKING

# rename for use with pism_regional.py
ncrename -d x1,x -d y1,y -v x1,x -v y1,y gris_5km.nc
ncrename -v usrf,usurf $WORKING

# remove time dimension
ncwa -O -a time $WORKING $WORKING

echo "done with cleaning file $WORKING"
echo "now do 'python ../pism_regional.py' and open $WORKING"
echo

