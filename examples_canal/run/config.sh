#!/bin/bash

# This file is writing to calculate Tokyo Bay run case.

export FVCOMNAME=fvcom_canal
export RUNDIR=$(pwd)
export FVCOMDIR=$HOME/github/fvcom4.0/FVCOM_source

# step 1

cp $FVCOMNAME.inc $FVCOMDIR/make.inc
cd $FVCOMDIR
make clean
make

# step 2

cd $RUNDIR
cp $FVCOMDIR/fvcom $RUNDIR/$FVCOMNAME
chmod 754 $FVCOMNAME