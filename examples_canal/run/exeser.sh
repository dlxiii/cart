#!/bin/bash

# This file is writing to calculate Tokyo Bay run case.

export FVCOMNAME=fvcom_canal

export CASENAME=test

export LOGFILE=$CASENAME'_ser.log'

cp $CASENAME'_'$FVCOMNAME'_run.nml' $CASENAME'_run.nml'

export NUM_CORES=36
export OMP_NUM_THREADS=$NUM_CORES

# ./fvcom --casename=$CASENAME
# mpirun -np $NUM_CORES ./$FVCOMNAME --casename=$CASENAME --logfile=$LOGFILE
# mpirun -np $NUM_CORES ./$FVCOMNAME --casename=$CASENAME
./$FVCOMNAME --casename=$CASENAME --dbg=0 --logfile=$LOGFILE