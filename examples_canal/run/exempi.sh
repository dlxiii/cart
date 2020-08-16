#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-s"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=48:00:00"
#PJM -j
#PJM -X

export FVCOMNAME=fvcom_canal

export CASENAME=test

export LOGFILE=$CASENAME'_mpi.log'

cp $CASENAME'_'$FVCOMNAME'_run.nml' $CASENAME'_run.nml'

module load intel/2018

NUM_NODES=${PJM_VNODES}
NUM_CORES=36
NUM_PROCS=36

export I_MPI_HYDRA_BOOTSTARP=rsh
export I_MPI_HYDRA_BOOTSTART_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./$FVCOMNAME --casename=$CASENAME --logfile=$LOGFILE
