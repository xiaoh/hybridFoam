#!/bin/bash

module load gcc/4.3.3
module load open_mpi/1.4.1
module load qt/4.3.5
module load cmake
.  $INSTALLATION_DIR/etc/bashrc


# --------------------------------------------------------------------------
# Only edit this section: 
EXE="channelFoam"

#wtime="-W 7:59" # Requested wall time (serious computing)
#NR=1  # Number of chained runs

wtime=" "      # For test purpose only
NR=0   # Number of chained runs


# ----------------------------------------------------------------------------
# Usually there is no need to edit the section below
# ----------------------------------------------------------------------------
# Number of processors read from system/decomposeParDict
NP=`sed -n '/numberOfSubdomains/s/[^0-9]//gp' system/decomposeParDict `   

VARS="-x WM_PROJECT_DIR -x MPI_BUFFER_SIZE"

echo "Number of processors: $NP"

# Name of the chain set to current dir name
chainName=`basename $PWD`

"Request time: $wtime on $NP processors with $NR chained jobs."

bsub $wtime -n $NP   -J $chainName mpirun -np $NP $VARS $EXE -parallel

for i  in `seq $NR`; do
    echo "Increment # " $i
    bsub $wtime -n $NP -J $chainName -w "ended($chainName)" "~/bin/setNewEndTime.sh; mpirun -np $NP $VARS $EXE -parallel"
done
