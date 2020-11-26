#!/bin/bash

# A script to clean a case (serial or parallel)

cleanInterval=120

if [ x$1 = "x-h" ]; then
    echo "Clean the given case(s), but keep the latest 2 time directories. "
    echo "If no argument is given, the current dir is cleaned."
    echo "Usage: purgeFoam [case]"
fi

nKeep=2

if [ $nKeep -lt 1 ]; then
    echo At least keep the latest time dir.
    nKeep=1
fi

purgeResultsInDir() 
{
    if [ $# -lt 1 ]; then
        operateDir="."
    else
        operateDir=$1
    fi
    
#   dirs=$(ls -l $operateDir | grep ^d | awk '{print $9}' | grep -v -e '[-#a-zA-Z]' | grep -v -e '^0$' | sort -n -r)
    dirs=$(ls -l $operateDir | grep ^d | awk '{print $9}' | grep -E '^[0-9]+(\.[0-9]+)?$'  | grep -v -e '^0$' | sort -n -r)
    # Note: still dangerous, what if 0 is followed spaces? Need to refine.

    dirArray=( $dirs )

    N=${#dirArray[@]}
    if [ $N -le $nKeep ]; then
        # echo Less than $nKeep dirs. Do nothing.
        echo purged not long ago. pass: $N, $nKeep, $dirs
        return 0
    else
        echo  Deleting these time dirs: 
    fi
    
    M=`expr $N - 1`

    for i in `seq $nKeep $M`; do
        echo -n "processing dir: ${dirArray[i]} "
        timeDir="$operateDir/${dirArray[i]}"

        gunzip $timeDir/impactRestart.gz
        restartNo=$(awk '/restartNumber/{ printf "%03d\n", $2 }' $timeDir/impactRestart)
        # remove hdf5 file corresponding to this time
        if [ -f $operateDir/../pre_restart${restartNo}.h5 ]; then
            echo "removing hdf5 files: " $operateDir/../*restart${restartNo}*
            rm -f $operateDir/../*restart${restartNo}*
        fi
        
        rm -fr $timeDir
    done
    echo " "
}

if [ $# -lt 1 ]; then
    allCases="."
else
    allCases="$@"
fi

for case in $allCases; do
    # echo "Cleaning directory '$case'"
    if [ ! -d $case/system ]; then
        echo "Warning: '$case' is not an OpenFOAM case."
        continue
    else
        # Check if case modified in since last purge
        Ncf=`find $case -mmin -$cleanInterval | wc -l`
        echo "Ncf = " $Ncf
        # At leate 3 files should have been changed if a case is active
        if [ $Ncf -lt 0 ]; then  
            echo "Ignore old case: $case"; echo ""
            continue
        else            
            echo Purge case: $case
            purgeResultsInDir $case yes
            if [ -d $case/processor0 ]; then
                echo `basename $case` is a parallel case. Clean processor dir as well.
                for procDir in $case/processor*; do
                    # echo "Purging dir: $procDir"
                    purgeResultsInDir $procDir
                done
            fi
            echo "Finishing cleaning case: $case"; echo ""
        fi

    fi
done

echo "done!"