#!/bin/bash

# A script to make sure all processors having the same starting time

if [ x$1 = "x-h" ]; then
    echo "Remove the latested time dir if some other CPUs do not have this timedir "
    echo "If no argument is given, the current dir is cleaned."
    echo "Usage: uniformizeFoamTime [case]"
fi


uniformResultsInDir() 
{
    if [ $# -lt 1 ]; then
        local operateDir="."
    else
        operateDir=$1
    fi
    
    local dirs=$(ls -l $operateDir | grep ^d | awk '{print $9}' | grep -E '^[0-9]+(\.[0-9]+)?$'  | grep -v -e '^0$' | sort -n -r)
    # Note: still dangerous, what if 0 is followed spaces? Need to refine.

    local dirArray=( $dirs )
    timeDir="$operateDir/${dirArray[0]}"
        
    echo "Removing $timeDir "
    /bin/rm -fr $timeDir
 
}

scanDir() 
{
    if [ $# -lt 1 ]; then
        local operateDir="."
    else
        operateDir=$1
    fi
    
    local dirs=$(ls -l $operateDir | grep ^d | awk '{print $9}' | grep -E '^[0-9]+(\.[0-9]+)?$'  | grep -v -e '^0$' | sort -n -r)
    local dirArray=( $dirs )

    echo ${dirArray[0]}
}



if [ $# -lt 1 ]; then
    allCases="."
else
    allCases="$@"
fi

for case in $allCases; do
    if [ ! -d $case/processor1 ] || [ ! -d $case/system ]; then
        echo "Warning: '$case' is not a parallel OpenFOAM case (or not a case at all)."
        continue
    else
        echo uniformize case: $case
        #scan dir once to determine uniformness
        tmpfile=tmp.$case
        if [ -f $tmpfile ]; then
            /bin/rm -fr $tmpfile
            touch $tmpfile
        fi

        for procDir in $case/processor*; do
            maxdir=$( scanDir $procDir )
            echo $maxdir >> $tmpfile
        done

        nt=$(cat $tmpfile | sort -n | uniq | wc -l)
        maxval=$(cat $tmpfile | sort -n -r | head -n 1)

        /bin/rm $tmpfile

        if [ $nt -gt 1 ]; then
            for procDir in $case/processor*; do
                if [ -d $procDir/$maxval ]; then
                    uniformResultsInDir $procDir
                fi
            done
        else
            echo "Case $case is uniform."
        fi

        echo "Finished case: $case"; echo ""
            
    fi
done
