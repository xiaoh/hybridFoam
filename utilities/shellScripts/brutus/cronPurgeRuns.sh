#!/bin/bash
# This script calls "purgeFoam" to clean the data on Brutus
# This script is called by a cron job from etzel

inputDir=$1

if [ -z "$inputDir" ]; then
    runsDir='/cluster/work/mavt/xiaoh/OF-runs' # Runs on Brutus
else
    echo "input: $inputDir"
    runsDir=$inputDir
fi

cd $runsDir

if [ -f $runsDir"/PurgeLock" ]; then
    echo "Another instance of purgeFoam is running. Exit at:" >> ~/.purgeFoamLog 2>&1
    date >> ~/.purgeFoamLog 2>&1; echo "" >> ~/.purgeFoamLog 2>&1
    echo "-------------------------------------------------------" >> ~/.purgeFoamLog 2>&1
    echo "" >> ~/.purgeFoamLog 2>&1
    exit 1
fi

touch $runsDir"/PurgeLock"

echo "*** Starting to purge foam OF runs at:" >> ~/.purgeFoamLog 2>&1
date >> ~/.purgeFoamLog 2>&1; echo "" >> ~/.purgeFoamLog 2>&1

# Actually purging of directory
~/bin/purgeFoam.sh  $runsDir/*  >> ~/.purgeFoamLog 2>&1

echo "@@@ Finished purge OF runs directory at:"  >> ~/.purgeFoamLog 2>&1
date >> ~/.purgeFoamLog 2>&1
echo "-------------------------------------------------------" >> ~/.purgeFoamLog 2>&1
echo "" >> ~/.purgeFoamLog 2>&1
/bin/rm -f $runsDir"/PurgeLock"