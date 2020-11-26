#!/bin/bash
# For Foam cases: Increment the endTime by s specified number

dict="system/controlDict";

incTime=$( sed -n '/[Ii]ncrement/s/[^0-9]//gp' $dict )
endTime=$( sed -n '/endTime*[0-9]*/s/[^0-9]//gp' $dict )

echo End time increased by $incTime from $endTime

newEndTime=`expr $incTime + $endTime`

echo "New endTime: " $newEndTime

sed -e '/endTime*[0-9]*/s/\([0-9]\+\)/'"$newEndTime"'/'  $dict > ${dict}-a

str=$(sed -n '/startFrom/p' system/controlDict | sed -n '/latestTime/p')

if [ "x" = "x$str" ]; then
    echo You probably forget to set \"startFrom\" to \"latestTime\"! I will do it for you.
    sed '/startFrom/c startFrom    latestTime;' ${dict}-a > ${dict}-b
    mv ${dict}-b  ${dict}-a
else
    sed -n '/startFrom/p' ${dict}-a
    sed -n '/stopAt/p' ${dict}-a
fi

mv ${dict}-a $dict



