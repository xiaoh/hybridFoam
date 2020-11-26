#!/bin/bash

# send reminder, and calculate withdraw week duty
# Fetch the doodle vote results, email to everybody

players=( Haixin Heng Liang Xuan Yousef )

tmpFile="tennis.tmp"
wn=`date +%U`  # compute week number
drn=$(( $wn - 24 )) # use offset to choose who to start

nn=$(( $drn + 1 ))
pn=$(( $drn % 5 ))
nwn=$(( $nn % 5 ))

withdraw=${players[$pn]}
next=${players[$nwn]}

echo Dear friends:  > $tmpFile
echo " "   >> $tmpFile

~/bin/extractVote.sh $tmpFile

echo If your choice listed above is not correct, please update it on Doodle at:   >> $tmpFile
echo http://www.doodle.com/wguhq4aq49utp439  >> $tmpFile
echo " "  >> $tmpFile

echo "---------Week $wn of 2010-------------"   >> $tmpFile
echo Withdraw week for: $withdraw  >> $tmpFile
echo "NEXT WEEK withdraw: $next"  >> $tmpFile
echo "-------------------------------------"   >> $tmpFile
echo " "  >> $tmpFile
echo  -e "To other players: \nIf you are not playing this week, please let $withdraw know, so that he/she can play!"  >> $tmpFile
echo " "  >> $tmpFile
echo "This is an auto-generated email, sent to you every Tuesday and Thursday @ 11:30 AM."   >> $tmpFile
echo "Report any errors to xiao.eth@gmail.com.  Have fun on court!" >> $tmpFile

mail -s "Tennis Irchel for Week $wn"  -c xiao.princeton@gmail.com  ray.enough@gmail.com  xuanshao1980@gmail.com  Yousef.Najafi@usz.ch tian.liang@mol.biol.ethz.ch  < $tmpFile
rm $tmpFile

echo "Finished!"