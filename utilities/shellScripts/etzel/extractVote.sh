
info='infoFromWeb.txt'
list='processedList.txt'
final='final.txt'

python ~/bin/fetchVotes.py > $info

grep "6:00 PM" $info | grep "$( date +'%B %-d'  -d "this Thursday" )," > $list

echo "The Doodle voting results as of `date`:" > $final
echo "-------------------------------------" >> $final
sed -e 's/[a-z ]*=\"//g' -e 's/[0-9]*//g' -e 's/[,\":]//g' -e 's/-/NO/' -e 's/OK/YES/' -e 's/ifneedbe/TBD/' $list | awk '{printf "%-10s : %-5s\n", toupper($1) ,$NF}' | sort >> $final
echo "-------------------------------------" >> $final

nyes=$(grep YES $final | wc -l)
nno=$(grep NO $final | wc -l)
ntbd=$(grep TBD $final | wc -l)

echo "Summary    : YES=$nyes, NO=$nno, TBD=$ntbd" >> $final
if [ $nyes -gt 4 ]; then
    echo -e "\nNOTE: More than 4 people playing. Rotation withdraw rule applies this week.\n" >> $final
    else
    echo -e "\nNOTE:  Only 4 or fewer are playing. Withdraw rule is deactivated this week.\n" >> $final
fi

cat $final >> $1

rm $list $info $final
