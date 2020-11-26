#!/bin/bash

new=$1

cp -r Durbin $new
rm $new/*.dep
cd $new
rename Durbin $new *.[CH]
mv fWallI.H fWallI$new.H
sed -i "s/Durbin/$new/g" $new.[CH]
sed -i "s/fWallI/fWallI$new/g" $new.C
echo "$new/$new.C" | cat - >> ../Make/addfile
cd ..

