#! /bin/bash

for f in $(ls *.000E) ; do
    dip=$(head -n 2 $f | tail -n 1 | sed "s/\;//g")
    pt=$(tail -n 1 $f | awk '{print $3,$5}')
    echo $dip $pt
done