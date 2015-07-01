#! /bin/bash

cd ..
rm -v gmx2xyz
make
cd -

#tndx=ndx.ndx
tndx=nosol.ndx
#echo "[ test ]" > t.ndx
#for i in `seq 1 150`; do echo $i >> t.ndx; done
find . -name "*.xyz" -print | xargs rm
if [ -z $1 ] ; then
time ../gmx2xyz -s tpr.tpr -f xtc.xtc -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm 
rm -f conf*.000E ;
for i in $( ls *.xyz ) ; do
    ~ritchie/software/tinker-6.3.3-chgpen/bin/field_parts-2 $i -k ~ritchie/tmp/amoeba_tests/fast.key 372 373 368 371
done
else
for i in `seq 0 5`; do
find . -name "*.xyz" -print | xargs rm
time ../gmx2xyz -s tpr.tpr -f gro.gro -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm  &> j
rm j
done
fi
