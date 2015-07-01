#! /bin/bash

cd ..
rm -v gmx2amoeba_cp
make
cd -

tndx=ndx.ndx
#tndx=t.ndx
#echo "[ test ]" > t.ndx
#for i in `seq 372 373`; do echo $i >> t.ndx; done
#for i in `seq 1 10`; do echo $i >> t.ndx ; done
#for i in `seq 340 380`; do echo $i >> t.ndx ; done
rm -f *.xvg
if [ -z $1 ] ; then
time ../gmx2amoeba_cp -s tpr.tpr -f xtc.xtc -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm -a1 372 -a2 373 -exclude "368 369 370 371" -e 15
else
for i in `seq 0 5`; do
find . -name "*.xyz" -print | xargs rm
time ../gmx2amoeba_cp -s tpr.tpr -f gro.gro -p top.top -n $tndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm  &> j
rm j
done
fi
