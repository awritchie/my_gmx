#! /bin/bash

clear
for i in `seq 100`; do echo ; done
make clean
make
if [ -f g_CNC-HOH ] ; then
    clear
    time ./g_CNC-HOH -f test/production-1.xtc -s test/production-1.tpr -o test/hbonds.xvg -op test/persist.xvg -oa test/residues.log
    rm test/*#
    fi
