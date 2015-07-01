#! /bin/bash

clear
rm g_nitrile_hbond examples/*{xvg,log} *xvg *log
make -j 8
if [ ! -f g_nitrile_hbond ] ; then
    exit
    fi
if [ ! -z $1 ] ; then
    ./g_nitrile_hbond -h
else
cd=$(grep CNC examples/gro.gro | grep CD | awk '{print $3}')
ne=$(grep CNC examples/gro.gro | grep NE | awk '{print $3}')
time ./g_nitrile_hbond -s examples/tpr.tpr -f examples/xtc.xtc -a1 $cd -a2 $ne -oa -op -select "resname SOL and same residue as within .5 of resname CNC and name NE"
    fi

