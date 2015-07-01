#! /bin/bash

clear
for i in `seq 100`; do echo ; done
make clean
make
if [ -f g_insert_dummy_atom ] ; then
    clear
    if [ -z $1 ] ; then
        ./g_insert_dummy_atom -f test_files/rap_aligned.gro -s test_files/rap_aligned.tpr -o test_files/dum.gro -a1 351 -a2 352
        rm test_files/*#
    else
        ./g_insert_dummy_atom -f test_files/rap_aligned.xtc -s test_files/rap_aligned.tpr -o test_files/dum.xtc -a1 351 -a2 352
        rm test_files/*#
        fi
    fi

