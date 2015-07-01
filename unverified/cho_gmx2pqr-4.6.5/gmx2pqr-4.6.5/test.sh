#1 /bin/bash

clear
rm cho_gmx2pqr examples/state0.pqr examples/field_projection.xvg
make -j 8
if [ ! -f cho_gmx2pqr ] ; then
    exit
    fi
if [ ! -z $1 ] ; then
    ./cho_gmx2pqr -h
else
    ./cho_gmx2pqr -s examples/tpr.tpr -f examples/gro.gro -o examples/state.pqr -of examples/field_projection.xvg -n examples/NoSol.ndx -a1 351 -a2 352 -exclude "350 347 348 349" -site "347 350"
#    cat examples/field_projection.xvg
#    grep "CNC\|DUM" examples/state0.pqr
    fi
