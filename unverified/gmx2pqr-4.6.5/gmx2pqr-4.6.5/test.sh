#1 /bin/bash

rm gmx2pqr examples/state0.pqr examples/field_projection.xvg
make -j 8
if [ ! -f gmx2pqr ] ; then
    exit
    fi
./gmx2pqr -s examples/tpr.tpr -f examples/gro.gro -o examples/state.pqr -of examples/field_projection.xvg -n examples/NoSol.ndx -a1 214 -a2 215
head examples/state0.pqr
