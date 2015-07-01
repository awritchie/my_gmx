#! /bin/bash

cd ..
rm -v gmx2xyz
make
cd -

../gmx2xyz -s tpr.tpr -f gro.gro -p top.top -n nosol.ndx -a ~ritchie/software/tinker/params/CNC_amoebabio09.prm