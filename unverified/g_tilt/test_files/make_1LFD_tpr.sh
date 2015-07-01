#! /bin/bash

grep " A \| B \| GNP B" 1LFD.pdb > 1LFD-2.pdb
if [ ! -f 1LFD.gro ] ; then
    echo 1 | pdb2gmx -f 1LFD-2.pdb -o 1LFD.gro -p 1LFD -i 1LFD -water tip3p -merge all
    fi
if [ ! -f 1LFD.tpr ] && [ -f 1LFD.gro ] ; then
    grompp -f junk.mdp -c 1LFD.gro -p 1LFD -o 1LFD.tpr
    rm mdout.mdp
    fi     
if [ ! -f rap_aligned.tpr ] ; then
    echo 1 | pdb2gmx -f rap_aligned.gro -o junk.gro -p junk -i junk -water tip3p -merge all
    grompp -f junk.mdp -c rap_aligned.gro -p junk -o rap_aligned.tpr
    rm junk* mdout.mdp
    fi
rm *#
