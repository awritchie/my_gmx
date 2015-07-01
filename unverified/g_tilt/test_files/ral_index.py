#! /anaconda/bin/python

import sys

LFD="1LFD.gro"
LFDfile = open(LFD)
LFDlines = LFDfile.readlines()
LFDfile.close()

ndxtxt = "[ 1LFD_Ral_CA ]"
n1atoms = 0
for line in LFDlines[2:-1] : 
    if "CA" in line : 
        resid = int(line[0:5])
        if resid > 13 and resid < 101 : 
            ndxtxt += "\n%i"%int(line[15:20])
            n1atoms += 1

try :
    gro = sys.argv[1]
    grofile = open(gro)
    grolines = grofile.readlines()
    grofile.close()
except :
    print "USAGE: %s <gro>"%sys.argv[0]
    sys.exit()

ndxtxt += "\n[ %s_Ral_CA ]"%gro
n2atoms = 0
for line in grolines[2:-1] : 
    if "CA" in line : 
        resid = int(line[0:5])
        if resid > 10 and resid < 98 :
            ndxtxt += "\n%i"%(int(line[15:20]))
            n2atoms += 1

if n1atoms != n2atoms : 
    print "ERROR:  Number of atoms in the 2 selections are different!"
    print " %10i : %s"%(n1atoms, LFD)
    print " %10i : %s"%(n2atoms, gro)
    sys.exit()
print ndxtxt

