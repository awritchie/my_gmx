#! /anaconda/bin/python

import os
import sys
sys.path.append("/Users/ritchie/Utilities/python_scripts")
from structure_manip import *

try :
    ref_name = sys.argv[1]
    ref = structure(ref_name)
    gro_name = sys.argv[2]
    gro = structure(gro_name)
except :
    print "Usage: %s <reference structure> <simulated structure>"%(sys.argv[0])
    sys.exit()

s1,s2,n1,n2 = align_sequence(ref,gro)
#print n1
#print n2