#! /usr/bin/env python

import os
import sys
sys.path.append("/Users/ritchie/Utilities/python_scripts")
from structure_manip import *

ref=structure("1LFD.gro")
tar=structure("rap_aligned.gro")
ndx="ral.ndx"

fndx=open(ndx)
ndxlines=fndx.readlines()
fndx.close()


ref_ndx = []
tar_ndx = []
intar = False
for line in ndxlines[1:] :
    if intar :
        tar_ndx.append(int(line))
    elif "[" not in line :
        ref_ndx.append(int(line))
    if "rap_aligned" in line :
        intar = True

kabsch_alignment(ref,tar,ref_ndx,tar_ndx,ndx=True)