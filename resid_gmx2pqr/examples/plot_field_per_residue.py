#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

try :
    input=sys.argv[1]
except :
    print "Usage: %s <xvg>\nWhere xvg obtained from resid_gmx2pqr -or <xvg>"
    sys.exit()
f = open(input,'r')
flines = f.readlines()
f.close()

maxF = 10
dat = []
for line in flines :
    if not line.startswith("#") and not line.startswith("@") and not line.startswith(";") :
        linedat = []
        for each in line.split()[1:] :
            linedat.append(float(each))
        dat.append(linedat)

nlines = len(dat)
nresid = len(dat[0])

dat = np.array(dat).T

H = np.array([[16,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])  # added some commas and array creation code

fig = plt.figure()

ax = fig.add_subplot(111)
plt.imshow(dat,vmax=maxF,vmin=-1*maxF,aspect=nlines/nresid,interpolation='nearest',cmap=plt.get_cmap('RdBu'))
plt.colorbar()
plt.savefig('%s'%input.replace(".xvg",".pdf"),format='pdf')