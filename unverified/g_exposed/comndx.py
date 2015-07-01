#! /usr/bin/python

import sys

def incoords ( line ) :
    try : 
        float(line.split()[-1])
        float(line.split()[-2])
        float(line.split()[-3])
        return True
    except :
        return False

def atomname( line, suffix ) :
    if suffix==".gro" :
        return line[11:15].split()[0]
    else :
        print "You have not yet set this up to look at %s files" %suffix
        sys.exit()
def index( line, suffix ) :
    if suffix==".gro" :
        return int(line[15:20].split()[0])
    else :
        print "You have not yet set this up to look at %s files" %suffix
        sys.exit()


try:
    input=sys.argv[1]
    suffix=input[-4:]
    output=input.replace(suffix, ".ndx")
    file=open(input)
    filelines=file.readlines()
    file.close()
    print output
except :
    print "Usage: %s <pdb or gro>" %sys.argv[0]

ndx=[]

for line in filelines :
    if incoords( line ) :
        if atomname(line,suffix) == "C" or \
           atomname(line,suffix) == "N" or \
           atomname(line,suffix) == "CA" or \
           atomname(line,suffix) == "O" or \:
            ndx.append( index(line,suffix) )

output=open(output, "w")
output.write("[ backbone, N-CA-C-O ]\n")
for i in range(len(ndx)) :
    output.write("%i\n" %ndx[i])
output.close()
        
