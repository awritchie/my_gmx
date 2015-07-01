#! /usr/bin/python

import sys
import os
from numpy import *
from math import cos, sin, atan2
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-r", "--reference", dest="reference", help="A reference pdb structure.\n Default: rapral.pdb", metavar="REFERENCE PDB", default="/home/ritchie/Desktop/Rap1a/Rap_Ral/rapral.pdb")
parser.add_option("-f", "--filelist", dest="filelist", help="File list containing <rotamer angle>,<.xtc, .trr, ...>,<.gro>,<reference atomselect>,<trajectory atomselect>,<atom 1>,<atom 2>", metavar="FILE LIST FILE")
parser.add_option("-o", "--output", dest="output", help="File output name.\m Default: output", metavar="OUTPUTNAME", default="output")
parser.add_option("-q", "--quick", default=False, help="Do not realign structures, just recalculate from a previously generated outname.coord.## file", dest="quicktest")
parser.add_option("-A", "--chain1", dest="chain1", help="Residues on chain 1 to be used to calculate a plane.  At least 3 residues must be selected, by residue number per 1LFD and 1GUA numbering, within quotations.\n Default: 25 33 34 35 36 37 38 40 41 42", metavar="CHAIN 1 RESIDUES", default="25 33 34 35 36 37 38 40 41 42")
parser.add_option("-B", "--chain2", dest="chain2", help="Residues on chain 2 to be used to calculate a plane.  At least 3 residues must be selected, by residue number per 1LFD and 1GUA numbering, within quotations.\n Default: 18 20 27 28 29 30 31 32 33 52 54", metavar="CHAIN 2 RESIDUES", default="18 20 27 28 29 30 31 32 33 52 54")
parser.add_option("-R", "--rotamer", dest="rotamer", help="Rotamer angles to look at.  If choosing more than one angle, place selection in quotes.\n Default: 0 60 120 180 240 300", metavar="ROTAMER", default="0 60 120 180 240 300")
parser.add_option("-p", "--plane", dest="plane", help="Print the reference plane and polar axis only.  DO NOT ADD THIS FLAG IF YOU WANT TO OUTPUT!!!", default=False)
(options, args) = parser.parse_args()

def norm( v ) :
    return v/mag(v)
def mag( v ) :
    return dist( v, array((0,0,0)) )
def diff( a1, a2 ):
    return a1-a2
def dist( a1, a2) :
    a = diff (a1, a2)
    return ( dot(a,a) )**.5
def periodicity( angle ) :
    return atan2( sin( angle*pi/180 ),cos( angle*pi/180 ) )*180/pi

# The input for this is a .pdb file in which there exists a "TER" in between the FIRST and SECOND chains.  
plane=options.plane

try :
    referencepdb = options.reference
    pdbfile = open( referencepdb )
    pdblines = pdbfile.readlines()
    pdbfile.close()
except :
    print 'Problem with reference file occurred...'
    sys.exit()
try :    
    filelist = options.filelist 
    filelistfile = open( filelist )
    listlines = filelistfile.readlines()
    filelistfile.close()
except :
    print 'Problem with file list file occurred...'
    if not plane :
        sys.exit()
try :
    outputname = options.output
    quicktest = options.quicktest
except :
    print 'Problem with the output name occurred...'
    if not plane :
        sys.exit()
try :    
    rotamers=options.rotamer
except :
    print 'Problem with rotamer angle list occurred...'
    if not plane:
        sys.exit()

RapResiduesRawInput=options.chain1
RalResiduesRawInput=options.chain2

RapResidue = RapResiduesRawInput.split()
Nresidues = len(RapResidue)
RalResidue = RalResiduesRawInput.split()
Mresidues = len(RalResidue)


# Building the surface plane of the first protein (Rap in my pdb)
in1 = True
Rapsurfacecoords = []
Ralsurfacecoords = []
CAcoords = []
RapD = []
RalD = []
CAd = []
rappoints = open( "rappoints.txt", "w")
ralpoints = open( "ralpoints.txt", "w")
CA1point = open("CA1points.txt", "w")
CA2point = open("CA2points.txt", "w")
print "\n"
for line in pdblines : 
    if "TER" in line :
        in1 = False
    if "ATOM" in line and " CA " in line :
        residuenumber = int(line[23:26])
        x = float(line[31:38])
        y = float(line[38:46])
        z = float(line[46:55])
        b = int(1)
        CAcoords.append((x,y,b))
        CAd.append(z)
        if in1 :
            CA1point.write(str(x)+" "+str(y)+" "+str(z)+"\n")
            for i in arange(Nresidues) :
                if residuenumber == int(RapResidue[i]) :
                    Rapsurfacecoords.append((x,y,b))
                    RapD.append(z)
                    rappoints.write(str(x)+" "+str(y)+" "+str(z)+"\n")
                    resname = line[17:21] 
                    resindex = int(line[6:11].strip())
                    print "FIRST chain residues selected: %s %d %i %f %f %f" % (resname, residuenumber, resindex, x, y, z) 
        elif not in1 :
            CA2point.write(str(x)+" "+str(y)+" "+str(z)+"\n")
            for i in arange(Mresidues) :
                if residuenumber == int(RalResidue[i]) :
                    b = int(1)
                    Ralsurfacecoords.append((x,y,b))
                    RalD.append(z)
                    ralpoints.write(str(x)+" "+str(y)+" "+str(z)+"\n")
                    resname = line[17:21]   
                    resindex = int(line[6:11].strip())
                    print "SECOND chain residues selected: %s %d %i %f %f %f" % (resname, residuenumber, resindex, x, y, z) 
rappoints.close()
ralpoints.close()
CA1point.close()
CA2point.close()


A1 = matrix(Rapsurfacecoords)
B1 = matrix(RapD).T
NormalVector1 = (A1.T*A1).I*(A1.T*B1)
if plane :
    print "A.T*A:\n",(A1.T*A1)
    print "(A.T*A).I:\n",(A1.T*A1).I
    print "A.T*B:\n",(A1.T*B1)
    print "\nFor Ax + By + D = z, the equation of the surface plane for the FIRST chain is: \n",NormalVector1[0],"x +",NormalVector1[1],"y +",NormalVector1[2],"= z\n"

A2 = matrix(Ralsurfacecoords)
B2 = matrix(RalD).T
NormalVector2 = (A2.T*A2).I*(A2.T*B2)
if plane :
    print "A.T*A:\n",(A2.T*A2)
    print "(A.T*A).I:\n",(A2.T*A2).I
    print "A.T*B:\n",(A2.T*B2)
    print "\nFor Ax + By + D = z, the equation of the surface plane for the SECOND chain is: \n",NormalVector2[0],"x +",NormalVector2[1],"y +",NormalVector2[2],"= z\n"

# Averaging the two planes
NormalVector3 = ((NormalVector1[0]+NormalVector2[0])/2,(NormalVector1[1]+NormalVector2[1])/2,(NormalVector1[2]+NormalVector2[2])/2)
print "\nFor Ax + By + D = z, the average surface plane is: \n",NormalVector3[0],"x +",NormalVector3[1],"y +",NormalVector3[2],"= z\n"
normalvector = (float(NormalVector3[0]),float(NormalVector3[1]),float(-1)) #the normal is the vector A, B, -C, where C is the coefficient for z
N=-norm(normalvector)

# Building the vertical plane from all CA
A3 = matrix(CAcoords)
B3 = matrix(CAd).T
VertPlane = (A3.T*A3).I*(A3.T*B3)
print "\nFor Ax + By + D = z, the vertical plane is: \n",VertPlane[0],"x +",VertPlane[1],"y +",VertPlane[2],"= z\n"

# Determine the line orthogonal to both the vertical and surface planes 
a = ((float(NormalVector3[0]),float(NormalVector3[1]),-1))
b = ((float(VertPlane[0]),float(VertPlane[1]),-1))
cx = cross(a,b)
cx = norm(cx)


A = ((float(NormalVector3[1]),int(-1)),(float(VertPlane[1]),int(-1)))
B = (float(-1*NormalVector3[2]),float(-1*VertPlane[2]))
amat = matrix(A)
bmat = matrix(B).T
result = amat.I*bmat
print "The surface axis in parametric coordinates:\n(0",float(result[0]),float(result[1]),") + T (",cx[0],cx[1],cx[2],")\n\n"


# Redefining the surface plane in term of the vectors which it spans
PV1=norm(cx)
PV=cross(cx,N)
PV2=norm(PV)
p1=((0,0,float(NormalVector3[0]*0+NormalVector3[1]*0+NormalVector3[2])))
p2=((1,0,float(NormalVector3[0]*1+NormalVector3[1]*0+NormalVector3[2])))
xaxis=norm(array(p2)-array(p1))

if PV1[1] >=0 :
    referenceangle=arccos(dot(PV1,xaxis))*180/pi
if PV1[1] < 0 :
    referenceangle=-1*arccos(dot(PV1,xaxis))*180/pi


# Exit if all we wanted was the axis properties.
if plane :
    sys.exit()

print "Looking at rotamers: %s" % rotamers
if quicktest!=False :
    for rot in rotamers.split():
        try:
            CNCvectors = open("%s.coords.%s" % (outputname, rot) )
            CNCvectors.close()
            quicktest=True
            print "You have indicated that the %s.coords.<##> files already exist, skipping trajectory alignment..." % outputname
        except:
            quicktest=False
            print "%s.coords.%s does not exist, despite you choosing a quick test.  Performing trajectory alignment..." % (outputname, rot)
        break

if quicktest==False : 
    # This will align the trajectories to the reference pdb file
    for rot in rotamers.split() :
        print "Aligning %s degree rotamer trajectories to reference..." % rot
        alignment = "%s.movescript.%s.tcl" % (outputname, rot)
        movescript = open( alignment, "w" ) 
        movescript.write('mol new %s\n' % referencepdb ) 
        Firstns = True
        count = 0
        for line in listlines :
            count+=1
            rotamer, trajectory, gro, reference, simulated, vector1, vector2 = line.split(",")
            if rotamer == rot :
                if Firstns :
                    Firstns = False
                    movescript.write('mol new %s\n' % gro )
                movescript.write('mol addfile %s waitfor all 1 \n' % trajectory)
        movescript.write('set outfile [open "%s.coords.%s" w]\nset sel1 [atomselect 0 %s]\nset sel2 [atomselect 1 %s]\nset all [atomselect 1 "all"]\nset v1 [atomselect 1 %s]\nset v2 [atomselect 1 %s]\nset nf [molinfo 1 get numframes]\nfor {set i 0} {$i < $nf} {incr i} {\n     $sel2 frame $i \n     $all frame $i \n     set M [measure fit $sel2 $sel1 weight mass] \n     $all move $M \n     $v1 frame $i\n     $v2 frame $i\n     set V1 [lindex [ $v1 get {x y z}] 0]\n     set V2 [lindex [ $v2 get {x y z}] 0]\n     set simdata($i.r) [vecsub $V1 $V2] \n     puts $outfile "$simdata($i.r)"\n     puts $simdata($i.r)\n     }\nputs "Done writing coordinates, please exit"\nclose $outfile\nquit\n\n' % ( outputname, rot, reference, simulated, vector1, vector2 ) ) 
        movescript.close()

    # This will align the trajectories to the reference pbd, assuming the file list file is correct.
        os.system("/Applications/VMD1.9.1.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -dispdev text -e %s.movescript.%s.tcl > vmdout.log " % (outputname, rot) )

# Now we calculate the azimuthal angle.  We do this by taking the dot product of the vector normal to the plane and the CN vector, and subtract that from 90 degrees.
for rot in rotamers.split() :
    try:
        CNCvectors = open("%s.coords.%s" % (outputname, rot) )
        CNClines = CNCvectors.readlines()[1:]
        CNCvectors.close()
        azimuthalfile = "%s.azimuthalangle.%s.xvg" % (outputname, rot )
        azimuthalangle = open( azimuthalfile, "w" )
    except:
        print "Error opening either %s.coords.%s or %s.azimuthalangle.%s.xvg, please consider checking your script execution input" % (outputname, rot, outputname, rot)
        sys.exit()

    for line in CNClines : 
        x, y, z = line.split()
        vector = ((float(x),float(y),float(z)))
        B = norm(vector)
        theta = float(90) - arccos((dot(N,B)))*180/pi
        azimuthalangle.write("%f\n" % theta )
    azimuthalangle.close()

# Now to calculate the polar angle.
for rot in rotamers.split() :
    CNCvectors = open("%s.coords.%s" % (outputname, rot) )
    CNClines = CNCvectors.readlines()[1:]
    CNCvectors.close()
    polarfile = "%s.polarangle.%s.xvg" % (outputname, rot)
    polarangle = open( polarfile, "w")
    for line in CNClines : 
        x, y, z = line.split()
        vector = ((float(x),float(y),float(z)))
        B = norm(vector)
        projwB = (dot(PV1,B)/dot(PV1,PV1))*PV1+(dot(PV2,B)/dot(PV2,PV2))*PV2
        ProjwB = norm(projwB)
        phi2 = arccos((dot(PV2,ProjwB)))*180/pi       
        if ProjwB[1] >= 0 :
            CNangle = arccos(dot(xaxis,ProjwB)/(mag(xaxis)*mag(ProjwB)))*180/pi
        if ProjwB[1] < 0 :
            CNangle = -1*arccos(dot(xaxis,ProjwB)/(mag(xaxis)*mag(ProjwB)))*180/pi
        phi = (CNangle - referenceangle) - 90
        #phi=periodicity(phi)
        
#        polarangle.write("%f %f\n" % (phi,phi2))
        polarangle.write("%f\n" % phi )
    polarangle.close()
