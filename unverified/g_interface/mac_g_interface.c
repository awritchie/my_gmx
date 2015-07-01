/*
 * $Id: template.c,v 1.5 2008/05/29 08:36:53 hess Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: template.c,v 1.5 2008/05/29 08:36:53 hess Exp $";

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include <math.h>
#include "atomprop.h"
#include <stdio.h>
#include <stdlib.h>

bool distance( double threshold, t_trxframe fr, int frame_index, int region_index, t_topology top, bool bVerbose )
{
    double separation = pow(
                            pow(fr.x[frame_index][XX]-fr.x[region_index][XX],2)
                            +
                            pow(fr.x[frame_index][YY]-fr.x[region_index][YY],2)
                            +
                            pow(fr.x[frame_index][ZZ]-fr.x[region_index][ZZ],2)
                            ,.5);

    if (separation <= threshold) {
        if (bVerbose) {
            printf("\nSeparation=%8.5f\n",separation);
            printf("Coordinates at %i%s %s : %8.5f %8.5f %8.5f\n",
                   top.atoms.atom[frame_index].resnr+1,
                   *top.atoms.resname[top.atoms.atom[frame_index].resnr],
                   *top.atoms.atomname[frame_index],
                   fr.x[frame_index][XX],
                   fr.x[frame_index][YY],
                   fr.x[frame_index][ZZ]);
            printf("Coordinates at %i%s %s : %8.5f %8.5f %8.5f\n",
                   top.atoms.atom[region_index].resnr+1,
                   *top.atoms.resname[top.atoms.atom[region_index].resnr],
                   *top.atoms.atomname[region_index],
                   fr.x[region_index][XX],
                   fr.x[region_index][YY],
                   fr.x[region_index][ZZ]);

        }
        return TRUE;
    }
    else {
        return FALSE;
    }
}

double compute_distance( rvec coords1, rvec coords2 ) 
{
    double separation = pow(
                            pow(coords1[0]-coords2[0],2)
                            +
                            pow(coords1[1]-coords2[1],2)
                            +
                            pow(coords1[2]-coords2[2],2)
                            ,.5);
    return separation;
}

int printMatrix2D( double array[][3], int length);
int printMatrix1D( double array[], int length);
void ATA( double input[][3], int length, double output[3][3]);
void ATB( double input[][3], int length, double output[3]);
void Minv( double input[3][3], double output[3][3]);
double detM( double input[3][3]);
void LeastSquares( double coords[][3], int length, double normal[3], bool bVerbose);
void cross( rvec a, rvec b, rvec result );
void projection( rvec a, rvec b, rvec result );
double angle_periodicity( double angle_in );


int main(int argc,char *argv[])
{
    static char *desc[] = {
      "    This utility is used to calculate an azimuthal and polar angle of a residue (first entry in index file, important only for cutoffs) sidechain (vector start and ending atom indices given with -a1 and -a2 flags, respectively) relative to an interfacial plane.  This plane is an average of two planes.  The first plane is determined by taking a least squares fit of group 1 atoms within x distance (-if1) of the residue.  The second plane is determined by taking a least squares fit of group 2 atoms within y distance (-if2) of the residue.  This is done in hopes of smoothing out the protein topology relative to its overall features.  I am hopeful that this proceedure will work out well enough...\n",
      "    The polar angle is determined by first taking a least squares fit of all CA in the structure file to form a vertical plane. The cross product of the vertical plane and the normal to the average plane is the the polar axis, rotated 90 degrees (to be consistant with Ensign et al.).  The polar angle then becomes the angle between the projection of the sidechain bondvector onto the surface plane and the rotated polar axis.",  
      "    The sign of the angle is determined suched that positive is pointed toward the first index group and negative is pointed away from the first index group.  Keep this in mind should you not want an average plane and instead use the same index group twice.\n",
      "    In addition, the number of water molecules in the first hydration shell (sum of -probe and -sphere flags) are also counted.\n"
      "    Despite what the -verbose flag says about the intermediate matrix calculations, it is not actually calculating A.T*A or A.T*B because matrix A it is refering to (to be consistant with azimop.py) is actually:\n",
      "  | x1 y1 1 |\n",
      "  | x2 y2 1 |\n",
      "  | x3 y3 1 |\n",
      "  |   ...   |\n",
      "and matrix B is :\n",
      "  | z1 |\n",
      "  | z2 |\n",
      "  | z3 |\n",
      "  | .. |\n",
      "While here, the only matrix given is :\n",
      "  | x1 y1 z1 |\n",
      "  | x2 y2 z2 |\n",
      "  | x3 y3 z2 |\n",
      "  |   ...    |\n",
      "     Instead, the functions used here actually computer A.T*A and A.T*B for the single matrix, looking at the elements as they would appear in the format of two matrices A and B.  The calculation is done in this way because we are doing a least-squares calculation fitting to the equation Ax+By+Cz+D=0, where we are setting C = -1 and solving for A, B, and D.  If verbose is flaggd, all of the middle-calculations will be shown and the text will slightly misleading without reading this. "
    };
  
    static real     solsize = 0.14, vdwsize=0.18;
    static real     if1_dist = 1.0;
    static real     if2_dist = 0.80;
    static int      n_atoms, n_interface1, n_interface2;
    atom_id         *ind_atoms, *ind_interface1, *ind_interface2;
    char            *gn_atoms, *gn_interface1, *gn_interface2;
    rvec            bondvector;
    t_atoms         *atoms=NULL, *interface1atoms=NULL, *interface2atoms=NULL;
    int             a1=-1, a2=-1;    // Initializing these negative as a check
    bool            bVerbose=FALSE, bVeryVerbose=FALSE, do_water_count=FALSE, do_gnuplot=FALSE, singleindex=FALSE, sep_files=FALSE;
    int             framenumber=0;
    static real     centering=2.5;
    double          phi, angle, reference_angle, nitrile_angle, phi2;

  
    t_pargs pa[] = {
        { "-a1", TRUE, etINT, 
            {&a1}, "Atom number 1 (must be set!)"},
        { "-a2", TRUE, etINT, 
            {&a2}, "Atom number 2 (must be set!)"},
        { "-if1", FALSE, etREAL, {&if1_dist},
            "Distance threshold for phase1 plane (nm)"},
        { "-if2", FALSE, etREAL, {&if2_dist},
            "Distance threshold for phase2 plane (nm)"},
        { "-probe", FALSE, etREAL, {&solsize},
            "Radius of the solvent probe (nm)"},
        { "-sphere", FALSE, etREAL, {&vdwsize},
            "Size of the hydration shell radius, less probe size"},
        { "-count", FALSE, etBOOL, {&do_water_count},
            "Count the number of waters within probe+sphere radius of selected group"},
        { "-v", FALSE, etBOOL, {&bVerbose},
            "Be slightly more verbose"},
        { "-vv", FALSE, etBOOL, {&bVeryVerbose},
            "Be much more verbose"},
        { "-gnuplot", FALSE, etBOOL, {&do_gnuplot},
            "Print out coordinate files and a file to load into gnuplot for visualization.  There is currently no smart way to change the (many) output file names.  This is intended purely for quick and easy visualization using gnuplot>load <filename>"},
        { "-center", FALSE, etREAL, {&centering},
            "Center +/- this value on all atoms contained in phase 2 (preferably the protein)"},
        { "-sep", FALSE, etBOOL, {&sep_files},
            "Write azimuthal and polar angles to separate .xvg files"}
    };
  
    if (solsize < 0) {
        solsize=1e-3;
        fprintf(stderr,"Probe size too small, setting it to %g\n",solsize);
    }
  
    t_topology top;
    int        ePBC;
    char       title[STRLEN];
    t_trxframe fr;
    rvec       *xtop;
    matrix     box;
    int        status;
    int        flags = TRX_READ_X;

    t_filenm fnm[] = {
      { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
      { efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
      { efXVG, "-o", "angles", ffWRITE },   /* azimuthal and polar angles xvg */
      { efXVG, "-o2", "polar", ffWRITE },   /* write polar angles to a separate file */
      { efXVG, "-c", "h2o_counts", ffOPTWR }, /* water count xvg */
      { efNDX, NULL, NULL, ffREAD } /* index file for knowing where to look */
    };
  
#define NFILE asize(fnm)

    CopyRight(stderr,argv[0]);

    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
                NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

    /* We don't need any topology information to write the coordinates,
     * but to show how it works we start by writing the name and
     * charge of the selected atom. It returns a boolean telling us
     * whether the topology was found and could be read
     */
  
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
    sfree(xtop);
    
    if (bVeryVerbose) {
        bVerbose=TRUE;
    }
    atoms=&top.atoms;

    printf("Select group containing atoms to compute cutoff distances from (to form interfaces and, if -count flagged, count water oxygens).  The normal vector to the plane should point towards this group if intelligently chosen:\n");
    /* n_atoms is the number of atoms, ind_atoms is the atoms indices */
    get_index(atoms,ftp2fn(efNDX,NFILE,fnm),1,&n_atoms,&ind_atoms,&gn_atoms);
    printf("Select group containing first group's atoms (above plane):\n");
    get_index(interface1atoms,ftp2fn(efNDX,NFILE,fnm),1,&n_interface1,&ind_interface1,&gn_interface1);
    printf("Select group containing second group's atoms (below plane):\n");
    get_index(interface2atoms,ftp2fn(efNDX,NFILE,fnm),1,&n_interface2,&ind_interface2,&gn_interface2);
    if (strncmp(gn_interface1,gn_interface2,20)==0) {
        singleindex=TRUE;
        printf("\nUsing protein CA Center of Mass to determine surface normal direction.\nNOTE: This assumes the CA Center of Mass is sufficiently buried relative to the locations of atoms -a1 and -a2.\n");
    }
  
    /* The first time we read data is a little special */
    read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
    
    /* check for required parameters */
    if( a1<1 || a2<1 ) 
        gmx_fatal(FARGS, "Atom numbers a1 and a2 defining the bond vector must be specified\n" );
    a1--;a2--; //internally, numbering starts at 0
    /* check that these atoms exists */
    if(a1<0 || a1>(top.atoms.nr)) 
    {
        printf("Error: Atom number %d is out of range.\n",a1);
        exit(1);
    }
    if(a2<0 || a2>(top.atoms.nr)) 
    {
        printf("Error: Atom number %d is out of range.\n",a2);
        exit(1);
    }
   
    FILE *wateroutputfile, *angleoutputfile, *polaroutputfile;
    angleoutputfile=ffopen( opt2fn("-o", NFILE,fnm),"w");
    if (angleoutputfile==NULL) {
        gmx_fatal(FARGS,"Failed to open output file, '%s'\n",opt2fn("-o",NFILE,fnm) );
    }
    if (sep_files) {
        polaroutputfile=ffopen( opt2fn("-o2",NFILE,fnm),"w");
        if (polaroutputfile==NULL) {
            gmx_fatal(FARGS,"Failed to open polar output file, '%s'\n", opt2fn("-o2",NFILE,fnm) );
        }
    }
    if (do_water_count) {
        wateroutputfile=ffopen( opt2fn("-c",NFILE,fnm),"w");
        if (wateroutputfile==NULL) {
            gmx_fatal(FARGS,"Failed to open water count output file, '%s'\n",opt2fn("-c",NFILE,fnm) );
        }
    }
    /* This is the main loop over frames */
    do {
        /* Bond vector we're looking at */
        rvec_sub(fr.x[a2],fr.x[a1],bondvector);
        unitv( bondvector, bondvector ); 
        
        double phase1_coords[n_interface1][3]; /* The number of kept coordinates should be less than or equal to the number of tested coordinates */
        memset( phase1_coords, 0, 3*n_interface1*sizeof(double) );
        double phase2_coords[n_interface2][3]; /* Same logic */ 
        memset( phase2_coords, 0, 3*n_interface2*sizeof(double) );
        int atoms1 = 0, atoms2 = 0, waters = 0;
        
        /* interface1 loop */
        double com[3] = {0};
        double totalM = 0;
        for (int n=0 ; n < n_interface1 ; n++ ) {
            for (int i=0 ; i < n_atoms ; i++ ) {
                if (distance( if1_dist,fr,ind_interface1[n],ind_atoms[i],top,bVeryVerbose) ) {
                    phase1_coords[atoms1][0]=fr.x[ind_interface1[n]][XX];
                    phase1_coords[atoms1][1]=fr.x[ind_interface1[n]][YY];
                    phase1_coords[atoms1][2]=fr.x[ind_interface1[n]][ZZ];
                    totalM += top.atoms.atom[ind_interface1[n]].m;
                    com[0] += phase1_coords[atoms1][0]*top.atoms.atom[ind_interface1[n]].m;
                    com[1] += phase1_coords[atoms1][1]*top.atoms.atom[ind_interface1[n]].m;
                    com[2] += phase1_coords[atoms1][2]*top.atoms.atom[ind_interface1[n]].m;
                    atoms1++;
                    break;
                }
            }
        }
        
        for (int i=0;i<3;i++) {
            com[i]=com[i]/totalM;
        }

        if (bVerbose) {
            printf("\n%i atoms as part of the first interface\n", atoms1);
        }
        
        /* interface2 loop */
        /* defining com2 and totalM2 for making gnuplot figures.  This gives us the center of mass to center the figure on. */
        double com2[3] = {0};
        double totalM2 = 0;
        for (int n=0 ; n < n_interface2 ; n++ ) {
            for (int i=0 ; i < n_atoms ; i++ ) {
                if (distance( if2_dist,fr,ind_interface2[n],ind_atoms[i],top,bVeryVerbose) ) {
                    phase2_coords[atoms2][0]=fr.x[ind_interface2[n]][XX];
                    phase2_coords[atoms2][1]=fr.x[ind_interface2[n]][YY];
                    phase2_coords[atoms2][2]=fr.x[ind_interface2[n]][ZZ];
                    totalM2 += top.atoms.atom[ind_interface2[n]].m;
                    com2[0] += fr.x[ind_interface2[n]][XX]*top.atoms.atom[ind_interface2[n]].m;
                    com2[1] += fr.x[ind_interface2[n]][YY]*top.atoms.atom[ind_interface2[n]].m;
                    com2[2] += fr.x[ind_interface2[n]][ZZ]*top.atoms.atom[ind_interface2[n]].m;
                    atoms2++;
                    break;
                }
            }
        }
        
        for (int i=0;i<3;i++) {
            com2[i]=com2[i]/totalM2;
        }

        if (bVerbose) {
            printf("\n%i atoms as part of the second interface\n", atoms2);
        }
        
        /* CA loop.  The polar axis is ALWAYS defined by a CA plane */
        double protein_coords[top.atoms.nr][3]; /* Should be less CA than total atoms*/
        memset( protein_coords, 0, 3*top.atoms.nr*sizeof(double) );
        int caatoms=0;
        double comca[3] = {0};
        double totalcaM = 0;
        for (int i=0 ; i < top.atoms.nr ; i++ ) {
            if (strncmp(*top.atoms.atomname[i],"CA",3) ==0 ) { /* 3 because I don't want something like CAG matching*/
                protein_coords[caatoms][0]=fr.x[i][XX];
                protein_coords[caatoms][1]=fr.x[i][YY];
                protein_coords[caatoms][2]=fr.x[i][ZZ];
                totalcaM += top.atoms.atom[i].m;
                comca[0] += fr.x[i][XX]*top.atoms.atom[i].m;
                comca[1] += fr.x[i][YY]*top.atoms.atom[i].m;
                comca[2] += fr.x[i][ZZ]*top.atoms.atom[i].m;
                caatoms++;
            }
        }
        for (int i=0;i<3;i++) {
            comca[i]=comca[i]/totalcaM;
        }

        /* counting waters */
        if (do_water_count) {
            if (bVeryVerbose) {printf("\nCounting Waters\n");}
            for (int n=0; n < top.atoms.nr ; n++ ) { 
                for (int i=0 ; i < n_atoms ; i++ ) {
                    if (strncmp(*top.atoms.atomname[n],"OW",3) ==0 ) {
                        if (distance( solsize+vdwsize,fr,n,ind_atoms[i],top,bVeryVerbose)) {
                            waters++;
                            break;
                        }
                    }
                }
            }
            fprintf(wateroutputfile,"%i\n",waters);
        }
        
        /* Compute a least squares fit plane for each set of coordinates */
        double normalvector1[3] = {0}, normalvector2[3] = {0}, systemnormal[3] = {0};

        LeastSquares(protein_coords, caatoms, systemnormal, bVerbose);
        LeastSquares(phase1_coords, atoms1, normalvector1, bVerbose);
        LeastSquares(phase2_coords, atoms2, normalvector2, bVerbose);

        /* Build the average plane and prepare to check the normal vector direction wrt the system */
        rvec average={(normalvector1[0]+normalvector2[0])/2,(normalvector1[1]+normalvector2[1])/2,-1};
        unitv( average, average );
        rvec naverage ;
        double daverage[3];
        
        for (int i=0; i<3 ; i++){
            naverage[i]=-average[i];
            daverage[i]=(normalvector1[i]+normalvector2[i])/2; /* I want the equation in the form Ax+By+D=z */
        }
        
        /* Build polar axis; I don't want this changing based on whether the +/- normal vector is used, so it's before that check */
        rvec vertplane = {systemnormal[0], systemnormal[1], -1};
        rvec polaraxis;
        cross(average, vertplane, polaraxis);
        unitv( polaraxis, polaraxis );
        
        /* I want the normal vector to point away from the protein surface */
        /* When the same index group is used for both plane-building groups, pick whichever is higher, || proteinCoM-(a2 +/- normalvector) || */
        if ( singleindex ) { 
            rvec proteinCoM = {comca[0], comca[1], comca[2]};
            rvec nplus, nminus;
            rvec_add( fr.x[a2], average, nplus);
            rvec_sub( fr.x[a2], average, nminus);
            if ( compute_distance(proteinCoM, nplus) < compute_distance(proteinCoM, nminus) ) {
                for (int i=0; i<3 ; i++){
                    average[i]=naverage[i];
                }
                if (bVerbose) {
                    printf("\nUsing negative of calculated normal vector\n");
                }
            }
        }
        else {
            /* When two different index groups are used for plane-building, pick whichever is higher, || bottomCoM-(topCoM +/- normalvector) || */
            rvec centerofmassbottom = {com2[0], com2[1], com2[2]};
            rvec nplus = {com[0]+average[0], com[1]+average[1], com[2]+average[2]};
            rvec nminus = {com[0]-average[0], com[1]-average[1], com[2]-average[2]};
            if ( compute_distance(centerofmassbottom, nplus) < 
                compute_distance(centerofmassbottom, nminus) ) {
                for (int i=0; i<3 ; i++){
                    average[i]=naverage[i];
                }
                if (bVerbose) {
                    printf("\nUsing negative of calculated normal vector\n");
                }
            }
        }
        
        if (bVerbose) {
            printf("\nAverage plane:\n");
            printMatrix1D(daverage,3);
            printf("\nNormal Vector:\n");
            for (int i=0; i<3;i++){
                printf("%f\n",average[i]);
            }
            printf("\nPolar Axis:\n");
            for (int i=0;i<3;i++){
                printf("%f\n",polaraxis[i]);
            }
            printf("\nBond Vector:\n");
            for (int i=0; i<3;i++){
                printf("%f\n",bondvector[i]);
            }
        }

        /* Elevation angle (azimuthal) = 90 - arccos(dot(plane normal, bond vector)) * 180/pi */
        unitv( average, average );
        double cosine = iprod( bondvector, average );
        double angle;
        if (cosine>=1) {
            printf("\n\nWARNING: FRAME %i, COSINE = %f\n\n",framenumber, cosine);
            angle=0.0;
        }
        else if (cosine <=-1){
            printf("\n\nWARNING: FRAME %i, COSINE = %f\n\n",framenumber, cosine);
            angle=180.0;
        }
        else {
            angle=90-acos(cosine)*180/M_PI;
        }
        
        /* Polar angle is a bit trickier since arccos gives an absolute value.  We need a reference vector too now */
        /* I'm going to redefine the surface plane in terms of vectors which span it, one of which is the polar axis.  I need to use the average plane (daverage), not normalized, to make the vector along the x-axis (y=0) */
        rvec xaxis = { 1, 0, daverage[0] };
        unitv( xaxis, xaxis );
        rvec reference_vector ;
        cross( polaraxis, average, reference_vector );

        /* Determine the sign on the angle between the polar axis and the x-axis vector based on the y value of the polar axis */
        if (polaraxis[1] >= 0) { reference_angle = acos( iprod(polaraxis, xaxis) )*180/M_PI;}
        else { reference_angle = -1*acos( iprod(polaraxis, xaxis) )*180/M_PI;}

        /* Project bondvector onto the polar axis */
        rvec proj_polar_bond;
        projection(polaraxis, bondvector, proj_polar_bond );

        /* Project bondvector onto the reference_vector */
        rvec proj_ref_bond;
        projection(reference_vector, bondvector, proj_ref_bond);
        
        /* Sum the 2 projections then convert to unit vector.  Each projection vector is the bond vector along an "Axis" made by the polar axis and the reference vector (which are perpenticular to each other by definition of cross product. */
        rvec projected_vector;
        rvec_add( proj_ref_bond, proj_polar_bond, projected_vector );
        unitv( projected_vector, projected_vector );
        
        /* Determine the sign on the angle between the projected vector and the x-axis based on the y value of the projected vector */
        if (projected_vector[1] >=0) { nitrile_angle = acos( iprod(xaxis, projected_vector) / 
                                                            ( iprod(xaxis, xaxis) * iprod(projected_vector, projected_vector)) )*180/M_PI;}
        else { nitrile_angle = -1*acos( iprod(xaxis, projected_vector) / 
                                    ( iprod(xaxis, xaxis) * iprod(projected_vector, projected_vector)) )*180/M_PI;}
        
        /* Angle about the polar axis is the nitrile angle about the x-axis less the angle between the x-axis and the polar axis.  90 degrees is subtracted (rotated) from this to be consistant with Dan's paper. */
        phi = nitrile_angle - reference_angle - 90;
        
        /* Let's make sure this is between -180 and 180 */
        phi = angle_periodicity(phi);
        
        /* Output to screen */
        if (bVerbose) {
            printf("\nAzimuthal, Polar:");
            printf("\n%f %f,\n",angle, phi);
        }
        
        /* save output */
        if (sep_files) {
            fprintf(angleoutputfile, "%f\n", angle);
            fprintf(polaroutputfile, "%f\n", phi);
        }
        else {
            fprintf( angleoutputfile, "%f %f\n", angle, phi );        
        }
        
        if (do_gnuplot) {
            FILE *p1c, *p2c, *p1ct, *p2ct, *cnvector, *gnufile; /* phase n coordinates and phase n coordinates total */
            char p1coords[128], p2coords[128], p1accept[128], p2accept[128], nitrilevector[128], gnuload[128];

            /* frame index number as part of the name */
            sprintf(p1coords, "phase1.%i.coords", framenumber);
            sprintf(p1accept, "phase1.%i.accept", framenumber);
            sprintf(p2coords, "phase2.%i.coords", framenumber);
            sprintf(p2accept, "phase2.%i.accept", framenumber);
            sprintf(nitrilevector, "nitrile.%i.coords", framenumber);
            sprintf(gnuload, "frame%i.gnu", framenumber);
            
            /* all phase 1 coords */
            p1ct=ffopen(p1coords, "w");
            for (int i=0; i<n_interface1; i++) {
                fprintf(p1ct,"%f %f %f\n",fr.x[ind_interface1[i]][XX],fr.x[ind_interface1[i]][YY],fr.x[ind_interface1[i]][ZZ]);
            }
            fclose(p1ct);
            
            /* accepted phase 1 coords */
            p1c=ffopen(p1accept, "w");
            for (int i=0; i<atoms1; i++) {
                fprintf(p1c,"%f %f %f\n",phase1_coords[i][0], phase1_coords[i][1], phase1_coords[i][2]);
            }
            fclose(p1c);
            
            /* all phase 2 coords */
            p2ct=ffopen(p2coords, "w");
            for (int i=0; i<n_interface2; i++) {
                fprintf(p2ct,"%f %f %f\n",fr.x[ind_interface2[i]][XX],fr.x[ind_interface2[i]][YY],fr.x[ind_interface2[i]][ZZ]);
            }
            fclose(p2ct);
            
            /* accepted phase 2 coords */
            p2c=ffopen(p2accept, "w");
            for (int i=0; i<atoms2; i++) {
                fprintf(p2c,"%f %f %f\n",phase2_coords[i][0], phase2_coords[i][1], phase2_coords[i][2]);
            }
            fclose(p2c);
            
            /* nitrile vector */
            cnvector=ffopen(nitrilevector, "w");
            fprintf(cnvector, "#CD %f %f %f\n", fr.x[a1][XX], fr.x[a1][YY], fr.x[a1][ZZ]);
            fprintf(cnvector, "#NE %f %f %f\n", fr.x[a2][XX], fr.x[a2][YY], fr.x[a2][ZZ]);
            fprintf(cnvector, "#nBV %f %f %f\n", bondvector[0], bondvector[1], bondvector[2]);
            fprintf(cnvector, "%f %f %f %f %f %f\n", fr.x[a1][XX], fr.x[a1][YY], fr.x[a1][ZZ], bondvector[0], bondvector[1], bondvector[2]);
            fclose(cnvector);
            
            /* gnuplot stuff */
            /* A check to make sure the vector atoms are within the centering range given */
            gnufile=ffopen(gnuload, "w");
            double d1=pow(
                         pow(com2[0]-fr.x[a1][XX],2)+
                         pow(com2[1]-fr.x[a1][YY],2)+
                         pow(com2[2]-fr.x[a1][ZZ],2)
                         ,.5);
            double d2=pow(
                         pow(com2[0]-fr.x[a2][XX],2)+
                         pow(com2[1]-fr.x[a2][YY],2)+
                         pow(com2[2]-fr.x[a2][ZZ],2)
                         ,.5);
            if ( d1 >= centering || d2 > centering ) {
                printf("\nWARNING: The first atom in the atom vector, %s, or the second atom in the atom vector, %s, is outside of your centering range.\nThey are %f and %f nm away, respectively.  Letting gnuplot decide on ranges.\n", *top.atoms.atomname[a1], *top.atoms.atomname[a2], d1, d2);
                fprintf(gnufile, "#angle=%.3f\nset term jpeg; \nset output 'frame%i.jpeg'; \nset ticslevel 0; \nset xrange[*:*]; \nset yrange[*:*]; \nset zrange[*:*]; \nset pointsize 1; \nset view 100,85; \nsplot \\\n%f*x+%f*y+%f lt 5 title 'average plane', \\\n%f*x+%f*y+%f lt 3 title 'phase 1 plane', \\\n%f*x+%f*y+%f lt 1 title 'phase 2 plane', \\\n\"%s\" with lines lw 3 lt 5 title 'protein', \\\n\"%s\" with points pt 7 lt 3 title 'accepted phase 1 atoms', \\\n\"%s\" with points pt 7 lt 1 title 'accepted phase 2 atoms', \\\n\"%s\" with vector lw 3 lt 4 title '%.3f'\n#set term x11; replot", \
                        angle, framenumber, daverage[0], daverage[1], daverage[2], normalvector1[0], normalvector1[1], normalvector1[2], normalvector2[0], normalvector2[1], normalvector2[2], p2coords, p1accept, p2accept, nitrilevector, angle);
            }
            else {
                fprintf(gnufile, "#angle=%.3f\nset term jpeg; \nset output 'frame%i.jpeg'; \nset ticslevel 0; \nset xrange[%f:%f]; \nset yrange[%f:%f]; \nset zrange[%f:%f]; \nset pointsize 1; \nset view 100,85; \nsplot \\\n%f*x+%f*y+%f lt 5 title 'average plane', \\\n%f*x+%f*y+%f lt 3 title 'phase 1 plane', \\\n%f*x+%f*y+%f lt 1 title 'phase 2 plane', \\\n\"%s\" with lines lw 3 lt 5 title 'protein', \\\n\"%s\" with points pt 7 lt 3 title 'accepted phase 1 atoms', \\\n\"%s\" with points pt 7 lt 1 title 'accepted phase 2 atoms', \\\n\"%s\" with vector lw 3 lt 4 title '%.3f'\n#set term x11; replot", \
                    angle, framenumber, com2[0]-centering, com2[0]+centering, com2[1]-centering, com2[1]+centering, com2[2]-centering, com2[2]+centering, daverage[0], daverage[1], daverage[2], normalvector1[0], normalvector1[1], normalvector1[2], normalvector2[0], normalvector2[1], normalvector2[2], p2coords, p1accept, p2accept, nitrilevector, angle);
            }
            fclose(gnufile);
        }
        framenumber++;
    } while(read_next_frame(status,&fr));
    
    if (do_water_count) {
        fclose(wateroutputfile);
    }
    
    
    fclose(angleoutputfile);
    thanx(stderr);
  
    return 0;
}

int printMatrix1D( double array[], int length)
{
    for (int i=0; i<length; i++) {
        printf("%f\n", array[i]);
    }
    return 0;
}

int printMatrix2D( double array[][3], int length)
{
    for (int i=0; i<length; i++) {
        for (int j=0; j<3; j++){
            printf("%f ",array[i][j]);
        }
        printf("\n");
    }
    return 0;
}

void ATA( double input[][3], int length, double output[3][3]) 
{
    for (int i=0;i<length;i++){
        output[0][0] += input[i][0]*input[i][0];
        output[0][1] += input[i][0]*input[i][1];
        output[0][2] += input[i][0]*1;
        output[1][1] += input[i][1]*input[i][1];
        output[1][2] += input[i][1]*1;
        output[2][2] += 1*1;
    }
    /* matrix is symmetric */
    output[1][0] = output[0][1];
    output[2][0] = output[0][2];
    output[2][1] = output[1][2];
    
    return;
}

void ATB( double input[][3], int length, double output[3] )
{
    for (int i=0 ; i < length ; i++) {
        output[0] += input[i][0]*input[i][2];
        output[1] += input[i][1]*input[i][2];
        output[2] += input[i][2]*1;
    }
    return;
}

double detM( double input[3][3] )
{
    double determinant = 0;
    determinant += input[0][0]*(input[1][1]*input[2][2]-input[1][2]*input[2][1]);
    determinant -= input[0][1]*(input[1][0]*input[2][2]-input[1][2]*input[2][0]);
    determinant += input[0][2]*(input[1][0]*input[2][1]-input[1][1]*input[2][0]);
    return determinant;
}

void Minv( double input[3][3], double output[3][3] )
{
    double determinant = detM( input );
    if (determinant == 0 ){
        gmx_fatal(FARGS,"Determinant of coordinates is zero; consider changing the number of atoms in index file or the distance cutoffs" );
    }
    double transpose[3][3]={{0}};
    
    for (int i=0; i<3 ; i++){
        for (int j=0; j<3 ; j++){
            transpose[i][j]=input[j][i];
        }
    }
    
    double adjoint[3][3]={{0}};
    adjoint[0][0] += transpose[1][1]*transpose[2][2]-transpose[1][2]*transpose[2][1];
    adjoint[0][1] -= transpose[1][0]*transpose[2][2]-transpose[1][2]*transpose[2][0];
    adjoint[0][2] += transpose[1][0]*transpose[2][1]-transpose[1][1]*transpose[2][0];
    adjoint[1][0] -= transpose[0][1]*transpose[2][2]-transpose[0][2]*transpose[2][1];
    adjoint[1][1] += transpose[0][0]*transpose[2][2]-transpose[0][2]*transpose[2][0];
    adjoint[1][2] -= transpose[0][0]*transpose[2][1]-transpose[0][1]*transpose[2][0];
    adjoint[2][0] += transpose[0][1]*transpose[1][2]-transpose[0][2]*transpose[1][1];
    adjoint[2][1] -= transpose[0][0]*transpose[1][2]-transpose[0][2]*transpose[1][0];
    adjoint[2][2] += transpose[0][0]*transpose[1][1]-transpose[0][1]*transpose[1][0];
    
    for (int i=0; i<3 ; i++){
        for (int j=0; j<3 ; j++){
            output[i][j]=adjoint[i][j]/determinant;
        }
    }
    return;
}

void LeastSquares( double coords[][3], int length, double normal[3], bool bVerbose )
{
    /* I'm not actually calculating A.T*A or A.T*B because matrix A in the previous python script is actually:
            | x1 y1 1 |
            | x2 y2 1 |
            | x3 y3 1 |
            |   ...   |
     and matrix B is :
            | z1 |
            | z2 |
            | z3 |
            | .. |
     While here, the only matrix given is :
            | x1 y1 z1 |
            | x2 y2 z2 |
            | x3 y3 z2 |
            |   ...    |
     Instead, the functions used here actually computer A.T*A and A.T*B for the single matrix, looking at the elements as they would appear in the format of two matrices A and B.  If verbose is flaggd, all of the middle-calculations will be shown and the text will slightly misleading without reading this. */
    if ( length > 3 ) {
        double xy1Txy1[3][3] = {{0}};
        double xy1Tz[3] = {0};
        double xy1Txy1Inv[3][3] = {{0}};
    
        ATA( coords, length, xy1Txy1 );
        ATB( coords, length, xy1Tz );
        Minv( xy1Txy1, xy1Txy1Inv );

        for (int i=0; i<3 ; i++) {
            for (int j=0; j<3 ; j++){
                normal[i] += xy1Txy1Inv[i][j]*xy1Tz[j];
            }
        }
        if (bVerbose) {
            printf("\nLooking at matrix containing %i atoms :\n", length);
            printMatrix2D(coords,length);
            printf("\n  A.T*A :\n");
            printMatrix2D(xy1Txy1,3);
            printf("\n  (A.T*A).I :\n");
            printMatrix2D(xy1Txy1Inv,3);
            printf("\n  A.T*B :\n");
            printMatrix1D(xy1Tz,3);
            printf("\n  Plane Equation Ax+By+D=z :\n");
            printMatrix1D(normal,3);
        }
    }
    else {
                gmx_fatal(FARGS,"Fewer than 4 atoms were selected for forming a plane.  A minimum of 4 coordinates are required for least-square fitting; consider changing your cutoff settings." );
    }
    return;
}

void cross( rvec a, rvec b, rvec result ) {
    result[0] = a[1]*b[2]-a[2]*b[1];
    result[1] = a[2]*b[0]-a[0]*b[2];
    result[2] = a[0]*b[1]-a[1]*b[0];
    return;
}

void projection( rvec a, rvec b, rvec result ) {
    double scalarvalue = iprod(a,b)/(iprod(a,a));
    for (int i=0; i<3; i++) {
        result[i]=a[i]*scalarvalue;
    }
    return;
}
    
double angle_periodicity( double angle_in ) {
    double x = cos( angle_in*M_PI/180 );
    double y = sin( angle_in*M_PI/180 );
    return atan2(y,x)*180/M_PI;
}
