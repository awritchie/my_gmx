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


int printMatrix2D( double array[][3], int length);
int printMatrix1D( rvec array, int length);



int main(int argc,char *argv[])
{
    static char *desc[] = {
      " This utility calculates the angle of solvent exposure for a side chain, defined by a vector between two atoms.  The angle is the angle of elevation between the protein plane and the vector.  The protein plane is defined as the normal vector from the protein center of mass (determined by atoms given in the index file) and the given atom, -a3, of the side chain vector (intended to be the side chain's CA).  As such, this is treating the protein exterior topology in a very coarse manner, although that should be sufficient for out purposes here."
    };
  

//    static int      n_atoms;
//    atom_id         *ind_atoms;
//    char            *gn_atoms;
    rvec            bondvector, normalvector;
    t_atoms         *atoms=NULL;
    int             a1=-1, a2=-1, a3=-1; // Initializing these negative as a check
    bool            bVerbose=FALSE, do_gnuplot=FALSE;
    int             framenumber=0;
    static real     centering=3;


  
    t_pargs pa[] = {
        { "-a1", TRUE, etINT, 
            {&a1}, "Atom number 1 (must be set!)"},
        { "-a2", TRUE, etINT, 
            {&a2}, "Atom number 2 (must be set!)"},
        { "-a3", TRUE, etINT,
            {&a3}, "Atom Surface Normal points to (must be set!)"},
        { "-v", FALSE, etBOOL, {&bVerbose},
            "Be slightly more verbose"},
        { "-gnuplot", FALSE, etBOOL, {&do_gnuplot},
            "Print out coordinate files and a file to load into gnuplot for visualization.  There is currently no smart way to change the (many) output file names.  This is intended purely for quick and easy visualization using gnuplot>load <filename>"},
        { "-center", FALSE, etREAL, {&centering},
            "Center +/- this value on all atoms given index file.  Set this to 0 to let gnuplot decide on ranges."}
    };
  
  
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
      { efXVG, "-o", "angles", ffWRITE }//,   /* azimuthal and polar angles xvg */
//      { efNDX, NULL, NULL, ffREAD } /* index file for knowing where to look */
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

    
    atoms=&top.atoms;

    printf("Select group for calculating center of mass on:\n");
    /* n_atoms is the number of atoms, ind_atoms is the atoms indices */
//    get_index(atoms,ftp2fn(efNDX,NFILE,fnm),1,&n_atoms,&ind_atoms,&gn_atoms);

  
    /* The first time we read data is a little special */
    read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
    
    /* check for required parameters */
    if( a1<1 || a2<1 || a3<1 ) 
        gmx_fatal(FARGS, "Atom numbers a1, a2, and a3 defining the bond vector and surface point must be specified.\n" );
    a1--;a2--;a3--; //internally, numbering starts at 0
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
    if(a3<0 || a3>(top.atoms.nr)) 
    {
        printf("Error: Atom number %d is out of range.\n",a3);
        exit(1);
    }
   
    FILE *angleoutputfile;
    angleoutputfile=ffopen( opt2fn("-o", NFILE,fnm),"w");
    if (angleoutputfile==NULL) {
        gmx_fatal(FARGS,"Failed to open output file, '%s'\n",opt2fn("-o",NFILE,fnm) );
    }

    /* This is the main loop over frames */
    do {
        /* get the CoM */
        int atoms = 0;
        rvec com = {0,0,0};
        double totalM = 0;
        for (int i=0 ; i < top.atoms.nr ; i++ ) {
            if (strncmp(*top.atoms.atomname[i],"CA",3) ==0 ) { /* 3 because I don't want something like CAG matching*/            totalM += top.atoms.atom[i].m;
                com[0] += fr.x[i][XX]*top.atoms.atom[i].m;
                com[1] += fr.x[i][YY]*top.atoms.atom[i].m;
                com[2] += fr.x[i][ZZ]*top.atoms.atom[i].m;
                atoms++;
                if (bVerbose) {
                    printf("Index %i, residue %s %s %f %f %f\n", i+1, *top.atoms.atomname[i], *top.atoms.resname[i], fr.x[i][XX], fr.x[i][YY],fr.x[i][ZZ]);
                }
            }
        }
        
        for (int i=0;i<3;i++) {
            com[i]=com[i]/totalM;
        }
        
        if (bVerbose) {
            printf("\n%i CA found\n", atoms);
        }

        
        /* Bond vector we're looking at */
        rvec_sub(fr.x[a2],fr.x[a1],bondvector);
        unitv(bondvector, bondvector);
        /* Vector normal to protein plane */
        rvec_sub(fr.x[a3], com, normalvector);
        unitv(normalvector, normalvector);
        
        /* I still want to draw the plan going through the first atom in the vector */
        rvec daverage = {-normalvector[0]/normalvector[2], -normalvector[1]/normalvector[2], -1/normalvector[2]*(normalvector[0]*-fr.x[a1][XX]+normalvector[1]*-fr.x[a1][YY]+normalvector[2]*-fr.x[a1][ZZ])};
        
        if (bVerbose) {
            printf("\nAverage plane:\n");
            printMatrix1D(daverage,3);
        }


        double cosine = iprod( bondvector, normalvector );
        double angle;
        if (cosine>=1) {
            printf("\n\nWARNING: FRAME %i, COSINE = %f\n\n",framenumber, cosine);
            angle=90-0.0;
        }
        else if (cosine <=-1){
            printf("\n\nWARNING: FRAME %i, COSINE = %f\n\n",framenumber, cosine);
            angle=90-180.0;
        }
        else {
            angle=90-acos(cosine)*180/M_PI;
        }

        if (bVerbose) {
            printf("Cos Angle\n");
            printf("%f %f\n",cosine,angle);
        }
        
        /* save output */
        fprintf( angleoutputfile, "%f\n", angle );        
        
        if (do_gnuplot) {
            FILE *backbone, *cnvector, *gnufile;
            char bbcoords[128], nitrilevector[128], gnuload[128];

            /* frame index number as part of the name */
            sprintf(bbcoords, "backbone.%i.coords", framenumber);
            sprintf(nitrilevector, "nitrile.%i.coords", framenumber);
            sprintf(gnuload, "frame%i.gnu", framenumber);
          
            /* backbone coords */
            backbone=ffopen(bbcoords, "w");
            for (int i=0; i<top.atoms.nr; i++) {
                if (strncmp(*top.atoms.atomname[i],"CA",3)==0) {
                    fprintf(backbone, "%f %f %f\n", fr.x[i][XX],fr.x[i][YY],fr.x[i][ZZ]);
                }
            }
            fclose(backbone);
            
            /* nitrile vector */
            cnvector=ffopen(nitrilevector, "w");
            fprintf(cnvector, "#CD %f %f %f\n", fr.x[a1][XX], fr.x[a1][YY], fr.x[a1][ZZ]);
            fprintf(cnvector, "#NE %f %f %f\n", fr.x[a2][XX], fr.x[a2][YY], fr.x[a2][ZZ]);
            fprintf(cnvector, "#nBV %f %f %f\n", bondvector[0], bondvector[1], bondvector[2]);
            fprintf(cnvector, "%f %f %f %f %f %f\n", fr.x[a1][XX], fr.x[a1][YY], fr.x[a1][ZZ], bondvector[0], bondvector[1], bondvector[2]);
            //fprintf(cnvector, "%f %f %f %f %f %f\n", fr.x[a1][XX], fr.x[a1][YY], fr.x[a1][ZZ], normalvector[0], normalvector[1], normalvector[2]);
            fclose(cnvector);
            
            /* gnuplot stuff */
            /* A check to make sure the vector atoms are within the centering range given */
            
            gnufile=ffopen(gnuload, "w");
            double d1=pow(
                         pow(com[0]-fr.x[a1][XX],2)+
                         pow(com[1]-fr.x[a1][YY],2)+
                         pow(com[2]-fr.x[a1][ZZ],2)
                         ,.5);
            double d2=pow(
                         pow(com[0]-fr.x[a2][XX],2)+
                         pow(com[1]-fr.x[a2][YY],2)+
                         pow(com[2]-fr.x[a2][ZZ],2)
                         ,.5);
            if ( d1 >= centering || d2 > centering ) {
                printf("\nWARNING: The first atom in the atom vector, %s, or the second atom in the atom vector, %s, is outside of your centering range.\nThey are %f and %f nm away, respectively.  Letting gnuplot decide on ranges.\n", *top.atoms.atomname[a1], *top.atoms.atomname[a2], d1, d2);
                fprintf(gnufile, "#angle=%.3f\n\
                        set term jpeg; \n\
                        set output 'frame%i.jpeg'; \n\
                        set label 1 \"%i\" at %f,%f; \n\
                        set ticslevel 0; \n\
                        set xrange[*:*]; \n\
                        set yrange[*:*]; \n\
                        set zrange[*:*]; \n\
                        set pointsize 1; \n\
                        set view 100,85; \n\
                        splot \\\n%f*x+%f*y+%f lt 1 title 'protein plane', \\\n\
                        \"%s\" with lines lw 3 lt 3 title 'protein', \\\n\
                        \"%s\" with vector lw 3 lt 4 title '%.3f'\n#set term x11; replot", \
                        angle, framenumber, framenumber, daverage[0], daverage[1], daverage[0], daverage[1], daverage[2], bbcoords, nitrilevector, angle);
            }
            else {
                fprintf(gnufile, "#angle=%.3f\n\
                        set term jpeg; \n\
                        set output 'frame%i.jpeg'; \n\
                        set ticslevel 0; \n\
                        set xrange[%f:%f]; \n\
                        set yrange[%f:%f]; \n\
                        set zrange[%f:%f]; \n\
                        set pointsize 1; \n\
                        set view 100,85; \n\
                        splot \\\n%f*x+%f*y+%f lt 1 title 'protein plane', \\\n\
                        \"%s\" with lines lw 3 lt 3 title 'protein', \\\n\
                        \"%s\" with vector lw 3 lt 4 title '%.3f'\n#set term x11; replot", \
                        angle, framenumber, com[0]-centering, com[0]+centering, com[1]-centering, com[1]+centering, com[2]-centering, com[2]+centering, daverage[0], daverage[1], daverage[2], bbcoords, nitrilevector, angle);
            }
            fclose(gnufile);
        }

        framenumber++;
    } while(read_next_frame(status,&fr));
    
    fclose(angleoutputfile);
    thanx(stderr);
  
    return 0;
}

int printMatrix1D( rvec array, int length)
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


    