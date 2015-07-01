//
//  main.cpp
//  gmx2pqr
//
//  Created by Andrew Ritchie on 5/3/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "tpxio.h"
#include "index.h"
#include <math.h>
#include "atomprop.h"
#include <stdio.h>
#include <stdlib.h>
#include "my_structs.h"
#include "string.h"

gmx2amb_dat read_DAT( char * filename, bool bVerbose );
int gmx_gmx2pqr(int argc, const char * argv[]);
bool h_bonding( rvec nitrile, rvec atom, rvec oxygen, float hbond_cutoff, float angle_cutoff, bool bVerbose );

int gmx_gmx2pqr(int argc, const char * argv[])
{
    static char *desc[] = {
        "   This is supposed to turn a .xtc file into a series of .pqr and .in files for apbs.\n"
        "   Note: Changing the angle cutoff to < about 100 degrees will allow both hydrogens on a water to hydrogen bond.  This duplication of h-bonding residues is NOT watched for and as such changing this parameter should be done at your own risk."
    };
    
    rvec            bondvector;
    t_atoms         *atoms=NULL;
    int             a1=-1, a2=-1;   // Initializing these negative as a check
    float           hbond_cutoff=1, angle_cutoff=150.0;
    bool            bVerbose=0;     // Remember, 0=False, 1=True
    
    
    t_pargs pa[] = {
        { "-a1", TRUE, etINT,
            {&a1}, "Atom number 1 (CNC CD, must be set!)"},
        { "-a2", TRUE, etINT,
            {&a2}, "Atom number 2 (CNC NE, must be set!)"},
        { "-v", FALSE, etBOOL, {&bVerbose},
            "Be slightly more verbose"},
        { "-hb_length", FALSE, etREAL,
            {&hbond_cutoff}, "Maximum N-H h-bonding distance in nm"},
        { "-hb_angle", FALSE, etREAL,
            {&angle_cutoff}, "Minimum N-H-O h-bonding angle in degrees"}
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
        { efTRX, "-f", NULL, ffREAD },      /* This is for the trajectory */
        { efDAT, "-d", "AMBER.DAT", ffREAD}
    };
    
#define NFILE asize(fnm)
    
    CopyRight(stderr,argv[0]);
    
    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
    sfree(xtop);
    
    atoms=&top.atoms;
    
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
    
    fprintf(stderr,"\nWill read reference charges and radii from %s\n",(opt2fn("-d",NFILE,fnm)));
    read_DAT(opt2fn("-d",NFILE,fnm),bVerbose);

    /* This is the main loop over frames */
    
    int framen=0;
    hbond_array * hbonding=(hbond_array*)malloc(sizeof(hbond_array)*0);
    do {
        hbonding=realloc(hbonding,sizeof(hbond_array)*framen+1);
        hbonding[framen].n=0;
        hbonding[framen].resid=(int*)malloc(sizeof(int)*hbonding[framen].n);
        rvec_sub(fr.x[a2],fr.x[a1],bondvector); // Nitrile ond vector we're looking at 
        
//        fprintf(stderr,"%1.12f %1.12f %1.12f\n",bondvector[0],bondvector[1],bondvector[2]);
        for (int i=0 ; i < top.atoms.nr ; i++ )
        {
            if (strncmp(*top.atoms.atomname[i],"HW1",3) ==0) {
                if (h_bonding(fr.x[a2],fr.x[i],fr.x[i-1],hbond_cutoff,angle_cutoff, bVerbose)){
                    if (bVerbose) {fprintf(stderr,"%i %i%s %i%s",i,top.atoms.atom[i].resnr, *top.atoms.atomname[i],top.atoms.atom[i-1].resnr,*top.atoms.atomname[i-1]);}
                    hbonding[framen].resid=realloc(hbonding[framen].resid,sizeof(int)*hbonding[framen].n+1);
                    hbonding[framen].resid[hbonding[framen].n]=top.atoms.atom[i].resnr;
                    hbonding[framen].n++;
                }
            }
            else if (strncmp(*top.atoms.atomname[i],"HW2",3) ==0 ) {
                if (h_bonding(fr.x[a2],fr.x[i],fr.x[i-2],hbond_cutoff,angle_cutoff,bVerbose)){
                    if (bVerbose) {fprintf(stderr,"%i %i%s %i%s",i,top.atoms.atom[i].resnr, *top.atoms.atomname[i],top.atoms.atom[i-2].resnr,*top.atoms.atomname[i-2]);}
                    hbonding[framen].resid=realloc(hbonding[framen].resid,sizeof(int)*hbonding[framen].n+1);
                    hbonding[framen].resid[hbonding[framen].n]=top.atoms.atom[i].resnr;
                    hbonding[framen].n++;
                }
            }
            // I need to make HW1 and HW2 indistinguishable; ie: if it's hbonding to HW1 in 1 frame then HW2 in
            // the next frame, I want to count that as "persistant" since my time step is long enough that it's
            // more-or-lesss resonable to allow that flip.
            // I can do this by keeping residue IDs instead of atom indices.
        }

        if (bVerbose) {
            fprintf(stderr,"\nHBonding(%d): %d events, resids: ",framen,hbonding[framen].n);
            for (int i=0;i<hbonding[framen].n;i++) {
                fprintf(stderr,"%d ",hbonding[framen].resid[i]);
            }
        }
        framen++;
    } while(read_next_frame(status,&fr));

    hbond_array * persistant_hbonds=(hbond_array*)malloc(framen*sizeof(hbond_array));
    int total_bonding=0;
    bool first_occurance;
    for (int i=0; i<framen; i++) {
        first_occurance=TRUE;
        persistant_hbonds[i].n=0;
        persistant_hbonds[i].resid=(int*)malloc(sizeof(int)*persistant_hbonds[i].n);
        if (bVerbose) {fprintf(stderr,"\nFrame %d presistant h-bonds: ",i);}
        if (i == 0 ) { // I can only look at the next frame
            for (int j=0; j<hbonding[i].n ; j++) {
                for (int k=0; k<hbonding[i+1].n; k++) {
                    if (hbonding[i].resid[j] == hbonding[i+1].resid[k]){
                        if (bVerbose) {fprintf(stderr,"%i ",hbonding[i].resid[j]);}
                        persistant_hbonds[i].resid=realloc(persistant_hbonds[i].resid,sizeof(int)*persistant_hbonds[i].n+1);
                        persistant_hbonds[i].resid[persistant_hbonds[i].n]=hbonding[i].resid[j];
                        persistant_hbonds[i].n++;
                        if (first_occurance){total_bonding++;first_occurance=FALSE;}
                    }
                }
                
            }
        }
        else if ( i == framen -1 ) { // I can only look at the previous frame
            first_occurance=TRUE;
            for (int j=0; j<hbonding[i].n ; j++) {
                for (int k=0; k<hbonding[i-1].n; k++) {
                    if (hbonding[i].resid[j] == hbonding[i-1].resid[k]){
                        if (bVerbose) {fprintf(stderr,"%i ",hbonding[i].resid[j]);}
                        persistant_hbonds[i].resid=realloc(persistant_hbonds[i].resid,sizeof(int)*persistant_hbonds[i].n+1);
                        persistant_hbonds[i].resid[persistant_hbonds[i].n]=hbonding[i].resid[j];
                        persistant_hbonds[i].n++;
                        if (first_occurance){total_bonding++;first_occurance=FALSE;}
                    }
                }
            }
        }
        else { // I want to look forward as well as behind
            first_occurance=TRUE;
            for (int j=0; j<hbonding[i].n ; j++) {
                for (int k=0; k<hbonding[i+1].n; k++) {
                    if (hbonding[i].resid[j] == hbonding[i+1].resid[k]){
                        if (bVerbose) {fprintf(stderr,"%i ",hbonding[i].resid[j]);}
                        persistant_hbonds[i].resid=realloc(persistant_hbonds[i].resid,sizeof(int)*persistant_hbonds[i].n+1);
                        persistant_hbonds[i].resid[persistant_hbonds[i].n]=hbonding[i].resid[j];
                        persistant_hbonds[i].n++;
                        if (first_occurance){total_bonding++;first_occurance=FALSE;}
                    }
                }
                for (int k=0; k<hbonding[i-1].n; k++) {
                    if (hbonding[i].resid[j] == hbonding[i-1].resid[k]){
                        if (bVerbose) {fprintf(stderr,"%i ",hbonding[i].resid[j]);}
                        persistant_hbonds[i].resid=realloc(persistant_hbonds[i].resid,sizeof(int)*persistant_hbonds[i].n+1);
                        persistant_hbonds[i].resid[persistant_hbonds[i].n]=hbonding[i].resid[j];
                        persistant_hbonds[i].n++;
                        if (first_occurance){total_bonding++;first_occurance=FALSE;}
                    }
                }
            }
        }
    }
    free(hbonding);
    printf("Total frames with persistent h-bonding waters: %d\n",total_bonding); //printf instead of fprintf so i can bash > to a file
    // Now I have a list of water residues that participate in hydrogen bonding
    // Next I need to put together the water index.  I think this part is going to be tedious.

    return 0;
}

bool h_bonding( rvec nitrile, rvec atom, rvec oxygen, float hbond_cutoff, float angle_cutoff, bool bVerbose )
{
    double x2,y2,z2,distance;
    x2=(nitrile[XX]-atom[XX])*(nitrile[XX]-atom[XX]);
    y2=(nitrile[YY]-atom[YY])*(nitrile[YY]-atom[YY]);
    z2=(nitrile[ZZ]-atom[ZZ])*(nitrile[ZZ]-atom[ZZ]);
    distance=pow(x2+y2+z2,0.5);
    if ( distance <= hbond_cutoff ) {
        rvec v1,v2;
        rvec_sub(nitrile,atom,v1);
        rvec_sub(oxygen,atom,v2);
        unitv(v1,v1);
        unitv(v2,v2);
        double cosine = iprod(v1,v2);
        double angle;
        if ( cosine >=1 ) { angle=0.0;}
        else if (cosine <= -1) {angle=180.0;}
        else { angle = acos(cosine)*180/M_PI;}
        if ( angle >= angle_cutoff ) {
            if (bVerbose) {fprintf(stderr,"\n%f %f ",distance, angle);}
            return TRUE;
        }
    }
    return FALSE;
}

gmx2amb_dat read_DAT( char * filename, bool bVerbose )
{
    char line[1024];
    int nlines=0;
    FILE *f=fopen(filename,"r");
    while (fgets(line,1024,f))
    {
        nlines +=1 ;
    }
    gmx2amb_dat dat;
    gmx2amb* mygmx2amb = (gmx2amb*)malloc(sizeof(gmx2amb)*nlines);
    if (mygmx2amb == NULL) {
        fprintf(stderr,"\nError allocating memory in gmx2amb!");
        exit(1);
    }
    if ( f == NULL ) {
        fprintf(stderr,"\nError opening %s",filename);
        exit(1);
    }
    
    fclose(f);
    FILE *F=fopen(filename,"r");
    int n=0;
    while (fgets(line,1024,f))
    {
        char escapechar[1];
        char resname[10];
        char ambername[10];
        char gmxname[10];
        float radius, charge;
        sscanf(line,"%s",escapechar);
        if (strncmp(escapechar,"#",1) != 0) {
            sscanf(line,"%s %s %f %f %s", resname, gmxname, &charge, &radius, ambername);
            memcpy(mygmx2amb[n].gmxname,gmxname,strlen(gmxname)+1);
            memcpy(mygmx2amb[n].resname,resname,strlen(resname)+1);
            memcpy(mygmx2amb[n].ambername,ambername,strlen(ambername)+1);
            mygmx2amb[n].radius=radius;
            mygmx2amb[n].charge=charge;
            n++;
        }
    }
    if (bVerbose) {
        for (int i=0;i<n;i++){
            fprintf(stderr,"\n%i out of %i: gmxname = %s radius = %f charge = %f",i,n,mygmx2amb[i].gmxname,mygmx2amb[i].radius, mygmx2amb[i].charge);
        }
    }
    fclose(f);
    dat.dat=mygmx2amb;
    dat.n=nlines;
    return dat;
}

int main(int argc, const char * argv[])
{
    gmx_gmx2pqr(argc,argv);
    return 0;
}
