//
//  main.cpp
//  gmx2pqr
//
//  Created by Andrew Ritchie on 5/3/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>

#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/vec.h>
#include <gromacs/tpxio.h>
#include <gromacs/index.h>
#include <gromacs/txtdump.h>
#include "string.h"
#include "stdlib.h"

#include <math.h>
#include <stdio.h>
#include "my_structs.h"

gmx2amb_dat read_DAT( const char * filename, const char * atpname, gmx2amb_dat dat, gmx_bool bVerbose );

double dot( double a[3], double b[3] ) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double calculate_field( int istart, int iend, t_topology top, t_trxframe fr, rvec midpoint, rvec bondvector, double bondlength )
{
    double field[3],xyz[3];
    double r2, r, rr3, pfield=0;
    double permitivity_constant = 35950207149.4727056;
    double coulombs = 6.2415096516E18;
    double kbt = 298 * 1.380648813e-23;
    double cfac = 2.5e9 * permitivity_constant / ( coulombs * coulombs * kbt );
    for (int i=0;i<3;i++) { field[i] = 0; }
    
//    fprintf(stderr,"\n%.4f %.4f %.4f",midpoint[0],midpoint[1],midpoint[2]);
//    fprintf(stderr,"\nLooking at the field due to atoms %i to %i.\n",istart,iend);
    for (int i=istart;i<=iend;i++) {
        xyz[0] = midpoint[0] - fr.x[i][XX];
        xyz[1] = midpoint[1] - fr.x[i][YY];
        xyz[2] = midpoint[2] - fr.x[i][ZZ];
        for (int j=0;j<3;j++) { xyz[j] *=10; } // convert from nm to Angstroms
//        fprintf(stderr,"%.4f %.4f %.4f\n",xyz[0],xyz[1],xyz[2]);
        r2 = dot(xyz,xyz);
        r = pow(r2,.5);
        std::cout << r << std::endl;
        rr3 = 1.0 / ( r * r2 );
        for (int j=0;j<3;j++) {
            field[j] += rr3 * top.atoms.atom[i].q * xyz[j];
        }
    }

    // Project the field along the bond vector
    for (int i=0;i<3;i++) {
//        fprintf(stderr,"%.4f ",field[i]);
        pfield += field[i]*bondvector[i]/bondlength;
    }
    pfield *= cfac;
//    fprintf(stderr,"\n%.4f\n",pfield);
    return pfield;
}

int gmx_gmx2pqr(int argc, const char * argv[])
{
    const char *desc[] = {
        "\tThis convert from gromacs to APBS pqr files.",
        "Use the -a1 flag for the start of the bond vector",
        "and the -a2 flag for the end of the bond vector."
        "\n"
    };
    t_topology      top;
    t_atoms         *atoms = NULL, useatoms;
    int             ePBC;
    char            title[STRLEN];
    t_trxframe      fr;
    rvec            *xtop;
    matrix          box;
    t_trxstatus     *status;
    int             flags = TRX_READ_X;
    const char      *parm_name, *atp_name;
    output_env_t    oenv;
    const char      *tpr_file, *xtc_file, *index_file, *out_file = NULL;
    char            tpr_title[256], pqr[256];
    rvec            bondvector,midpoint;
    double          bondlength;
    int             a1=-1, a2=-1;   // Initializing these negative as a check
    gmx_bool        nn=FALSE;
    float           pdie=1;
    double          pfield,pcfield,pnfield,psfield;
    gmx_bool        bVerbose=FALSE;     // Remember, 0=False, 1=True
    static char     *a3=NULL;
    
    t_pargs pa[] = {
        { "-nn", TRUE, etBOOL,
            {&nn}, "Use a separate index file for each frame, IE: dynamic index files created through g_select."},
        { "-a1", TRUE, etINT,
            {&a1}, "Atom number 1 (CNC CD, must be set!)"},
        { "-a2", TRUE, etINT,
            {&a2}, "Atom number 2 (CNC NE, must be set!)"},
        { "-a3", TRUE, etSTR,
            {&a3},"String, in quotations, for atoms you want to know the field due to"},
        { "-pdie", TRUE, etREAL,
            {&pdie}, "Protein dielectric"},
        { "-v", FALSE, etBOOL, {&bVerbose},
            "Be slightly more verbose"}
    };
    
    t_filenm fnm[] = {
        { efTPS, NULL, NULL, ffREAD },
        { efTRX, "-f", NULL, ffREAD },
        { efDAT, "-d", "/Users/ritchie/Utilities/apbs/AMBER.DAT", ffREAD },
        { efPQR, "-o", NULL, ffWRITE },
        { efXVG, "-of", "field_projection", ffWRITE },
        { efNDX, NULL, NULL, ffREAD }
    };
    
#define NFILE asize(fnm)
#define NPA asize(pa)
    
    CopyRight(stderr,argv[0]);
    
    parse_common_args(&argc, argv,
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW |
                      PCA_TIME_UNIT | PCA_BE_NICE,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);
    tpr_file = ftp2fn(efTPS, NFILE, fnm);
    xtc_file = opt2fn("-f", NFILE, fnm);
    parm_name = opt2fn("-d", NFILE, fnm);
    atp_name = opt2fn("-a", NFILE, fnm);
    out_file = opt2fn("-o", NFILE,fnm);
    
    read_tps_conf(tpr_file, tpr_title, &top, &ePBC, &xtop, NULL, box, TRUE);
    sfree(xtop);
    atoms = &top.atoms;

    printf("Select group for to save atom coordinates:\n");
    index_file = ftp2fn(efNDX, NFILE, fnm);
//    int             *i_index = (int*)malloc(sizeof(int)*nn);
//    atom_id         **ind_index = (atom_id*)malloc(sizeof(atom_id)*nn);
//    char            **gn_index = (char*)malloc(sizeof(char)*nn);
    int i_index;
    atom_id *ind_index;
    char *gn_index;
    get_index(atoms,index_file,1,&i_index,&ind_index,&gn_index);

    
    /* The first time we read data is a little special */
    read_first_frame(oenv, &status, xtc_file, &fr, flags);
    
    /* check for required parameters */
    if( a1<1 || a2<1 )
        gmx_fatal(FARGS, "Atom numbers a1 and a2 defining the bond vector must be specified\n" );
    a1--;a2--; //internally, numbering starts at 0
    /* check that these atoms exists */
    if(a1<0 || a1>(top.atoms.nr))
    {
        fprintf(stderr,"Error: Atom number %d is out of range.\n",a1);
        exit(1);
    }
    if(a2<0 || a2>(top.atoms.nr))
    {
        fprintf(stderr,"Error: Atom number %d is out of range.\n",a2);
        exit(1);
    }
    
    /* Get all of the atoms indices for finding the field of */
    char *sep = " ,;:/'\"";
    char *word, *brkt;
    int field_due_to_atom[1000], n_field_atoms=0;
    if (a3) {
        for (word = strtok_r(a3, sep, &brkt); word; word = strtok_r(NULL, sep, &brkt))
        {
            int field_atom = atoi(word);
            if (field_atom>(top.atoms.nr))
            {
                fprintf(stderr,"Error: Atom number %d is out of range.\n",field_atom);
                exit(1);
            }
            field_due_to_atom[n_field_atoms] = field_atom - 1;
            n_field_atoms++;
        }
        printf("\nWill include field due to atoms: ");
        for (int i=0;i<n_field_atoms;i++){printf("%i ",field_due_to_atom[i]+1);}
        printf("\n");
    }
    
    fprintf(stderr,"\nWill read reference charges and radii from %s\n",(opt2fn("-d",NFILE,fnm)));
    gmx2amb_dat dat;
    dat = read_DAT(parm_name,atp_name,dat,bVerbose);
    
    /* This is the main loop over frames */
    int framen=0;
    int resid=0;
    double charge=0.0,radius;
    double radii[100000]; // lazy, I should make this dynamic
    char *atomname, *resname;
    char *pqr_file = out_file;
    pqr_file[strlen(pqr_file)-4] = 0;
    FILE *fout = ffopen( opt2fn("-of",NFILE,fnm),"w");
    fprintf(fout,";%13i-%i %18i %18i ",ind_index[0]+1,ind_index[i_index-1]+1,a1+1,a2+1);
    if (a3) {
        for (int i=0;i<n_field_atoms;i++){fprintf(fout,"%18i ",field_due_to_atom[i]+1); }
    }
    fprintf(fout,"%18s","F-excpt");
    do {
        if (nn && framen > 0) {
            get_index(atoms,index_file,1,&i_index,&ind_index,&gn_index);
        }
        
        sprintf(pqr,"%s%i.pqr",pqr_file,framen);
        FILE *pqrout = ffopen(pqr,"w");
        if ( pqrout == NULL ) {
            gmx_fatal(FARGS,"\nFailed to open output file, '%s'\n",pqr);
        }
        
        rvec_sub(fr.x[a2],fr.x[a1],bondvector); // Nitrile ond vector we're looking at
        bondlength = pow(iprod(bondvector,bondvector),.5);
        midpoint[0] = fr.x[a1][XX] + 0.5 * bondvector[0];
        midpoint[1] = fr.x[a1][YY] + 0.5 * bondvector[1];
        midpoint[2] = fr.x[a1][ZZ] + 0.5 * bondvector[2];

//        fprintf(stderr,"\n%1.12f %1.12f %1.12f : %1.12f\n",bondvector[0],bondvector[1],bondvector[2],bondlength);
        
        double x,y,z;
        int npoints = 50;
        for (int n=0; n < npoints+1 ; n++ )
        {
            x = fr.x[a1][XX]*10 + n*bondvector[0]*10/npoints;
            y = fr.x[a1][YY]*10 + n*bondvector[1]*10/npoints;
            z = fr.x[a1][ZZ]*10 + n*bondvector[2]*10/npoints;
            fprintf(pqrout,"ATOM %6i %4s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,"DUM","CNC",0,x,y,z,0.0,0.0);
        }

        pcfield = calculate_field(a1,a1,top,fr,midpoint,bondvector,bondlength) / pdie;
        pnfield = calculate_field(a2,a2,top,fr,midpoint,bondvector,bondlength) / pdie;
        pfield = 0;
//        fprintf(stderr,"\nfields, pdie=%.2f: %.4f %.4f\n",pdie,calculate_field(a1,a1,top,fr,midpoint,bondvector,bondlength)/pdie,calculate_field(a1,a1,top,fr,midpoint,bondvector,bondlength));

        // Let's get the radii from the first frame, then use the same radii for every other frame.
        // I don't see a significant performance increase relative to re-matching every frame.
        if ( framen == 0 ) {
            for (int n=0; n < i_index ; n++ ) {
                radii[n] = -1000. ;
                int i = ind_index[n];
                resid = top.atoms.atom[i].resind;
                charge = top.atoms.atom[i].q;
                atomname = *top.atoms.atomname[i];
                resname = *top.atoms.resinfo[resid].name;
                for ( int j=0; j < dat.n ; j++ ) {
                    if ( bVerbose) {fprintf(stderr,"\nTrying to match: %s %s to %s (%d/%d)",resname,*top.atoms.atomtype[i],dat.dat[j].ambername,j,dat.n);}
                    if (strncmp(*top.atoms.atomtype[i],dat.dat[j].ambername,(strlen(*top.atoms.atomtype[i])+1)) == 0 ) {
                        radii[n] = dat.dat[j].radius;
                        if (bVerbose){fprintf(stderr,"\nMatched %s on %s to %s! (%d/%d)\n",resname,*top.atoms.atomtype[i],dat.dat[j].ambername,j,dat.n);}
                        break;
                    }
                }
                if (radii[n] == -1000.0) {
                    fprintf(stderr,"\nFrame %d : Could not match %s %s %s %i\n",framen,atomname,resname,*top.atoms.atomtype[i],n+1);
                    exit(1);
                }
            }
        }
        for (int n=0 ; n < i_index ; n++ )
        {
            radius = -1000.;
            int i = ind_index[n];
            //
            pfield += calculate_field(i,i,top,fr,midpoint,bondvector,bondlength) / pdie;
            //
            resid = top.atoms.atom[i].resind;
            charge = top.atoms.atom[i].q;
            atomname = *top.atoms.atomname[i];
            resname = *top.atoms.resinfo[resid].name;
/*            for (int j=0;j<dat.n;j++){
//                if (strncmp(resname,dat.dat[j].resname,strlen(resname)+1) == 0 ) {
                    if ( bVerbose) {fprintf(stderr,"\nTrying to match: %s %s to %s (%d/%d)",resname,*top.atoms.atomtype[i],dat.dat[j].ambername,j,dat.n);}
//                    if ( bVerbose) {fprintf(stderr,"\nTrying to match: %s %s to %s %s (%d/%d)",resname,*top.atoms.atomtype[i],dat.dat[j].resname,dat.dat[j].ambername,j,dat.n);}

                    if (strncmp(*top.atoms.atomtype[i],dat.dat[j].ambername,(strlen(*top.atoms.atomtype[i])+1)) == 0 ) {
                        radius = dat.dat[j].radius;
//                      fprintf(stderr,"\n%s =? %s -> %s", *top.atoms.atomtype[i], dat.dat[j].gmxname, dat.dat[j].ambername);
                            break;
                    }
//                }
            }
            if (radius == -1000.0) {
                fprintf(stderr,"\nFrame %d : Could not match %s %s %s %i\n",framen,atomname,resname,*top.atoms.atomtype[i],n+1);
                exit(1);
            }
 */
            fprintf(pqrout,"ATOM %6i %4s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",i+1,atomname,resname,resid+1,fr.x[i][XX]*10,fr.x[i][YY]*10,fr.x[i][ZZ]*10,charge,radii[n]);
        }
        fclose(pqrout);
        framen++;
        //
        fprintf(fout,"\n %18.6e %18.6e %18.6e ",pfield,pcfield,pnfield);
        double efield = pfield - pcfield - pnfield;
        if (a3) {
            for (int i=0;i<n_field_atoms;i++) {
                double afield = calculate_field(field_due_to_atom[i],field_due_to_atom[i],top,fr,midpoint,bondvector,bondlength)/pdie;
                fprintf(fout,"%18.6e ",afield);
                efield -= afield;
            }
        }
        fprintf(fout,"%18.6e\n",efield);

    } while(read_next_frame(oenv, status, &fr));
    fclose(fout);
    return 0;
}


gmx2amb_dat read_DAT( const char * filename, const char * atpname, gmx2amb_dat dat, gmx_bool bVerbose )
{
    char line[1024];
    int nlines = 0;
    FILE *f = fopen(filename,"r");
    while (fgets(line,1024,f))
    {
        nlines +=1 ;
    }
    gmx2amb *mygmx2amb = (gmx2amb*)malloc(sizeof(gmx2amb)*nlines*2);
    if (mygmx2amb == NULL) {
        fprintf(stderr,"\nError allocating memory in gmx2amb!");
        exit(1);
    }
    if ( f == NULL ) {
        fprintf(stderr,"\nError opening %s",filename);
        exit(1);
    }
//    printf("%i\n",nlines);

    fclose(f);
    FILE *F=fopen(filename,"r");
    int n=0,j;
    while (fgets(line,1024,f))
    {
        char escapechar[3];
        char resname[20];
        char ambername[20];
        char gmxname[20];
        float radius, charge;
        gmx_bool add_to_library = TRUE;
        sscanf(line,"%s",escapechar);
//        fprintf(stderr,"\n%s",line);
        if (strncmp(escapechar,"#",1) != 0)
        {
            sscanf(line,"%s %s %f %f %s", resname, gmxname, &charge, &radius, ambername);
            // GLY H0 in gromacs is GLY H1 in AMBER.DAT
            if ( strncmp(resname,"GLY",strlen(resname)+1) == 0 ) {
                if ( strncmp(ambername,"H1",strlen(ambername)+1) == 0 ){
//                    fprintf(stderr,"\n%s %s\n",resname,ambername);
                    memcpy(ambername,"H0",strlen("H0")+1);
                }
            }
            for (int i=0;i<n;i++) {
                if ( strncmp(ambername,mygmx2amb[i].ambername,strlen(ambername)+1) == 0 ) {
                    add_to_library = FALSE;
                }
            }
            if ( add_to_library ) {
                memcpy(mygmx2amb[n].ambername,ambername,strlen(ambername)+1);
                mygmx2amb[n].radius=radius;
                n++;
            }

//            memcpy(mygmx2amb[n].resname,resname,strlen(resname)+1);
//            memcpy(mygmx2amb[n].gmxname,gmxname,strlen(gmxname)+1);
//            memcpy(mygmx2amb[n].ambername,ambername,strlen(ambername)+1);
//            mygmx2amb[n].radius=radius;
//            mygmx2amb[n].charge=charge;
//            n++;
        }
    }
//    bVerbose=TRUE;
    if (bVerbose) {
        for (int i=0;i<n;i++){
            fprintf(stderr,"\n%i out of %i: ambername = %s, radius = %f ",i,n,mygmx2amb[i].ambername,mygmx2amb[i].radius);
        }
        fprintf(stderr,"\n%d values",n);
    }
    fclose(f);
    dat.dat=mygmx2amb;
    dat.n=n;
//    free(atpfile);
    return dat;
}

int main(int argc, const char * argv[])
{
    gmx_gmx2pqr(argc,argv);
    return 0;
}
