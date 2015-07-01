//
//  main.cpp
//  gmx2pqr
//
//  Created by Andrew Ritchie on 5/3/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

/* 
 * I want the option parser to have a flag for:
 *  -single site potentials
 *  -chromophore (CD and NE)
 * If either are unassigned, make a pqr without dummy atoms
 */


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

#include <stdio.h>
#include "my_structs.h"

gmx2amb_dat read_DAT( const char * filename, const char * atpname, gmx2amb_dat dat, gmx_bool bVerbose );

double dot( double a[3], double b[3] ) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cross( double *a, double *b, double *ab )
{
    ab[0] = a[1] * b[2] - a[2] * b[1];
    ab[1] = a[2] * b[0] - a[0] * b[2];
    ab[2] = a[0] * b[1] - a[1] * b[0];
    return;
}

void ring_points( rvec point, double *xvec, double *yvec, double *zvec, float ring_dist, double (*each_point)[3], int *total_points )
{
    int j = total_points[0];
    // plus x
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] + ring_dist * xvec[i];
    }
    j++;
    // minus x
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] - ring_dist * xvec[i];
    }
    j++;
    // plus y
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] + ring_dist * yvec[i];
    }
    j++;
    // minus y
    for (int i=0; i<3; i++){
        each_point[j][i] = point[i] - ring_dist * yvec[i];
    }
    j++;
    double pos[3];
    // normalized[ (plus x) + (plus y) ] -> j-4, j-2
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] + ring_dist * (xvec[i] + yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (plus x) + (minus y) ] -> j-4, j-1
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] + ring_dist * (xvec[i] - yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (minus x) + (plus y) ] -> j-3, j-2
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] - ring_dist * (xvec[i] + yvec[i]) / sqrt(2);
    }
    j++;
    // normalized[ (minus x) + (minus y) ] -> j-3, j-1
    for (int i=0; i<3; i++)
    {
        each_point[j][i] = point[i] - ring_dist * (xvec[i] - yvec[i]) / sqrt(2);
    }
    j++;
    
    total_points[0] = j;
}

void pm_delta_pqrline( FILE *pqrout, char *name, rvec point, double *xvec, double *yvec, double *zvec, float delta )
{
    char * resname = "DUM";
    double x,y,z;
    // Midpoint, z-dz direction
    x = point[0] - delta/2 * zvec[0];
    y = point[1] - delta/2 * zvec[1];
    z = point[2] - delta/2 * zvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"z",resname,0,x*10,y*10,z*10,0.0,0.0);
    // Midpoint, z+dz direction
    x = point[0] + delta/2 * zvec[0];
    y = point[1] + delta/2 * zvec[1];
    z = point[2] + delta/2 * zvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"z",resname,0,x*10,y*10,z*10,0.0,0.0);
    // Midpoint, y-dy direction
    x = point[0] - delta/2 * yvec[0];
    y = point[1] - delta/2 * yvec[1];
    z = point[2] - delta/2 * yvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"y",resname,0,x*10,y*10,z*10,0.0,0.0);
    // Midpoint, y+dy direction
    x = point[0] + delta/2 * yvec[0];
    y = point[1] + delta/2 * yvec[1];
    z = point[2] + delta/2 * yvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"y",resname,0,x*10,y*10,z*10,0.0,0.0);
    // Midpoint, x-dx direction
    x = point[0] - delta/2 * xvec[0];
    y = point[1] - delta/2 * xvec[1];
    z = point[2] - delta/2 * xvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"x",resname,0,x*10,y*10,z*10,0.0,0.0);
    // Midpoint, x+dx direction
    x = point[0] + delta/2 * xvec[0];
    y = point[1] + delta/2 * xvec[1];
    z = point[2] + delta/2 * xvec[2];
    fprintf(pqrout,"ATOM %6i %3s%s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,name,"x",resname,0,x*10,y*10,z*10,0.0,0.0);
    return;
}

double calculate_atom_r( int i, int j, t_topology top, t_trxframe fr)
{
    double xyz[3];
    double r2, r;
    for (int k=0;k<3;k++)
    {
        xyz[k] = (fr.x[j][k] - fr.x[i][k]);
    }
    r2 = dot(xyz,xyz);
    return pow(r2,0.5);
}

double calculate_r( double *point, int i, t_topology top, t_trxframe fr)
{
    double xyz[3];
    double r2, r;
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    return pow(r2,0.5);
}

double calculate_potential( double r, int i, t_topology top, t_trxframe fr)
{
    double rr;
    rr = 1./r;
    return rr * top.atoms.atom[i].q * cfac;
}

/*double calculate_potential( double *point, int i, t_topology top, t_trxframe fr)
{
    double xyz[3];
    double r2, rr;
    // don't forget to convert from nm to Angstroms
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    rr = pow(r2,-0.5);
    return rr * top.atoms.atom[i].q * cfac;
}*/

void calculate_field( rvec point, int i, t_topology top, t_trxframe fr, double *field)
{
    double xyz[3];
    double r, r2, rr3;
    // don't forget to convert from nm to Angstroms
    for (int j=0;j<3;j++)
    {
        xyz[j] = 10 * (point[j] - fr.x[i][j]);
    }
    r2 = dot(xyz,xyz);
    r = pow(r2, 0.5);
    rr3 = 1.0 / (r * r * r);
    for (int j=0; j<3; j++)
    {
        field[j] += rr3 * top.atoms.atom[i].q * xyz[j] * cfac;
    }
    return;
}

double project_field(rvec vector, double *field)
{
    double projection = 0;
    double vlen = pow(iprod(vector,vector),.5);
    for (int i=0; i<3; i++)
    {
        projection += field[i] * vector[i] / vlen;
    }
    return projection;
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
    float           pdie=1,delta=0.005,ring_dist=0.07,cutoff=0.5;
    double          pfield,pcfield,pnfield,psfield;
    gmx_bool        bVerbose=FALSE;     // Remember, 0=False, 1=True
    static char     *exclude=NULL, *site=NULL;
    
    t_pargs pa[] = {
        { "-nn", TRUE, etBOOL,
            {&nn}, "Use a separate index file for each frame, IE: dynamic index files created through g_select."},
        { "-a1", TRUE, etINT,
            {&a1}, "Atom number 1 (CNC CD, must be set!)"},
        { "-a2", TRUE, etINT,
            {&a2}, "Atom number 2 (CNC NE, must be set!)"},
        { "-site", TRUE, etSTR,
            {&site},"String, in quotations, for atoms you want to include as sites for the Cho multivariate treatment.  Do not need to include atoms specified with -a1 and -a2"},
        { "-exclude", TRUE, etSTR,
            {&exclude},"String, in quotations, for atoms you want to exclude from the field calculation"},
        { "-pdie", TRUE, etREAL,
            {&pdie}, "Protein dielectric"},
        { "-delta", TRUE, etREAL,
            {&delta}, "deltaR (nm?) for dummy atom positioning"},
        { "-cutoff", TRUE, etREAL,
            {&cutoff}, "Only include water molecules within -cutoff (nm) of at least one site"},
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
    if (exclude) {
        for (word = strtok_r(exclude, sep, &brkt); word; word = strtok_r(NULL, sep, &brkt))
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
    int site_atom[1000], n_sites=0;
    if (site) {
        for (word = strtok_r(site, sep, &brkt); word; word = strtok_r(NULL, sep, &brkt))
        {
            int field_atom = atoi(word);
            if (field_atom>(top.atoms.nr))
            {
                fprintf(stderr,"Error: Atom number %d is out of range.\n",field_atom);
                exit(1);
            }
            site_atom[n_sites] = field_atom - 1;
            n_sites++;
        }
        site_atom[n_sites] = a1;
        n_sites++;
        site_atom[n_sites] = a2;
        n_sites++;
        printf("\nWill make potential sites at atoms: ");
        for (int i=0;i<n_sites;i++){printf("%i ",site_atom[i]+1);}
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
    
    // Give the xvg file a header to help it make sense without reading the code
    fprintf(fout,";%13i-%i %15i%-3s %15i%-3s ",ind_index[0]+1,ind_index[i_index-1]+1,a1+1,*top.atoms.atomname[a1],a2+1,*top.atoms.atomname[a2]);
    if (exclude) {
        for (int i=0;i<n_field_atoms;i++){
            fprintf(fout,"%15i%-3s ",field_due_to_atom[i]+1,*top.atoms.atomname[field_due_to_atom[i]]);
        }
    }
    fprintf(fout,"%18s ","F-excpt");
    fprintf(fout,"%18s %18s %18s ","MIDX","MIDY","MIDZ");
    fprintf(fout,"%17s%s %17s%s %17s%s ",*top.atoms.atomname[a1],"X",*top.atoms.atomname[a1],"Y",*top.atoms.atomname[a1],"Z");
    fprintf(fout,"%17s%s %17s%s %17s%s ",*top.atoms.atomname[a2],"X",*top.atoms.atomname[a2],"Y",*top.atoms.atomname[a2],"Z");
    if (site){
        for (int n=0; n<n_sites; n++){
            fprintf(fout,"%16s%-2i %16s%-2i %16s%-2i ","coulomb_protein_p",n,"coulomb_solvent_p",n,"coulomb_total_p",n);}
            //fprintf(fout,"%16s%-2i %16s%-2i %16s%-2i ","-d coulomb_protein_p -u kbt/e",n,"-d coulomb_solvent_p -u kbt/e",n,"-d coulomb_total_p -u kbt/e",n);}
        for (int n=0; n<17; n++){ // 8 per atom (2*16) plus 1 out along the bond vector
            fprintf(fout,"%16s%-2i %16s%-2i ","coulomb_protein_p",n_sites+n,"coulomb_solvent_p",n_sites+n);
            //fprintf(fout,"%16s%-2i %16s%-2i %16s%-2i ","-d coulomb_protein_p -u kbt/e",n_sites+n,"-d coulomb_solvent_p -u kbt/e",n_sites+n,"-d coulomb_total_p -u kbt/e",n_sites+n);
        }
    }

    do {
        if (nn && framen > 0) {
            get_index(atoms,index_file,1,&i_index,&ind_index,&gn_index);
        }
        
        sprintf(pqr,"%s%i.pqr",pqr_file,framen);
        FILE *pqrout = ffopen(pqr,"w");
        if ( pqrout == NULL ) {
            gmx_fatal(FARGS,"\nFailed to open output file, '%s'\n",pqr);
        }
        
        // Calculate the bond vector and midpoint
        rvec_sub(fr.x[a2],fr.x[a1],bondvector); // Nitrile bond vector we're looking at
        bondlength = pow(iprod(bondvector,bondvector),.5);
        for (int i=0; i<3; i++){
            midpoint[i] = fr.x[a1][i] + 0.5 * bondvector[i];
        }
        
        /* 
         *  Find the points we need to make dummy atoms :
         *      I will consider the bond vector the z axis.  I need to
         *      find an x and y axis perpendicular.
         */
        double zvec[3] = { bondvector[0] / bondlength, bondvector[1] / bondlength, bondvector[2] / bondlength };
        double yvec[3] = { 1,zvec[2],-zvec[1] };
        yvec[0] = -1*(zvec[1]*yvec[1]+zvec[2]*yvec[2])/zvec[0];
        double ylen = pow(dot(yvec,yvec),0.5);
        for (int i=0;i<3;i++){
            yvec[i] = yvec[i] / ylen;
        }
        double xvec[3];
        cross(yvec,zvec,xvec);

        // Calculate the field at the midpoint
        pm_delta_pqrline(pqrout, "MID", midpoint, xvec, yvec, zvec, delta);
        double midpoint_field[3] = { 0., 0., 0. };
        double a1point_field[3] = { 0., 0., 0. };
        double a2point_field[3] = { 0., 0., 0. };
        for (int n=0; n<i_index; n++) // Loop over atoms in index
        {
            int i = ind_index[n]; //
            calculate_field(midpoint, i, top, fr, midpoint_field);
            if ( i != a1 && i != a2 ) // Not looking at the chromophore
            {
                gmx_bool isIncluded = 1;
                for (int j=0; j<n_field_atoms; j++)
                {
                    if (i == field_due_to_atom[j]) // Not looking at exception atoms
                    {
                        isIncluded = 0;
                        break;
                    }
                }
                if (isIncluded)
                {
                    calculate_field(fr.x[a1],i,top,fr,a1point_field);
                    calculate_field(fr.x[a2],i,top,fr,a2point_field);
                }
            }
        }
        double proj_midpoint_field = project_field(bondvector,midpoint_field);
        
        // Calculate the field due to a1
        pm_delta_pqrline(pqrout, *top.atoms.atomname[a1], fr.x[a1], xvec, yvec, zvec, delta);
        double a1_field[3] = { 0., 0., 0. };
        calculate_field(midpoint, a1, top, fr, a1_field);
        double proj_a1_field = project_field(bondvector,a1_field);
        
        // Calculate the field due to a2
        pm_delta_pqrline(pqrout, *top.atoms.atomname[a2], fr.x[a2], xvec, yvec, zvec, delta);
        double a2_field[3] = { 0., 0., 0. };
        calculate_field(midpoint, a2, top, fr, a2_field);
        double proj_a2_field = project_field(bondvector,a2_field);
        
        // Write to xvg file
        fprintf(fout,"\n %18.6e %18.6e %18.6e ",proj_midpoint_field/pdie,proj_a1_field/pdie,proj_a2_field/pdie);
        fprintf(stderr,"\n%18.6f \n%18.6f \n%18.6f\n",proj_midpoint_field/pdie,proj_a1_field/pdie,proj_a2_field/pdie);


        double except_field = proj_midpoint_field - proj_a1_field - proj_a2_field;
        // Calculate the field due to each except[i] and write to xvg file
        if (exclude) {
            for (int i=0;i<n_field_atoms;i++) {
                double exclude_field[3] = { 0., 0., 0. };
                calculate_field(midpoint, field_due_to_atom[i], top, fr, exclude_field);
                double proj_exclude_field = project_field(bondvector,exclude_field);
                fprintf(fout,"%18.6e ",proj_exclude_field/pdie);
                fprintf(stderr,"%18.6f\n",proj_exclude_field/pdie);
                except_field -= proj_exclude_field;
                for (int j=0; j<3 ; j++)
                {
                    midpoint_field[j] -= exclude_field[j];
                }
            }
        }
        fprintf(fout,"%18.6e ",except_field/pdie);
        fprintf(stderr,"%18.6f\n",except_field/pdie);
        rvec rxvec, ryvec, rzvec;
        for (int i=0; i<3 ; i++)
        {
            midpoint_field[i] -= a1_field[i] + a2_field[i];
            rxvec[i] = xvec[i];
            ryvec[i] = yvec[i];
            rzvec[i] = zvec[i];
        }
        fprintf(fout,"%18.6e %18.6e %18.6e ",project_field(rxvec,midpoint_field)/pdie,project_field(ryvec,midpoint_field)/pdie,project_field(rzvec, midpoint_field)/pdie);
        fprintf(fout,"%18.6e %18.6e %18.6e ",project_field(rxvec, a1point_field)/pdie,project_field(ryvec, a1point_field)/pdie,project_field(rzvec,  a1point_field)/pdie);
        fprintf(fout,"%18.6e %18.6e %18.6e ",project_field(rxvec, a2point_field)/pdie,project_field(ryvec, a2point_field)/pdie,project_field(rzvec,  a2point_field)/pdie);
        fprintf(stderr,"%.6f %.6f %.6f\n",project_field(rxvec,midpoint_field)/pdie,project_field(ryvec,midpoint_field)/pdie,project_field(rzvec, midpoint_field)/pdie);
        fprintf(stderr,"%.6f %.6f %.6f\n",project_field(rxvec, a1point_field)/pdie,project_field(ryvec, a1point_field)/pdie,project_field(rzvec,  a1point_field)/pdie);
        fprintf(stderr,"%.6f %.6f %.6f\n",project_field(rxvec, a2point_field)/pdie,project_field(ryvec, a2point_field)/pdie,project_field(rzvec,  a2point_field)/pdie);
        
        /* 
         *  Now we generate the potentials at all the sites...  I could probably
         *  do some loop fusing here, but I'm not going to because it's already
         *  1) fast and 2) easier to read maybe?
         */
        int nwaters = 0;
        int * keep_waters;
        if ( NULL == (keep_waters = malloc(top.atoms.nr * sizeof(int))) ) {
            printf("malloc failed\n");
            return(-1);
        }
        if (site)
        {
            int total_points = 0;
            double each_point[100][3]; // 1st dimension is arbitrary in size b/c I'm lazy
            for (int n=0; n<n_sites; n++)
            {
                int i = site_atom[n];
                for (int j = 0; j<3; j++)
                {
                    each_point[total_points][j] = fr.x[i][j];
                }
                total_points++;
            }
            // Point between nitrile and hydrogen-water along bond vector
            for (int i=0; i<3; i++)
            {
                each_point[total_points][i] =  fr.x[a2][i] + ring_dist * bondvector[i] / bondlength;
            }
            total_points++;
            // 8 Points around C
            ring_points( fr.x[a1], xvec, yvec, zvec, ring_dist, &each_point[0], &total_points );
            
            // 8 points are N
            ring_points( fr.x[a2], xvec, yvec, zvec, ring_dist, &each_point[0], &total_points );
            
            for (int n=0; n<total_points; n++)
            {
                double all_potential = 0, protein_potential = 0, solvent_potential = 0;
                for (int j=0; j<i_index; j++)
                {
                    int k = ind_index[j];
                    gmx_bool isIncluded = 1;
                    for (int l=0; l<n_sites; l++)
                    {
                        int m = site_atom[l];
                        if ( k == m )
                        {
                            isIncluded = 0;
                            break;
                        }
                    }
                    for (int l=0; l<n_field_atoms; l++)
                    {
                        int m = field_due_to_atom[l];
                        if ( k == m )
                        {
                            isIncluded = 0;
                            break;
                        }
                    }
                    if (isIncluded)
                    {
                        double r = calculate_r(each_point[n],k,top,fr);
                        protein_potential += calculate_potential(r,k,top,fr);
                    }
                }
                //for (int jj=0; jj<i_index; jj++)
                if (cutoff > 0){
                    for (int j=0; j<top.atoms.nr ; j++)
                    {
                        //int j = ind_index[jj];
                        if (strncmp("OW",*top.atoms.atomname[j],strlen(*top.atoms.atomname[j])+1) == 0 )
                        {
                            gmx_bool keep = 0;
                            for (int m=0; m<n_sites; m++)
                            {
                                if (keep == 0)
                                {
                                    int k = site_atom[m];
                                    for (int o=0; o<3; o++)
                                    {
                                        if (calculate_atom_r(j+o,k,top,fr) <= cutoff)
                                        {
                                            //double x = calculate_atom_r(j+o,k,top,fr);
                                            //printf("%i %s, %i %s, %f\n",j,*top.atoms.atomname[j+o],k,*top.atoms.atomname[k],x);
                                            keep = 1;
                                            break;
                                        }
                                    }
                                }
                            }
                            if (keep)
                            {
                                solvent_potential += calculate_potential(calculate_r(each_point[n],j,top,fr),j,top,fr);
                                if (n==0){ // So we know what to keep in the pqr
                                    keep_waters[nwaters] = j;
                                    keep_waters[nwaters+1] = j+1;
                                    keep_waters[nwaters+2] = j+2;
                                    nwaters += 3;
                                }
                            }
                        }
                    }
                }
                fprintf(pqrout,"ATOM %6i %2s%-2i %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",0,"p",n,"DUM",0,each_point[n][0]*10,each_point[n][1]*10,each_point[n][2]*10,0.0,0.0);
                fprintf(fout,"%18.6e %18.6e ",protein_potential,solvent_potential);
                fprintf(stderr,"%18.6f %18.6f \n",protein_potential,solvent_potential);

                //printf("Keeping %i waters\n",nwaters);
            }
        }
        
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
        for (int n=0; n<nwaters;n++)
        {
            radius = -1000.;
            int i = keep_waters[n];
            
            resid = top.atoms.atom[i].resind;
            charge = top.atoms.atom[i].q;
            atomname = *top.atoms.atomname[i];
            resname = *top.atoms.resinfo[resid].name;
            for ( int j=0; j < dat.n ; j++ ) {
                if ( bVerbose) {fprintf(stderr,"\nTrying to match: %s %s to %s (%d/%d)",resname,*top.atoms.atomtype[i],dat.dat[j].ambername,j,dat.n);}
                if (strncmp(*top.atoms.atomtype[i],dat.dat[j].ambername,(strlen(*top.atoms.atomtype[i])+1)) == 0 ) {
                    radius = dat.dat[j].radius;
                    if (bVerbose){fprintf(stderr,"\nMatched %s on %s to %s! (%d/%d)\n",resname,*top.atoms.atomtype[i],dat.dat[j].ambername,j,dat.n);}
                    break;
                }
            }
            if (radius == -1000.0) {
                fprintf(stderr,"\nFrame %d : Could not match %s %s %s %i\n",framen,atomname,resname,*top.atoms.atomtype[i],n+1);
                exit(1);
            }
            fprintf(pqrout,"ATOM %6i %4s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",i+1,atomname,resname,resid+1,fr.x[i][XX]*10,fr.x[i][YY]*10,fr.x[i][ZZ]*10,charge,radius);
        }
        fclose(pqrout);
        free(keep_waters);
        framen++;
        //
    } while(read_next_frame(oenv, status, &fr));
    fprintf(fout,"\n");
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
