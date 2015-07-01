/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
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
 * For more info, check our website at http://www.gromacs.org
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


typedef struct {
    int bond[5];
    int nb;
} t_bonded;


t_bonded * read_topol( t_bonded *bonds, const char *top_file )
{
    char line[1024];
    gmx_bool has_bonding = FALSE;
    gmx_bool in_bonds = FALSE;

    
    FILE *f = fopen(top_file, "r");
    while (fgets(line,1024,f))
    {
        char char1[25];
        char char2[25];
        char junk[25];
        memset(char1,0,sizeof(char1));
        memset(char2,0,sizeof(char2));
        sscanf(line,"%s %s",char1,char2);
        if (in_bonds) {
            int atom1,atom2;
            atom1=-1;
            if (strncmp(char1,";",1) != 0 && strncmp(char1,"[",1) != 0)
            {
                atom1 = atoi(char1);
                atom2 = atoi(char2);
                bonds[atom1].bond[bonds[atom1].nb]=atom2;
                bonds[atom2].bond[bonds[atom2].nb]=atom1;
                bonds[atom1].nb++;
                bonds[atom2].nb++;
            }
        }
        if (strncmp(char1,"[",1) == 0) {
            if (strncmp(char2,"bond",4) == 0) {
                has_bonding = TRUE;
                in_bonds = TRUE;
            }
            else {in_bonds = FALSE;}
        }
    }
    if ( f == NULL ) {
        fprintf(stderr,"\nError opening %s",top_file);
        exit(1);
    }
    if ( !has_bonding) {
        fprintf(stderr,"\nError!  There are no bonding parameters in your topology file!  You may need to catenate your itp file into the topology.\n");
        exit(1);
    }
    
    fclose(f);
    return bonds;
}

typedef struct {
    char atomname[5];
    char resname[5];
    int  atomid;
} amoeba_parm;

typedef struct {
    amoeba_parm * parm;
    int nparm;
} amoeba_parms;

amoeba_parm ambernames( amoeba_parm para, char c1[25], char c2[25], char c3[25], char c4[25], char c5[25], char atomname[5])
{
    char resname[5];
    char atmname[5];
    int atomid=-1;
    memset(resname,0,sizeof(resname));
    
    if (strncmp(atomname,"HN",3) == 0) {
        strcpy(atmname,"H");
        memcpy(atomname,atmname,strlen(atmname)+1);
    }
    
    if (strncmp(c2,"\"",3) == 0) {
        if      (strncmp(c1,"\"Glycine",8)      == 0) {strcpy(resname,"GLY");}
        else if (strncmp(c1,"\"Alanine",8)      == 0) {strcpy(resname,"ALA");}
        else if (strncmp(c1,"\"Valine",7)       == 0) {strcpy(resname,"VAL");}
        else if (strncmp(c1,"\"Leucine",8)      == 0) {strcpy(resname,"LEU");}
        else if (strncmp(c1,"\"Isoleucine",8)   == 0) {strcpy(resname,"ILE");}
        else if (strncmp(c1,"\"Serine",7)       == 0) {strcpy(resname,"SER");}
        else if (strncmp(c1,"\"Threonine",8)    == 0) {strcpy(resname,"THR");}
        else if (strncmp(c1,"\"Proline",7)      == 0) {strcpy(resname,"PRO");}
        else if (strncmp(c1,"\"Phenylalanin",8) == 0) {strcpy(resname,"PHE");}
        else if (strncmp(c1,"\"Tyrosine",8)     == 0) {strcpy(resname,"TYR");}
        else if (strncmp(c1,"\"Tryptophan",8)   == 0) {strcpy(resname,"TRP");}
        else if (strncmp(c1,"\"Asparagine",8)   == 0) {strcpy(resname,"ASN");}
        else if (strncmp(c1,"\"Glutamine",8)    == 0) {strcpy(resname,"GLN");}
        else if (strncmp(c1,"\"Methionine",8)   == 0) {strcpy(resname,"MET");}
        else if (strncmp(c1,"\"Lysine",7)       == 0) {strcpy(resname,"LYP");}
        else if (strncmp(c1,"\"Arginine",8)     == 0) {strcpy(resname,"ARG");}
        else if (strncmp(c1,"\"Water",8)        == 0) {
            strcpy(resname,"SOL");
            if (strncmp(atomname,"O",1) == 0) {
                strcpy(atmname,"OW");
            }
            else {
                strcpy(atmname,"HW"); //1 or 2
            }
            memcpy(atomname,atmname,strlen(atmname)+1);
        }
        else if (strncmp(c1,"\"CNC",4)          == 0) {strcpy(resname,"CNC");}
        else if (strncmp(c1,"\"DCN",4)          == 0) {strcpy(resname,"DCN");}
        else if (strncmp(c1,"\"GNP",4)          == 0) {
            strcpy(resname,"GNP");
            if (strncmp(&atomname[3],"'",1) == 0){
                strcpy(atmname,&atomname[1]);
                strncat(atmname,&atomname[0],1);
                memcpy(atomname,atmname,strlen(atmname)+1);
            }
        }
        atomid = atoi(c3);    
    }
    else if (strncmp(c2,"(SH)",4) == 0) {
        strcpy(resname,"CYN");
        atomid = atoi(c4);
    }
    else if (strncmp(c2,"(SS)",4) == 0) {
        strcpy(resname,"CYM");
        atomid = atoi(c4);
    }
    else if (strncmp(c2,"(+)",4) == 0) {
        strcpy(resname,"HIP");
        atomid = atoi(c4);
    }
    else if (strncmp(c2,"(HD)",4) == 0) {
        strcpy(resname,"HID");
        atomid = atoi(c4);
    }
    else if (strncmp(c2,"(HE)",4) == 0) {
        strcpy(resname,"HIE");
        atomid = atoi(c4);
    }
    else if (strncmp(c2,"Acid",4) == 0) {
        atomid = atoi(c4);
        if (strncmp(c1,"\"Glutamic",9) == 0) {
            strcpy(resname,"GLU");
        }
        else if (strncmp(c1,"\"Aspartic",9) == 0) {
            strcpy(resname,"ASP");
        }
    }
    else if (strncmp(c2,"Ion",3) == 0) {
        atomid = atoi(c4);
        if (strncmp(c1,"\"Sodium",7) == 0) {
            strcpy(resname,"Na+");
            strcpy(atmname,"Na");
            memcpy(atomname,atmname,strlen(atmname)+1);
        }
        else if (strncmp(c1,"\"Potassium",7) == 0) {
            strcpy(resname,"K");
        }
        else if (strncmp(c1,"\"Magnesium",7) == 0) {
            strcpy(resname,"MG");
        }
        else if (strncmp(c1,"\"Calcium",7) == 0) {
            strcpy(resname,"Ca");
        }
        else if (strncmp(c1,"\"Chloride",7) == 0) {
            strcpy(resname,"Cl");
        }
        else if (strncmp(c1,"\"Zinc",5) == 0) {
            strcpy(resname,"Zn");
        }
    }
    else if (strncmp(c1,"\"N-Terminal",10) == 0) {
        strcpy(resname,"N");
        if (strncmp(c3,"\"",1) ==0) {
            strcat(resname,c2);
            atomid = atoi(c4);
        }
        else if (strncmp(c3,"(SH)",4) == 0) {
            strcat(resname,"CYN");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(SS)",4) == 0) {
            strcat(resname,"CYM");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(+)",4) == 0) {
            strcat(resname,"HIP");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(HD)",4) == 0) {
            strcat(resname,"HID");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(HE)",4) == 0) {
            strcat(resname,"HIE");
            atomid = atoi(c5);
        }
    }
    else if (strncmp(c1,"\"C-Terminal",10) ==0) {
        strcpy(resname,"C");
        if (strncmp(atomname,"OXT",4) == 0) {
            strcpy(atmname,"OC"); //1 or 2
            memcpy(atomname,atmname,strlen(atmname)+1);
        }
        if (strncmp(c3,"\"",1) ==0) {
            strcat(resname,c2);
            atomid = atoi(c4);}
        else if (strncmp(c3,"(SH)",4) == 0) {
            strcat(resname,"CYN");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(SS)",4) == 0) {
            strcat(resname,"CYM");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(+)",4) == 0) {
            strcat(resname,"HIP");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(HD)",4) == 0) {
            strcat(resname,"HID");
            atomid = atoi(c5);
        }
        else if (strncmp(c3,"(HE)",4) == 0) {
            strcat(resname,"HIE");
            atomid = atoi(c5);
        }
    }
    para.atomid=atomid;
    memcpy(para.resname,resname,strlen(resname)+1);
    memcpy(para.atomname,atomname,strlen(atomname)+1);
//    printf("%s %i : %s %s %s %s %s\n",resname,atomid,c1,c2,c3,c4,c5);
//    printf("%s %s %i\n",para.resname,para.atomname,para.atomid);
    return para;
}

amoeba_parm * read_amoeba_parameters( const char *parm_name )
{
    char line[1024];
    int n=0;
    FILE *f = fopen(parm_name,"r");
    amoeba_parm *param = (amoeba_parm*)malloc(sizeof(amoeba_parm)*0);
    while (fgets(line,1024,f)){
        char b[25],buffer[25];
        char atomname[5];
        int atomid;
        char char4[25],char5[25],char6[25],char7[25],char8[25];
        memset(b,0,sizeof(b));
        memset(atomname,0,sizeof(atomname));
        sscanf(line, "%s %s %s %s %s %s %s %s",b,buffer,atomname,char4,char5,char6,char7,char8);
        if (strncmp(b,"biotype",3) == 0 ) {
            param = realloc(param,sizeof(amoeba_parm)*n+1);
            param[n]=ambernames(param[n],char4,char5,char6,char7,char8,atomname);
//            printf("%s",line);
//            printf("outside: %s %s %i\n\n",param[n].resname,param[n].atomname,param[n].atomid);
            n++;
        }
    }
//    for (int i=0;i<n;i++){printf("%i %s %s %i\n",i,param[i].resname,param[i].atomname,param[i].atomid);}
    return param;
}

int * parm_order(int * order, t_topology top, amoeba_parm * parms, int nparams, int i_index, atom_id *index )
{
    int resid;
    int atomid;
    int matchlen;
    int i;
    char atomname[5];
    char resname[5];
    for (int n=0;n<i_index;n++) {
        i=index[n];
        resid = top.atoms.atom[i].resind;
        atomid=0;
        strcpy(resname,*top.atoms.resinfo[resid].name);
        strcpy(atomname,*top.atoms.atomname[i]);
        
        /*// I realized that since it's usually 1-2 characters shorter,
        // I can just successively search 1-2 characters shorter until I
        // get a match!
        if (strncmp(resname,"GLY",strlen(resname)+1) == 0) {
            if (strncmp(atomname,"HA1",strlen(atomname)+1) == 0 
                || strncmp(atomname,"HA2",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HA");                
            }                    
        }
        else if (strncmp(resname,"ALA",strlen(resname)+1) == 0) {
            if (strncmp(atomname,"HB1",strlen(atomname)+1) == 0 
                || strncmp(atomname,"HB2",strlen(atomname)+1) == 0 
                || strncmp(atomname,"HB3",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HB");
            }
        }
        else if (strncmp(resname,"VAL",strlen(resname)+1) == 0) {
            if (strncmp(atomname,"HG11",strlen(atomname)+1) == 0 
                || strncmp(atomname,"HG12",strlen(atomname)+1) == 0
                || strncmp(atomname,"HG13",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HG1");
            }
            if (strncmp(atomname,"HG21",strlen(atomname)+1) == 0
                || strncmp(atomname,"HG22",strlen(atomname)+1) == 0
                || strncmp(atomname,"HG23",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HG2");
            }
        }
        else if (strncmp(resname,"LEU",strlen(resname)+1) == 0) {
            if (strncmp(atomname,"HB1",strlen(atomname)+1) == 0
                || strncmp(atomname,"HB2",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HB");
            }
            else if (strncmp(atomname,"HD11",strlen(atomname)+1) == 0
                || strncmp(atomname,"HD12",strlen(atomname)+1) == 0
                || strncmp(atomname,"HD13",strlen(atomname)+1) == 0
                ) {
                strcpy(atomname,"HD1");
            }
            else if (strncmp(atomname,"HD21",strlen(atomname)+1) == 0
                     || strncmp(atomname,"HD22",strlen(atomname)+1) == 0
                     || strncmp(atomname,"HD23",strlen(atomname)+1) == 0
                     ) {
                strcpy(atomname,"HD2");
            }
        }
        */
        matchlen=strlen(atomname)+1;
        for (int ii=0 ; ii<4 ; ii++) {
            //printf("\n%5i %6i%s %6s(%i) ",i,resid,resname,atomname,matchlen-ii);
            for (int k=0;k<nparams;k++) {
                if (strncmp(resname, parms[k].resname, strlen(resname)+1) == 0 ) {
                    if (strncmp(atomname,parms[k].atomname,matchlen-ii) == 0) {
                        /*
                        printf("<%s ",parms[k].resname);
                        printf("%s ",parms[k].atomname);
                        printf("%5i> ",parms[k].atomid);
                        */
                        atomid=parms[k].atomid;
                        order[i]=atomid;
                        if (atomid != 0) { break; }
                    }
                    if (atomid != 0) { break; }
                }
                if (atomid != 0) { break; }
            }
            if (atomid != 0) { break; }
        }
        if ((strncmp(resname,"C",1) == 0 && strlen(resname) > 3) &&
            (
             strncmp(atomname,"OC",3) != 0 &&
             strncmp(atomname,"OC1",4) != 0 &&
             strncmp(atomname,"OC2",4) != 0 &&
             strncmp(atomname,"H",2) != 0 &&
             strncmp(atomname,"C",2) != 0 &&
             strncmp(atomname,"CA",3) != 0 &&
             strncmp(atomname,"N",2) != 0 &&
             strncmp(atomname,"HA",3) != 0 &&
             strncmp(atomname,"HA1",4) != 0 &&
             strncmp(atomname,"HA2",4) != 0 &&
             strncmp(atomname,"HA3",4) != 0
             ))
        {
            atomid = 0;
            strcpy(resname,&resname[1]);
            for (int ii=0 ; ii<matchlen-1 ; ii++) {
                //                    printf("\n%5i %6i%s %6s(%i) ",i,resid,resname,atomname,matchlen-ii);
                for (int k=0;k<nparams;k++) {
                    if (strncmp(resname, parms[k].resname, strlen(resname)+1) == 0 ) {
                        if (strncmp(atomname,parms[k].atomname,matchlen-ii) == 0) {
                            /*
                             printf("<%s ",parms[k].resname);
                             printf("%s ",parms[k].atomname);
                             printf("%5i> ",parms[k].atomid);
                             */
                            atomid=parms[k].atomid;
                            order[i]=atomid;
                            if (atomid != 0) { break; }
                        }
                        if (atomid != 0) { break; }
                    }
                    if (atomid != 0) { break; }
                }
                if (atomid != 0) { break; }
            }
        }
        else if ((strncmp(resname,"N",1) == 0 && strlen(resname) > 3) &&
                 (
                  strncmp(atomname,"O",2) != 0 &&
                  strncmp(atomname,"H",2) != 0 &&
                  strncmp(atomname,"H1",3) != 0 &&
                  strncmp(atomname,"H2",3) != 0 &&
                  strncmp(atomname,"H3",3) != 0 &&
                  strncmp(atomname,"C",2) != 0 &&
                  strncmp(atomname,"N",2) != 0 &&
                  strncmp(atomname,"HA",3) != 0 &&
                  strncmp(atomname,"HA1",4) != 0 &&
                  strncmp(atomname,"HA2",4) != 0 &&
                  strncmp(atomname,"HA3",4) != 0
                  ))
        {
            atomid = 0;
            strcpy(resname,&resname[1]);
            for (int ii=0 ; ii<matchlen-1 ; ii++) {
                //            printf("\n%5i %6i%s %6s(%i) ",i,resid,resname,atomname,matchlen-ii);
                for (int k=0;k<nparams;k++) {
                    if (strncmp(resname, parms[k].resname, strlen(resname)+1) == 0 ) {
                        if (strncmp(atomname,parms[k].atomname,matchlen-ii) == 0) {
                            /*
                             printf("<%s ",parms[k].resname);
                             printf("%s ",parms[k].atomname);
                             printf("%5i> ",parms[k].atomid);
                             */
                            atomid=parms[k].atomid;
                            order[i]=atomid;
                            if (atomid != 0) { break; }
                        }
                        if (atomid != 0) { break; }
                    }
                    if (atomid != 0) { break; }
                }
                if (atomid != 0) { break; }
            }
        }
        if (atomid ==0) {
            printf("\nMissing parameters");
            printf("\n%5i %6i%s %6s",i,resid,resname,atomname);
            exit(1);
        }
    }
    return order;
}

int gmx_gmx2xyz(int argc, char *argv[])
{
    const char *desc[] = {
        "\tThis convert from gromacs to Tinker .xyx.",
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
    int             i_index;
    atom_id         *ind_index;
    char            *gn_index;
    const char      *parm_name;
    output_env_t    oenv;
    const char      *tpr_file, *top_file, *in_file, *index_file, *out_file = NULL;
    char            tpr_title[256], xyz[256];
    t_bonded        *bonds;
    
    t_pargs pa [] = {
    };
    
    t_filenm fnm[] = {
        { efTPS, NULL, NULL, ffREAD },
        { efTRX, "-f", NULL, ffREAD },
        { efTOP, NULL, NULL, ffREAD },
        { efRND, "-a", "amoeba.prm", ffREAD },
        { efXYZ, "-o", NULL, ffWRITE },
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
    top_file = ftp2fn(efTOP, NFILE, fnm);
    init_top(&top);
    in_file = opt2fn("-f", NFILE, fnm);
    parm_name = opt2fn("-a", NFILE, fnm);
    out_file = opt2fn("-o", NFILE,fnm);

    read_tps_conf(tpr_file, tpr_title, &top, &ePBC, &xtop, NULL, box, TRUE);
    atoms = &top.atoms;
    
    printf("Select group for to save atom coordinates:\n");
    index_file = ftp2fn(efNDX, NFILE, fnm);
    get_index(atoms, index_file, 1, &i_index, &ind_index, &gn_index);
    
    snew(bonds,top.atoms.nr+1);
    bonds=read_topol(bonds,top_file);
    for (int i=0;i<top.atoms.nr;i++){
        if (strncmp(*top.atoms.atomname[i],"OW",2) == 0) {
            bonds[i+1].nb = 2;
            bonds[i+1].bond[0] = i+1;
            bonds[i+1].bond[1] = i+2;
        }
        else if (strncmp(*top.atoms.atomname[i],"HW1",3) ==0) {
            bonds[i+1].nb = 1;
            bonds[i+1].bond[0] = i-1;
        }
        else if (strncmp(*top.atoms.atomname[i],"HW2",3) ==0) {
            bonds[i+1].nb = 1;
            bonds[i+1].bond[0] = i-2;
        }
    }
    /*
    for (int i=i_index-10;i<=i_index+17;i++) {
        printf("\n%i <%i>",i,bonds[i+1].nb);
        for (int j=0;j<bonds[i+1].nb;j++) {
            printf(" %i",bonds[i+1].bond[j]);
        }
    }
     */
    int nparams = 967; /* I know this from playing around.  There's gotta be a smarter way? */
    amoeba_parm * atoms_types = (amoeba_parm*)malloc(sizeof(amoeba_parm)*nparams);
    atoms_types = read_amoeba_parameters( parm_name );
    
    int * parmid = (int*)malloc(sizeof(int)*top.atoms.nr);
    parm_order(parmid, top, atoms_types, nparams,i_index,ind_index);
     
    read_first_frame(oenv, &status, in_file, &fr, flags);
    int resid;
//    strncpy(xyz_file,out_file,strlen(out_file)-4);
    int framen=0;
    char *xyz_file = out_file;
    xyz_file[strlen(xyz_file)-4] = 0;
    do {
        memset(xyz,0,sizeof(xyz));
        memcpy(xyz,xyz_file,strlen(xyz_file));//-1);
        sprintf(xyz,"%s%i.xyz",xyz,framen);
        FILE * xyzout = ffopen(xyz,"w");
        if (xyzout==NULL) {
            gmx_fatal(FARGS,"Failed to open output file, '%s'\n",xyz );
        }
        fprintf(xyzout,"%7i  %s",i_index,xyz);
        for (int n=0;n<i_index;n++)
        {
            char name[10];
            int i = ind_index[n];
            resid = top.atoms.atom[i].resind;
            strcpy(name,*top.atoms.atomname[i]);
            strcat(name,".");
            strcat(name,*top.atoms.resinfo[resid].name);
            fprintf(xyzout,"\n%7i %9s %13.4f %13.4f %13.4f ",i+1,name,fr.x[i][XX]*10,fr.x[i][YY]*10,fr.x[i][ZZ]*10);
            fprintf(xyzout,"%7i ",parmid[i]);
            for (int j = 0; j < bonds[i+1].nb ; j++ ) {
                fprintf(xyzout,"%7i ",bonds[i+1].bond[j]);
            }
        }
        fclose(xyzout);
        framen++;
    } while(read_next_frame(oenv, status, &fr));

    printf("\n");
    thanx(stderr);
    return 0;
}


int main(int argc, char *argv[])
{
    gmx_gmx2xyz(argc, argv);
    return 0;
}
