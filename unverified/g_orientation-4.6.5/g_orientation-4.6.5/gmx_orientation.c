//
//  gmx_orientation.c
//  g_orientation-4.6.5
//
//  Created by Andrew Ritchie on 2/28/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#include "gmx_orientation.h"



int gmx_orientation(int argc, const char * argv[]){
    const char *desc[] = {
        "    This utility calculates the the angles of a dipole ",
        "vector off of the moment of inertia tensor's ",
        "eigenvectors.  An index file containing ",
        "only the atom indices to be used for the MoI tensor ",
        "is retured--dipole moment vector is a little more ",
        "complex than I anticipated and has been put on the ",
        "wayside for now.",
        "    The -sp flag causes the vector angles in spherical-",
        "polar coordinates to be returned, instead of the angle",
        "off each eigenvector.",
    };
    t_topology top;
    int        ePBC;
    char       title[STRLEN];
    t_trxframe fr;
    t_trxstatus *status;
    rvec       *xtop;
    matrix     box;
    int        i_inertia, i_dipole;
    int        flags = TRX_READ_X;
    atom_id    *ind_inertia, *ind_dipole;
    t_atoms    *atoms=NULL;
    char       *gn_inertia, *gn_dipole;
    int        a1=-1, a2=-1;    // Initializing these negative as a check
    rvec       bondvector;
    rvec       n_bondvector;
    FILE       *outputfile;
    output_env_t    oenv;
    const char      *tpr_file, *xtc_file, *index_file, *out_file = NULL;
    gmx_bool       bVerbose=FALSE, docorr=FALSE, spherical_polar=TRUE;
    
    
    //  We need to know which atoms to use as the bond vector
    t_pargs pa[] = {
        { "-a1", TRUE, etINT,
            {&a1}, "Atom number 1 (must be set!)"
        },
        { "-a2", TRUE, etINT,
            {&a2}, "Atom number 2 (must be set!)"
        },
        { "-v", FALSE, etBOOL,
            {&bVerbose}, "Be slightly more verbose"
        },
        { "-corr", FALSE, etBOOL,
            {&docorr}, "Calculate the dipole moment correlation function"
        },
        { "-sp", TRUE, etBOOL,
            {&spherical_polar}, "Print angles in spherical polar coordinates (Rather than angle off each eigenvector)"
        }
    };
    
    t_filenm fnm[] = {
        { efTPS,  NULL,  NULL, ffREAD },    /* this is for the topology */
        { efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
        { efXVG, "-o", "dipang", ffWRITE },     /* and this is for the output file */
        { efNDX, NULL, NULL, ffREAD },   /* and finally, the index file */
        { efXVG, "-c", "dipcorr", ffOPTWR },     /* If -corr is flagged, we need an output */
    };
    
#define NFILE asize(fnm)
#define NPA asize(pa)
    
    CopyRight(stderr,argv[0]);
    
    /* This is the routine responsible for adding default options,
     * calling the X/motif interface, etc. */
    parse_common_args(&argc,argv,
                      PCA_CAN_BEGIN | PCA_CAN_TIME | PCA_CAN_VIEW |
                      PCA_TIME_UNIT | PCA_BE_NICE,
                      NFILE,fnm,NPA,pa,asize(desc),desc,
                      0,NULL, &oenv);
    tpr_file = ftp2fn(efTPS,NFILE,fnm);
    xtc_file = opt2fn("-f",NFILE,fnm);
    out_file = opt2fn("-o",NFILE,fnm);
    
    read_tps_conf(tpr_file,title,&top,&ePBC,&xtop,NULL,box,TRUE);
    sfree(xtop);
    atoms=&top.atoms;
    
    //i_something is number of members in that index
    //ind_something is the atom indices in the index
    printf("Select group for generating the moment of inertia tensor:\n");
    get_index(atoms, ftp2fn(efNDX,NFILE,fnm),1,&i_inertia,&ind_inertia,&gn_inertia);
    
    //  Calculating the dipole moment vector needs connectivity information, which I'm not sure
    //  how gromacs handles.  I don't want to delve too deeply into this right now,
    //  especially since everything I need can be approximated very easily with intelligent
    //  choice of a1 and a2
    //    printf("Select group for generating the dipole vector:\n");
    //    get_index(atoms, ftp2fn(efNDX,NFILE,fnm),1,&i_dipole,&ind_dipole,&gn_dipole);
    
    /* The first time we read data is a little special */
    read_first_frame(oenv, &status,xtc_file,&fr,flags);
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
    
    /* output file */
    outputfile=ffopen( opt2fn("-o",NFILE,fnm),"w" );
    if(outputfile==NULL) {
        gmx_fatal(FARGS,"Failed to open output file '%s'\n", opt2fn("-o",NFILE,fnm) );
    }
    
    double dip_comp[100000][3]; // arbitrarily large
    int frame=0;
    
    do {
        // Build the tensor from group 0 selected
        double tensor[3][3]={{0}};
        
        for (int i=0 ; i<=i_inertia ; i++ ) {
            // look at the indexed atoms
            int n=ind_inertia[i];
            // diagonals first
            tensor[0][0]+=top.atoms.atom[n].m*(fr.x[n][YY]*fr.x[n][YY]+fr.x[n][ZZ]*fr.x[n][ZZ]);
            tensor[1][1]+=top.atoms.atom[n].m*(fr.x[n][XX]*fr.x[n][XX]+fr.x[n][ZZ]*fr.x[n][ZZ]);
            tensor[2][2]+=top.atoms.atom[n].m*(fr.x[n][XX]*fr.x[n][XX]+fr.x[n][YY]*fr.x[n][YY]);
            // off diagonals; tensor is symmetric
            tensor[0][1]-=top.atoms.atom[n].m*fr.x[n][XX]*fr.x[n][YY];
            tensor[0][2]-=top.atoms.atom[n].m*fr.x[n][XX]*fr.x[n][ZZ];
            tensor[1][2]-=top.atoms.atom[n].m*fr.x[n][YY]*fr.x[n][ZZ];
        }
        
        // find the moment of inertia tensor's eigenvectors
        real * eigenvectors;
        real * eigenvalues;
        real * full_tensor;
        snew(eigenvectors,9);
        snew(eigenvalues,3);
        snew(full_tensor,9);
        
        // Gromacs' eigensolver wants a 1D matrix such that m[i][j]=>m[i]=(m[i1][j1], m[i1][j2], m[i1][j3],...).
        // As such, this can only handle symmetric matrices.  Fortunately, the moment of inertia tensor IS symmetric
        int count=0;
        for (int i=0; i<3;i++)
        {
            for (int j=0; j<3;j++)
            {
                full_tensor[count]=tensor[i][j];
                count ++;
            }
        }
        eigensolver(full_tensor,3,0,3,eigenvalues,eigenvectors);
        // I should have created a structure here; instead I did it the long way.
        rvec e1={eigenvectors[0],eigenvectors[3],eigenvectors[6]};
        rvec e2={eigenvectors[1],eigenvectors[4],eigenvectors[7]};
        rvec e3={eigenvectors[2],eigenvectors[5],eigenvectors[8]};
        
        if (bVerbose)
        {
            printf("\nUpper-Triangular Tensor:\n");
            for (int i=0; i<3 ; i++ ) {
                printf("[");
                for (int j=0; j<3; j++) {
                    printf("\t%15.3f ",tensor[i][j]);
                }
                printf("]\n");
            }
            printf("Eigenvalues: ");
            
            for (int i=0;i<3;i++)
            {
                printf("%.3f ",eigenvalues[i]);
            }
            printf("\n");
            printf("Eigenvectors: ");
            for (int i=0;i<3;i++)
            {
                printf("[%.3f %.3f %.3f]",eigenvectors[i], eigenvectors[i+3], eigenvectors[i+6]);
            }
            printf("\n");
        }
        
        // calculate the dipole moment vector from group 1 selected
        //  Connectivity information is needed to calculate a dipole moment, so instead
        //  I'm using bond vectors, since the nitrile is linear and the side chains I've
        //  thus-far been interested in can be approximated very easily with intelligent
        //  choice of starting and ending atoms.
        //        rvec dipole_vector=fr.x[ind_dipole[0]];//, index_atom={0,0,0};
        //        double x=0, y=0, z=0;
        //        for (int i=1 ; i < i_dipole ; i++ ) {
        //            x=-1*fr.x[ind_dipole[i]][XX];//*top.atoms.atom[ind_dipole[i]].q;
        //            y=-1*fr.x[ind_dipole[i]][YY];//*top.atoms.atom[ind_dipole[i]].q;
        //            z=-1*fr.x[ind_dipole[i]][ZZ];//*top.atoms.atom[ind_dipole[i]].q;
        //            printf("%s: %f %f %f, q=%f\n", *(top.atoms.atomname[ind_dipole[i]]),x,y,z,top.atoms.atom[ind_dipole[i]].q);
        //            rvec index_atom={x,y,z};
        //            rvec_add( dipole_vector, index_atom, dipole_vector );
        //        }
        
        rvec_sub(fr.x[a2],fr.x[a1],bondvector);
        unitv( bondvector, n_bondvector);
        
        // Again, had I created a structure to hold the eigenvectors, this could be a loop.
        double angle1=my_acos(iprod(n_bondvector,e1));
        double angle2=my_acos(iprod(n_bondvector,e2));
        double angle3=my_acos(iprod(n_bondvector,e3));
        
        if (bVerbose)
        {
            printf("Looking at bond vector from %s to %s\n",*(top.atoms.atomname[a1]),*(top.atoms.atomname[a2]));
            printf("Bond unit vector is: [%f %f %f]\n", n_bondvector[0], n_bondvector[1], n_bondvector[2]);
        }
        
        double rho = pow(iprod(bondvector,bondvector),.5);
        double theta = angle1;
        double phi1, phi2, phi;
        if (spherical_polar) {
            rvec vector, nvector, v2, v3;
            projection(e2, bondvector, v2); //component along eigenvector 2
            projection(e3, bondvector, v3); //component along eigenvector 3
            rvec_add(v2,v3,vector); //sum of the components is the total vector
            unitv(vector, nvector);
            unitv(e2,e2);
            unitv(e3,e3);
            phi1 = my_acos( iprod(e2,nvector) );
            phi2 = my_acos( iprod(e3,nvector) );
            
            if ( phi2 > 90. ) { phi=-1*phi1; }
            else { phi = phi1; }
            
            if (bVerbose)
            {
                printf("rho(nm) = %6.3f\ttheta(degrees) = %6.3f\tphi(degrees) = %6.3f\n\n", rho, theta, phi);
            }
            
            fprintf( outputfile, "%9.4f %9.4f %9.4f\n", rho, theta, phi );
        }
        else {
            if (bVerbose)
            {
                printf("Angles in degrees:\t%6.3f\t%6.3f\t%6.3f\n\n", angle1, angle2, angle3);
            }
            //fprintf( outputfile, "%-10.4f\t %9.4f %9.4f %9.4f\n", fr.time, angle1, angle2, angle3 );
            fprintf( outputfile, "%9.4f %9.4f %9.4f\n", angle1, angle2, angle3 );
            
        }
        
        if (docorr)
        {
            if (spherical_polar) {
                dip_comp[frame][0]=rho;
                dip_comp[frame][1]=cos(theta*M_PI/180);
                dip_comp[frame][2]=cos(phi*M_PI/180);
                //printf("\nACF: %f %f %f", dip_comp[frame][0], dip_comp[frame][1], dip_comp[frame][2]);
            }
            else {
                dip_comp[frame][0]=cos(angle1*M_PI/180);
                dip_comp[frame][1]=cos(angle2*M_PI/180);
                dip_comp[frame][2]=cos(angle3*M_PI/180);
                //printf("\n%f %f %f", dip_comp[frame][0], dip_comp[frame][1], dip_comp[frame][2]);
            }
        }
        else free(dip_comp);
        frame++;
        
    } while(read_next_frame(oenv,status,&fr));
    fclose( outputfile );
    
    if (docorr && frame > 1)
    {
        /* autocorrelation output file */
        FILE *corroutputfile;
        corroutputfile=ffopen( opt2fn("-c",NFILE,fnm),"w" );
        if(corroutputfile==NULL)
        {
            gmx_fatal(FARGS,"Failed to open correlation output file '%s'\n", opt2fn("-c",NFILE,fnm) );
        }
        
        // Need to know the average for each dipole vector component
        double average[3]={0};
        for (int i=0; i<frame ; i++ )
        {
            for (int j=0; j<3; j++ )
            {
                average[j]+=dip_comp[i][j]/frame;
            }
        }
        
        if (bVerbose)
        {
            printf("\nAverage values: %3.8f %3.8f %3.8f\n",average[0],average[1],average[2]);
        }
        
        // Get the denominator, since it is the same for all lags
        // denominator = sum_i( (value_i-average)^2 )
        double denominator[3]={0};
        for (int j=0; j<3; j++)
        {
            for (int i=0; i<frame; i++)
            {
                denominator[j]+=(dip_comp[i][j]-average[j])*(dip_comp[i][j]-average[j]);
            }
            if (bVerbose)
            {
                printf("%3.8f ", denominator[j]);
            }
        }
        
        // Get the numerator for each lag and calculate the autocorrelation coefficient
        // numerator = sum_i( (value_i - average)*(value_i+n - average ) )
        // This appears to be fairly consistant with python's numpy.correlate:
        // def correlation(x):
        //     result=numpy.correlation(x,x,mode='full')
        //     return result[result.size/2:]
        // need to dynamically allocate the array
        double **numerator = malloc(sizeof *numerator * 3);
        if (numerator)
        {
            for (int i=0;i<3;i++)
            {
                numerator[i] = malloc(sizeof *numerator[i] * (frame+1));
            }
        }
        
        for (int j=0;j<4;j++)
        {
            for (int k=0;k<frame;k++)
            {
                numerator[j][k]=0.0 ;
            }
        }
        
        for (int k=0; k<frame; k++)
        {
            fprintf( corroutputfile, "%i ", k);
            for (int j=0; j<3 ; j++)
            {
                for (int i=0; i<frame-k ; i++)
                {
                    if (bVerbose)
                    {
                        printf("j=%i k=%i, numerator[j][k] = %f\n",j,k,numerator[j][k]);
                    }
                    numerator[j][k]+=(dip_comp[i][j]-average[j])*(dip_comp[i+k][j]-average[j]);
                    if (bVerbose)
                    {
                        printf("i=%i, j=%i, k=%i : (dip_comp[i][j]-average[j])*(dip_comp[%i+k][j]-average[j]) = ( %f - %f ) * ( %f - %f ) = %f \n ", i,j,k,i,dip_comp[i][j],average[j],dip_comp[i+k][j],average[j], numerator[j][k] );
                        printf("%f\n",(dip_comp[i][j]-average[j])*(dip_comp[i+k][j]-average[j]));
                    }
                }
                if (bVerbose)
                {
                    printf("numerator[j][k]/denominator[j] = %f / %f = %f \n", numerator[j][k],denominator[j],numerator[j][k]/denominator[j]);
                }
                fprintf( corroutputfile, "%f ", numerator[j][k]/denominator[j]);
            }
            fprintf(corroutputfile,"\n");
            if (bVerbose)
            {
                printf("\n");
            }
        }
        
        fclose(corroutputfile);
    }
    else if ( docorr && frame <= 1 )
    {
        printf("\nYou indicated that you would like to calculate the autocorrelation, but there ");
        printf("is only 1 frame... Either you submited a structure file instead of a trajectory, ");
        printf("or your trajectory is only a single frame.  Either way, reevaluate what you are ");
        printf("trying to do...\n");
    }
    
    if (bVerbose){printf("\n");}
    
    thanx(stderr);
    

    
    return 0;
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

double my_acos( double ip ) {
    if ( ip >= 1 ) {
        return 0.0;
    }
    else if ( ip <= -1 ) {
        return 180.0;
    }
    else {
        return acos( ip )*180/M_PI;
    }
}