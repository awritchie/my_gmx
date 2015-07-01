#include "gmx_tilt.h"

int gmx_tilt(int argc, char *argv[])
{
    const char      *desc[] = {
        "\tLook at the tilt angle of RalGDS when docked to different\n",
        "GTPases: p21Ras (1LFD.pdb) and Rap1a (1GUA.pdb).",
        "\tUse -s for the simulated molecule topology.  Use -fr for\n",
        "the reference molecule topology.\n",
        "\n",
        "Workflow:\n",
        "\t1) VMD STAMP Structural alignment of reference structure to\n",
        "\t   1LFD.pdb chain B (Ras).\n",
        "\t2) Make a .tpr of the reference structure (which has already\n",
        "\t   been aligned to 1LFD.pdb, per step 1)\n"
        "\t3) Kabsch alignment of a single frame of the simulation\n",
        "\t   trajectory to the GTPase of the reference structure.\n",
        "\t4) Use trjconv to align the GTPase of the simulated trajectory\n",
        "\t    to the (aligned) frame from step 3.\n",
        "\t5) Run this code.\n",
    };

    gmx_bool        bVerbose = FALSE;
    const char      *struct_file, *traj_file;
    const char      *ndx_file, *xvg_file;
    const char      *ref_file;
    t_tiltdata  data;
    t_pargs         pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be slightly more verbose"},
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        //{efGRO, "-fr", NULL, ffREAD},
        {efTPS, "-fr", NULL, ffREAD},
        {efXVG, "-o","tilt",ffOPTWR},
        {efNDX, NULL, NULL, ffREAD},
    };
#define NFILE asize(fnm)
#define NPA asize(pa)
    output_env_t    oenv;
    int             ngrps, nrefgrps;
    t_topology      struct_top, ref_top;
    t_atoms         *struct_atoms=NULL, *ref_atoms=NULL;
    t_trxframe      struct_fr, ref_fr;
    t_trxstatus     *struct_status, *ref_status;
    rvec            *struct_xtop, *ref_xtop;
    matrix          struct_box, ref_box;
    int             struct_ePBC, ref_ePBC;
    int             struct_flags=TRX_READ_X, ref_flags=TRX_READ_X;
    char            buffer[1024];

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc, argv, 
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | 
                      PCA_TIME_UNIT | PCA_BE_NICE | PCA_CAN_TIME,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);

    /* Get inputs */
    struct_file = ftp2fn(efTPS, NFILE, fnm);
    traj_file   = opt2fn( "-f", NFILE, fnm);
    ref_file    = opt2fn("-fr", NFILE, fnm); 
    ndx_file    = ftp2fn(efNDX, NFILE, fnm);
    
    /* Open inputs */
    read_tps_conf(struct_file, buffer, &struct_top, &struct_ePBC,
                  &struct_xtop, NULL, struct_box, TRUE);
    sfree(struct_xtop);
    struct_atoms = &struct_top.atoms;

    read_tps_conf(ref_file, buffer, &ref_top, &ref_ePBC,
                  &ref_xtop, NULL, ref_box, TRUE);
    ref_atoms = &ref_top.atoms;

    /* open xvg file */
    data.fp = NULL;
    data.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Center of Mass displacement","Time[ps]","Tilt Angle[degrees]", oenv);
    const char *flegend[] = {
                              "Time[ps]",
                              "tilt[deg]",
                              "rmsd[nm]",
                              "CoM displacement[nm]",
                              "Rotation axis x-component[nm]",
                              "Rotation axis y-component[nm]",
                              "Rotation axis z-component[nm]",
                              "X-axis rotation[deg]",
                              "Y-axis rotation[deg]",
                              "Z-axis rotation[deg]",
                            };
    xvgr_legend(data.fp,asize(flegend),flegend,oenv);

    /* Get index information */
    static int      n_structatoms, n_refatoms;
    atom_id         *ind_structatoms, *ind_refatoms;
    char            *gn_structatoms, *gn_refatoms;
    int             i_structatoms, i_refatoms;
    std::cout<< "\nSelect group from reference structure <" << ref_file << ">." << std::endl;
    get_index(ref_atoms, ndx_file, 1, &n_refatoms, &ind_refatoms, &gn_refatoms);
    std::cout<< "\nSelect group from simulation structure(s) <" << struct_file <<"> to compare to the reference structure, <" << ref_file << ">." << std::endl;
    get_index(struct_atoms, ndx_file, 1, &n_structatoms, &ind_structatoms, &gn_structatoms);

    if ( n_structatoms != n_refatoms ) 
    {
        gmx_fatal(FARGS, "\nThe number of atoms in the index groups are different.");
    }
     
    /* read the reference file and get coordinates */
    std::vector<float> ref_coords(n_refatoms*3,0);
    readCoords(n_refatoms, ind_refatoms, ref_xtop, &ref_top, ref_coords, bVerbose);

    /* read trajectory file and loop through frames */
    int nframes = 0;
    read_first_frame(oenv, &struct_status, traj_file, &struct_fr, struct_flags);
    do {
        std::vector<float> struct_coords(n_structatoms*3,0);
        readCoords(n_structatoms, ind_structatoms, &struct_fr, &struct_top, struct_coords, bVerbose);
        displacement(ref_coords, struct_coords, data, bVerbose);
        kabsch_alignment(ref_coords,struct_coords, data, bVerbose);
        data.time.push_back(struct_fr.time);
        nframes++;
    } while(read_next_frame(oenv, struct_status, &struct_fr));
    for (int i=0; i<nframes; i++)
    {
        fprintf(data.fp,"%8.2f %8.2f %8.3f %8.3f %12.6f %12.6f %12.6f %8.2f %8.2f %8.2f\n",data.time[i],data.rotation[i],data.rmsd[i],data.distance[i],data.x_rotation_axis[i],data.y_rotation_axis[i],data.z_rotation_axis[i],data.x_rotation[i],data.y_rotation[i],data.z_rotation[i]);
    }

    return 0;
}

void average_coordinate(std::vector<float> coords, std::vector<float> &xyz)
{
    int natoms = coords.size() / 3;
    for (int i=0; i<natoms;i++)
    {
        for (int j=0; j<3; j++)
        {
            xyz[j] += coords[i+j*natoms] / natoms;
        }
    }
    return;
}

void displacement( std::vector<float> reference, std::vector<float> frame, t_tiltdata &data, gmx_bool bVerbose )
{
    std::vector<float> reference_xyz(3,0), frame_xyz(3,0);
    average_coordinate(reference, reference_xyz);
    average_coordinate(frame, frame_xyz);
    float x, y, z, d;
    x = reference_xyz[0] - frame_xyz[0];
    y = reference_xyz[1] - frame_xyz[1];
    z = reference_xyz[2] - frame_xyz[2];
    d = sqrt(x*x+y*y+z*z);
    data.distance.push_back( d );
    //std::cout << "\n{ " << reference_xyz[0] << " " << reference_xyz[1] << " " << reference_xyz[2] << " } { " << frame_xyz[0] << " " << frame_xyz[1] << " " << frame_xyz[2] << " } -> " << d << "";
    if (bVerbose)
    {
        fprintf(stdout,"  Distance between center of mass of each selection: %8.4f nm\n",d);
    }
    return;
}

/* Reading a TPX */
void readCoords(int n_atoms, atom_id ind_atoms[], rvec *x, t_topology *top, std::vector<float> &coords, gmx_bool bVerbose)
{
    Matrix coords3N(n_atoms);
    float mass = 0;
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        std::vector<float> xyz(3,0);
        xyz[0] = x[n][0];
        xyz[1] = x[n][1];
        xyz[2] = x[n][2];
        mass += top->atoms.atom[n].m;
        coords3N[i] = xyz;
    }
    int nn = 0;
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        for (int j=0; j<3; j++)
        {
            //coords3N[i][j] *= top->atoms.atom[n].m;
            //coords3N[i][j] /= mass;
            coords[i+j*n_atoms] = coords3N[i][j];
        }
    }
    if (bVerbose)
    {
        //std::cout << "\n\n Mass Weighted Atom Coordinates: ";
        std::cout << "\n\n Coordinates: \tMass\tx\ty\tz";
        for (int i=0; i<n_atoms; i++)
        {
            std::cout << "\n\t" << *top->atoms.atomname[ind_atoms[i]];
            std::cout << "\t" << top->atoms.atom[ind_atoms[i]].m << "\t" << coords3N[i][0] << "\t" << coords3N[i][1] << "\t" << coords3N[i][2];
       }
        std::cout << "\n";
    }
    return;
}
/* Reading a trajectory */
void readCoords(int n_atoms, atom_id ind_atoms[], t_trxframe *fr, t_topology *top, std::vector<float> &coords, gmx_bool bVerbose)
{
    Matrix coords3N(n_atoms);
    float mass = 0;
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        std::vector<float> xyz(3,0);
        xyz[0] = fr->x[n][XX];
        xyz[1] = fr->x[n][YY];
        xyz[2] = fr->x[n][ZZ];
        mass += top->atoms.atom[n].m;
        coords3N[i] = xyz;
    }
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        for (int j=0; j<3; j++)
        {
            //coords3N[i][j] *= top->atoms.atom[n].m;
            //coords3N[i][j] /= mass;
            coords[i+j*n_atoms] = coords3N[i][j];
        }
    }
    if (bVerbose)
    {
        //std::cout << "\n\n Mass Weighted Atom Coordinates: ";
        std::cout << "\n\n Coordinates: \tMass\tx\ty\tz";
        for (int i=0; i<n_atoms; i++)
        {
            std::cout << "\n\t" << *top->atoms.atomname[ind_atoms[i]];
            std::cout << "\t" << top->atoms.atom[ind_atoms[i]].m << "\t" << coords3N[i][0] << "\t" << coords3N[i][1] << "\t" << coords3N[i][2];
        }
        std::cout << "\n";
    }
    return;
}

/* Kabsch alignment */
void kabsch_alignment( std::vector<float> ref, std::vector<float> tar, t_tiltdata &data, gmx_bool bVerbose)
{
    if (ref.size() != tar.size())
    {
        std::cerr << "\nError! Sizes of reference coordinate matrix and simulated structure coordinate matrices do not match!" << std::endl;
        std::exit(1);
    }
    int ncoords = ref.size();
    int natoms = ncoords/3;
    // Center the two selections
    std::vector<float> stsel1(ncoords,0), stsel2(ncoords,0), stsel2T(ncoords,0);
    std::vector<float> ref_com(3,0), tar_com(3,0);
    average_coordinate(ref, ref_com);
    average_coordinate(tar, tar_com);
    for (int i=0; i<natoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            stsel1[i+j*natoms] = ref[i+j*natoms] - ref_com[j];
            stsel2[i+j*natoms] = tar[i+j*natoms] - tar_com[j];
        }
    }
    // Initial residual
    float E0 = sdot(ncoords,&stsel1[0],1,&stsel1[0],1)+sdot(ncoords,&stsel2[0],1,&stsel2[0],1) ;
    // dot(target_transpose,reference)
    std::vector<float> T1_dot_2(3*natoms,0);
    sgemm('T','N',3,natoms,natoms,1,&stsel2[0],natoms,&stsel1[0],natoms,1,&T1_dot_2[0],3);
    // SVD of the dot product
    std::vector<float> U(9,0), S(3,0), V(9,0), work(5*9,0);
    int info;
    sgesvd('A','A',3,3,&T1_dot_2[0],3,&S[0],&U[0],3,&V[0],3,&work[0],9*5,info);
    /*std::cout << "\n S: ";
    for (int i=0;i<3;i++)
    {
        std::cout << S[i] << " ";
    }
    std::cout << "\n U: ";
    for (int i=0;i<9;i++)
    {
        std::cout << U[i] << " ";
    }*/
    float reflect = det3x3(&U[0]) * det3x3(&V[0]);
    if ( 1 - reflect > 1e-5)
    {
        S[2] = -S[2];
        U[6] = -U[6];
        U[7] = -U[7];
        U[8] = -U[8];
    }
    float rmsd = sqrt(fabs(
                           E0
                           - (2.0 *
                              (S[0]+S[1]+S[2])
                              )
                           )
                      /natoms);
    // Rotation matrix is dot(U,V)
    std::vector<float> M(9,0);
    sgemm('N','N',3,3,3,1,&U[0],3,&V[0],3,1,&M[0],3);
    /*
     M = [ 0 3 6 ] = [ 00 01 02 ]
         [ 1 4 7 ]   [ 10 11 12 ]
         [ 2 5 8 ]   [ 20 21 22 ]
     */
    float trace = M[0]+M[4]+M[8];
    float angle = acos((trace-1)/2)*RAD2DEG;
    float rx,ry,rz,ux,uy,uz;
    rx = atan2(M[5],M[8])*RAD2DEG;
    ry = atan2(-M[2],sqrt(M[5]*M[5]+M[8]*M[8]))*RAD2DEG;
    rz = atan2(M[1],M[0])*RAD2DEG;
    float zeta = sqrt(
                        (M[5]-M[7])*(M[5]-M[7])
                      + (M[6]-M[2])*(M[6]-M[2])
                      + (M[3]-M[1])*(M[3]-M[1])
                      );
    //std::cout << "\n" << M[5] << " - " << M[7] << " = " << M[5]-M[7];
    //std::cout << "\n" << M[6] << " - " << M[2] << " = " << M[6]-M[2];
    //std::cout << "\n" << M[3] << " - " << M[1] << " = " << M[3]-M[1] << std::endl;
    ux = (M[5]-M[7])/zeta;
    uy = (M[6]-M[2])/zeta;
    uz = (M[3]-M[1])/zeta;
    //std::cout << zeta << " { " << ux << " " << uy << " " << uz << " }" << sqrt(ux*ux+uy*uy+uz*uz) << std:: endl;
    if (bVerbose)
    {
        fprintf(stdout,"%12s%12s%12s%12s%12s%12s%12s%12s\n","Angle(deg)","rmsd(nm)","x(deg)","y(deg)","z(deg)","ux(nm)","uy(nm)","uz(nm)");
        fprintf(stdout,"%12.3f%12.6f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n",angle,rmsd,rx,ry,rz,ux,uy,uz);
    }
    data.rotation.push_back(angle);
    data.rmsd.push_back(rmsd);
    data.x_rotation.push_back(rx);
    data.y_rotation.push_back(ry);
    data.z_rotation.push_back(rz);
    data.x_rotation_axis.push_back(ux);
    data.y_rotation_axis.push_back(uy);
    data.z_rotation_axis.push_back(uz);

    return;
}

