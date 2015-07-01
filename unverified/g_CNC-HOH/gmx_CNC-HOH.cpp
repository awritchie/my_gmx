#include "gmx_CNC-HOH.h"

int gmx_CNCHOH(int argc, char *argv[])
{
    const char      *desc[] = {
        "Determine whether or not water is hydrogen bonding to the",
        "nitrile and how much water is persistent",
    };

    const char      *flegend[] = {
        "Persistent water", "Persistent h-bond",
        "Total persistent water", "Total persistent h-bond"};
    gmx_bool        bVerbose = FALSE;
    const char      *file, *traj_file;
    const char      *ndx_file, *xvg_file;
    t_hbdata        data;
    /*
     Hydrogen-bond acceptor properties of nitriles: a combined crystallographic and ab initio theoretical investigation
     Jean-Yves Le Questel, Michel Berthelot and Christian Laurence
    r average = .205(.020)
    theta1 = CNH, average = 145(23)
    theta2 = NHO, average = 156(18)
    95% of data within +/- 2 * STD
    */
    parms           cutoffs;
    cutoffs.r       = 0.245;
    cutoffs.cnh     = 145 - (2*23);
    cutoffs.nho     = 156 - (2*18);
    cutoffs.rmax    = 0.5;
    t_pargs         pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be slightly more verbose"},
        { "-NHO_cutoff", FALSE, etREAL, {&cutoffs.nho},
          "Minimum N-H-O bond angle to count as a hydorgen bond.  Default=(156 - 2*18) degrees"},
        { "-CNH_cutoff", FALSE, etREAL, {&cutoffs.cnh},
          "Minimum N-H-O bond angle to count as a hydorgen bond.  Default=(145 - 2*23) degrees"},
        { "-NH_cutoff", FALSE, etREAL, {&cutoffs.r},
          "Minimum N-H bond distance to count as a hydrogen bond.  Default=(0.205 + 2*.020) nm"},
        { "-cutoff", FALSE, etREAL, {&cutoffs.rmax},
          "Minimum distance from CNC nitrogen to count as \"near\" the probe.  Default=0.50 nm"},
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        {efXVG, "-o","CNC-HOH",ffWRITE},
        {efXVG, "-op","persistent_CNC-HOH",ffWRITE},
        {efLOG, "-oa","residues",ffWRITE},
        {efNDX, NULL, NULL, ffREAD},
    };
    gmx_bool doLog = false, doPersistent = false;
#define NFILE asize(fnm)
#define NPA asize(pa)
    output_env_t    oenv;
    int             ngrps, nrefgrps;
    t_topology      top;
    t_atoms         *atoms=NULL;
    t_trxframe      fr;
    t_trxstatus     *status;
    matrix          box;
    int             ePBC;
    rvec            *xtop;
    int             flags=TRX_READ_X;
    char            buffer[1024];

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc, argv, 
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | 
                      PCA_TIME_UNIT | PCA_BE_NICE | PCA_CAN_TIME,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);

    /* Get inputs */
    file        = ftp2fn(efTPS, NFILE, fnm);
    traj_file   = opt2fn( "-f", NFILE, fnm);
    ndx_file    = ftp2fn(efNDX, NFILE, fnm);
    
    /* Open inputs */
    read_tps_conf(file, buffer, &top, &ePBC,
                  &xtop, NULL, box, TRUE);
    sfree(xtop);
    atoms = &top.atoms;

    /* open xvg files */
    data.fp = NULL;
    data.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Trajectory Hydrogen Bonds","Time[ps]", "Number of Hydrogen Bonds", oenv);
    
    data.fpp = NULL;
    if (opt2bSet("-op", NFILE, fnm))
    {
        data.fpp = xvgropen(opt2fn("-op", NFILE, fnm), "Persistent Waters and Hydrogen Bonds", "Number of Frames", "Number Water Molecules", oenv);
        xvgr_legend(data.fpp,asize(flegend), flegend, oenv);
        doPersistent = true;
    }
    
    data.fpa = NULL;
    std::ofstream outputfile;
    if (opt2bSet("-oa", NFILE, fnm))
    {
        data.fpa = opt2fn("-oa", NFILE, fnm);
        outputfile.open(data.fpa);
        doLog = true;
    }
    /* read trajectory file and loop through frames */
    int i=0;
    std::vector<t_water> water;
    data.nhb = 0;
    int nhb = 0;
    read_first_frame(oenv, &status, traj_file, &fr, flags);
    do {
        int frame_hbond = 0; // per frame number of hbonds, 0 or 1
        gmx_bool hbond = false;
        water.push_back(analyze_frame(&fr, &top, cutoffs));
        // Count the number of waters within cutoff and hydrogen bonding
        outputfile << " " << traj_file << "\n [ frame" << i << " ]\n";
        
        data.nps = std::vector<int> ((int)water[i].size(),0);
        data.nhbs = std::vector<int> ((int)water[i].size(),0);
/*
 * Making threads is too expensive...
 */
//#pragma omp parallel for num_threads(8) reduction(+ : nhb)
        for (int j=0; j<(int)water[i].size(); j++)
        {
            if (i > 0) // Not persistent in first frame
            {
                if (water[i][j].rstatus > 0)
                {
                    if (water[i-1][j].rstatus > 0)
                    {
                        data.nps[j]++;
                        if (doLog)
                        {
                            outputfile << water[i][j].resid << " ";
                        }
                    }
                    if (water[i][j].rstatus == 2)
                    {
                        data.nhbs[j]++;
                        if (doLog)
                        {
                            outputfile << "(" << water[i][j].resid << ") ";
                        }
                        if (hbond == false) // only 1 hbond per frame while we still allow persistent checks
                        {
                            nhb++;
                            frame_hbond = 1;
                            hbond = true;
                        }
                    }
                }
            }
            else
            {
                // Make sure the array size of data.nps and .nhbs is the number of water molecules
                //data.nps.push_back(0);
                //data.nhbs.push_back(0);
                if (water[i][j].rstatus > 0)
                {
                    data.nps[j] = 1;
                    {
                        //outputfile << water[i][j].resid << " ";
                    }
                    if (water[i][j].rstatus == 2)
                    {
                        data.nhbs[j] = 1;
                        {
                            outputfile << "(" << water[i][j].resid << ") ";
                        }
                        if (hbond == false) // only 1 h-bond per frame while we still allow persistent checks
                        {
                            nhb = 1;
                            frame_hbond = 1;
                            hbond = true;
                        }
                    }
                }
            }
        }
        if (data.fp)
        {
            fprintf(data.fp,"%10.3f %5i\n",fr.time,frame_hbond);
        }
        i++;
    } while(read_next_frame(oenv, status, &fr));

    data.nhb = nhb;

    double nframes = i;
    double percent = data.nhb/nframes * 100.;
    
    std::vector<int> nps(i,0);
    std::vector<int> tnps(i,0);
    std::vector<int> nhbs(i,0);
    std::vector<int> tnhbs(i,0);
    for (int i=0;i<data.nps.size();i++) // for each water molecule, check the number of times it was persistant, and add that to the appropriate entrie in nps and nhbs
    {
        nps[data.nps[i]]++;
        nhbs[data.nhbs[i]]++;
    }
    int sump = 0, sumhb = 0;
    for (int i=nps.size()-1; i>=0; i--){
        sump += nps[i];
        sumhb += nhbs[i];
        tnps[i] = sump;
        tnhbs[i] = sumhb;
    }
    if (doPersistent)
    {
        for (int i=0;i<nps.size();i++)
        {
            fprintf(data.fpp,"%10i %10i %10i %10i %10i\n",i,nps[i],nhbs[i],tnps[i],tnhbs[i]);
        }
    }
    fprintf(stdout,"Total hydrogen bonding: %10i (%10.3f%% of frames)\n",data.nhb,percent);
    if (data.fp)
    {
        ffclose(data.fp);
    }
    if (doLog)
    {
        outputfile << "\n";
        outputfile.close();
    }
    return 0;
}

std::vector<double> vecsub( std::vector<double> r1, std::vector<double> r2 )
{
   std::vector<double> result(3,0);
   result[0] = r1[0] - r2[0];
   result[1] = r1[1] - r2[1];
   result[2] = r1[2] - r2[2];
   return result;
}

double veclen( std::vector<double> r)
{
   return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
}

double vecangle( std::vector<double> r1, std::vector<double> r2, std::vector<double> r3)
{
    std::vector<double> r21 = vecsub(r1,r2);
    std::vector<double> r23 = vecsub(r3,r2);
    double l_r21 = veclen(r21);
    double l_r23 = veclen(r23);
    for (int i=0; i<3; i++)
    {
        r21[i] /= l_r21;
        r23[i] /= l_r23;
    }
    double xx = r21[0]*r23[0];
    double yy = r21[1]*r23[1];
    double zz = r21[2]*r23[2];
    return acos(xx+yy+zz)*180./M_PI; 
}

int is_Hbond(std::vector<double> CD, std::vector<double> NE, std::vector<double> OW, std::vector<double> HW, parms cutoffs)
{
    double CNH = vecangle(CD,NE,HW);
    double NHO = vecangle(NE,HW,OW);
    double NH = veclen(vecsub(NE,HW));
    double NO = veclen(vecsub(NE,OW));
    if (NH < cutoffs.rmax || NO < cutoffs.rmax){
        if (NH < cutoffs.r && NHO > cutoffs.nho && CNH > cutoffs.cnh) {
            return 2;
        }
        return 1;
    }
    return 0;
}

int parse_water(std::vector<double> CD, std::vector<double> NE, std::vector<double> OW, std::vector<double> HW1, std::vector<double> HW2, parms cutoffs)
{
    int hb2HW1 = is_Hbond(CD,NE,OW,HW1,cutoffs);
    int hb2HW2 = is_Hbond(CD,NE,OW,HW2,cutoffs);
    if (hb2HW1 > 0 || hb2HW2  > 0)
    {
        if (hb2HW1 == 2 || hb2HW2 == 2)
        {
            return 2;
        }
        return 1;
    }
    
    return 0;
}

t_water analyze_frame(t_trxframe *fr, t_topology *top, parms cutoffs)
{
    std::vector<double> CD(3);
    std::vector<double> NE(3);
    std::vector<double> OW(3);
    std::vector<double> HW1(3);
    std::vector<double> HW2(3);
    t_water water;
    for (int i=0; i<top->atoms.nr; i++)
    {
        if (strncmp(*top->atoms.resinfo[ top->atoms.atom[i].resind ].name, "CNC", 4) == 0)
        {
            if (strncmp(*top->atoms.atomname[i],"CD", 3) == 0)
            {
                CD[0] = fr->x[i][XX];
                CD[1] = fr->x[i][YY];
                CD[2] = fr->x[i][ZZ];
            }
            else if (strncmp(*top->atoms.atomname[i],"NE", 4) == 0)
            {
                NE[0] = fr->x[i][XX];
                NE[1] = fr->x[i][YY];
                NE[2] = fr->x[i][ZZ];
            }
        }
        else if (strncmp(*top->atoms.resinfo[ top->atoms.atom[i].resind ].name, "SOL", 4) == 0)
        {
            if (strncmp(*top->atoms.atomname[i], "OW", 3) == 0)
            {
                t_mol mol;
                OW[0] = fr->x[i][XX];
                OW[1] = fr->x[i][YY];
                OW[2] = fr->x[i][ZZ];
                /* Assuming H1 is i+1 and H2 is i+2 */
                HW1[0] = fr->x[i+1][XX];
                HW1[1] = fr->x[i+1][YY];
                HW1[2] = fr->x[i+1][ZZ];
                
                HW2[0] = fr->x[i+2][XX];
                HW2[1] = fr->x[i+2][YY];
                HW2[2] = fr->x[i+2][ZZ];
                
                mol.resid = top->atoms.atom[i].resind + 1;
                mol.rstatus = parse_water(CD,NE,OW,HW1,HW2,cutoffs);
                water.push_back(mol);
            }
        }
    }
    return water;
}
    
