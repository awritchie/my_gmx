/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "g_nitrile_hbond.hpp"

/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata      *d = (t_analysisdata *)data;
    t_water             water;
    int sigma_nhb = 0, pi_nhb = 0, lQ_nhb = 0;
    float a1[3] = { fr->x[d->a1][0], fr->x[d->a1][1], fr->x[d->a1][2] };
    float a2[3] = { fr->x[d->a2][0], fr->x[d->a2][1], fr->x[d->a2][2] };
    /* Loop through everything */
    for (int g = 0; g < nr; ++g) { // for each group
        if (g > 1) {
            fprintf(stderr,"\nWarning: More than 1 group was specified.  Cowardly ignoring all additional groups\n");
            break;
        }
        for (int i = 0; i < sel[g]->p.nr; ++i) { // for each atom in the group
            int atom_ndx = sel[g]->g->index[i]; // how we get to the atom index
            /* Get only the water oxygen atoms since we know how water molecules are numbered */
            int resid = top->atoms.atom[atom_ndx].resind; // increment by +1 to match gro file
            char *resname = *top->atoms.resinfo[resid].name;
            if (strncmp(resname, "SOL", 4) == 0 || strncmp(resname, "HOH", 4) == 0)
            {
                char *atomname = *top->atoms.atomname[atom_ndx];
                if (strncmp(atomname, "OW", 3) == 0)
                {
                    t_mol mol;
                    mol.resid = resid;
                    mol.is_hb = 0;
                    float OW[3] = { fr->x[atom_ndx][0], fr->x[atom_ndx][1], fr->x[atom_ndx][2] };
                    float H1[3] = { fr->x[atom_ndx+1][0], fr->x[atom_ndx+1][1], fr->x[atom_ndx+1][2] };
                    float H2[3] = { fr->x[atom_ndx+2][0], fr->x[atom_ndx+2][1], fr->x[atom_ndx+2][2] };
                    if (parse_water(a1, a2, OW, H1, H2, d->cr, mol))
                    {
                        d->nhb++;
                        if (mol.nh < d->r && mol.cnh > d->cnh && mol.nho > d->nho )
                        {
                            if (d->doGeo) {
                                fprintf(d->fpg, "%10i %10i %12.4f %12.4f %12.4f\n", d->framen, atom_ndx, 10*mol.nh, mol.cnh, mol.nho);
                            }
                            lQ_nhb++;
                            mol.is_hb++;
                        }
                        if (mol.cnh >= 120.)
                        {
                            sigma_nhb++;
                        }
                        else
                        {
                            pi_nhb++;
                        }
                    }
                    water.push_back(mol);
                    d->nwater++;
                }
            }
            
        }
    }
    d->water.push_back(water);
    d->frame_lQ_nhb.push_back(lQ_nhb);
    d->frame_sigma_nhb.push_back(sigma_nhb);
    d->frame_pi_nhb.push_back(pi_nhb);
    /* increment the frame number */
    d->framen++;
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

int parse_water(const float *a1, const float *a2, const float *ow, const float *h1, const float *h2, const float &r, t_mol &mol)
{
    float v1 = vlen(a2,h1);
    float v2 = vlen(a2,h2);
    float t1 = 0, t2 = 0, d = 0;
    if (v1 <= r || v2 <= r )
    {
        if (v1 < v2)
        {
            mol.nh = v1;
            mol.cnh = vangle(a1,a2,h1);
            mol.nho = vangle(a2,h1,ow);
        }
        else
        {
            mol.nh = v2;
            mol.cnh = vangle(a1,a2,h2);
            mol.nho = vangle(a2,h2,ow);
        }
        return 1;
    }
    return 0;
}

void vsub( const float *a1, const float *a2, const int &size, float *d )
{
    for (int i=0; i<size; i++) {
        d[i] = a2[i] - a1[i];
    }
}

float vlen_( const float *v )
{
    return pow( v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0.5);
}

float vlen(const float *a1, const float *a2 )
{
    float vec[3];
    vsub(a1,a2,3,vec);
    return vlen_(vec);
}

float ndot( const float *a1, const float *a2 )
{
    float rv1 = 1./vlen_(a1);
    float rv2 = 1./vlen_(a2);
    return a1[0]*rv1*a2[0]*rv2 + a1[1]*rv1*a2[1]*rv2 + a1[2]*rv1*a2[2]*rv2;
}

float vangle(const float *a1, const float *a2, const float *a3)
{
    float v21[3], v23[3];
    vsub(a2,a1,3,v21);
    vsub(a2,a3,3,v23);
    float d = ndot(v21,v23);
    if (d > 1) { d = 1; };
    return RAD2DEG(acos(d));
}


/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */

int gmx_nitrile_hbond(int argc, char *argv[])
{
    const char *desc[] = {
        "\tThis convert from gromacs to APBS pqr files.",
        "Use the -a1 flag for the nitrile carbon",
        "and the -a2 flag for the nitrile nitrogen.",
        "\n\tUse the -select flag to specify which atoms you",
        "want to parse through.  This can have a large effect",
        "on performance, and as such should be used.",
        "Using -select",
        "'resname SOL and same residue as within 0.5 of resname",
        "CNC and name NE' is recommended because it would loop over",
        "only waters near the cyanocysteine nitrogen.",
        "\n"
    };
    /* Command-line arguments */
    /*
     Hydrogen-bond acceptor properties of nitriles: a combined crystallographic and ab initio theoretical investigation
     Jean-Yves Le Questel, Michel Berthelot and Christian Laurence
     r average = .205(.020)
     theta1 = CNH, average = 145(23)
     theta2 = NHO, average = 156(18)
     95% of data within +/- 2 * STD
     */
    t_analysisdata      d;
    d.a1        = -1;
    d.a2        = -1;
    d.bVerbose  = false;
    d.cr        = 0.30;
    d.r         = 0.205 + (2*0.02);
    d.cnh       = 145 - (2*23);
    d.nho       = 156 - (2*18);
    d.nhb       = 0;
    
    t_pargs         pa[] = {
        { "-a1", TRUE, etINT,
            {&d.a1}, "Starting atom for bond vector--ie: CD in CNC"},
        { "-a2", TRUE, etINT,
            {&d.a2}, "Ending atom for bond vector--ie: NE in CNC"},
        { "-NHO_cutoff", FALSE, etREAL,
            {&d.nho}, "Minimum N-H-O bond angle to count as a hydorgen bond (degrees)"},
        { "-CNH_cutoff", FALSE, etREAL,
            {&d.cnh}, "Minimum N-H-O bond angle to count as a hydorgen bond (degrees)"},
        { "-angle_NH_cutoff", FALSE, etREAL,
            {&d.r}, "Minimum N-H bond distance to count as a hydrogen bond include angle criteria (nm)"},
        { "-distance_NH_cutoff", FALSE, etREAL,
            {&d.cr}, "Minimum N-H bond distance to count as a hydrogen bond with only distance criteria (nm)"},
        { "-v", FALSE, etBOOL,
            {&d.bVerbose}, "Be slightly more verbose"}
    };
    
    t_filenm        fnm[] = {
         { efXVG, "-o",  "frame_hb", ffWRITE }
        ,{ efXVG, "-op", "persistent", ffWRITE }
        ,{ efXVG, "-oa", "hb_count", ffWRITE }
        ,{ efXVG, "-or", "geometry", ffWRITE }
    };
    d.doLog = false, d.doPersistent = false, d.doGeo = false;
    
#define NFILE asize(fnm)
    
    gmx_ana_traj_t      *trj;
    t_topology          *top;
    output_env_t        oenv;
    int                 ngrps;
    gmx_ana_selection_t **sel;
    int                 g;
    int                 rc;
    
    CopyRight(stderr, argv[0]);
    
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nanagrps(trj, -1);
    parse_trjana_args(trj, &argc, argv, 0,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    gmx_ana_get_topology(trj, FALSE, &top, NULL);

    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);
    gmx_ana_init_coverfrac(trj, CFRAC_SOLIDANGLE);
    
    /* open xvg file */
    d.fp = NULL;
    d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Trajectory Hydrogen Bonds","Step", "Number of Water Molecules", oenv);
    const char * legend[] = {   "Frame Number",
                                "Nearby Water Molecules",
                                "Le Questel HBonded Water Molecules",
                                "Cho Sigma HBonded Water Molecules",
                                "Cho Pi HBonded Water Molecules"
                            };
    xvgr_legend(d.fp,asize(legend),legend,oenv);
    
    d.fpp = NULL;
    if (opt2bSet("-op", NFILE, fnm))
    {
        d.doPersistent = true;
        const char *flegend[] = {   "Number of Hydrogen Bonds",
                                    "Persistent Water Molecules",
                                    "Persistent Hydrogen Bonds",
                                    "Total Persistent Water Molecules",
                                    "Total Persistent Hydrogen Bonds"
                                };
        d.fpp = xvgropen(opt2fn("-op", NFILE, fnm), "Persistent Waters and Hydrogen Bonds", "Number of Frames", "Number Water Molecules", oenv);
        xvgr_legend(d.fpp,asize(flegend), flegend, oenv);
    }
    
    d.fpa = NULL;
    if (opt2bSet("-oa", NFILE, fnm))
    {
        d.doLog = true;
        const char *alegend[] = {   "N Hydrogen Bonds",
                                    "Number of Frames (Le Questel Hydrogen Bonds)",
                                    "Number of Frames (Cho Sigma Hydrogen Bonds)",
                                    "Number of Frames (Cho Pi Hydrogen Bonds)"
                                };
        d.fpa = xvgropen(opt2fn("-oa", NFILE, fnm), "Frames with types of hydrogen bonds", "Number of Hydrogen Bonds", "Number of Frames", oenv);
        xvgr_legend(d.fpa,asize(alegend),alegend,oenv);
    }
    d.fpg = NULL;
    if (opt2bSet("-or", NFILE, fnm))
    {
        d.doGeo = true;
        const char *glegend[] = {   "Frame",
                                    "Water Index",
                                    "NH Distance (Angstrom)",
                                    "CNH Angle (Degrees)",
                                    "NHO Angle (Degrees)"
                                };
        d.fpg = xvgropen(opt2fn("-or", NFILE, fnm), "Hydrogen Bonding Geometry", "Step", "Geometry", oenv);
        xvgr_legend(d.fpg,asize(glegend), glegend, oenv);
    }
    
    /* Make sure -a1 and -a2 are included and increment by -1 to match internal numbering */
    if ( d.a1<0 || d.a2<0 ) {
        gmx_fatal(FARGS, "\nAtom numbers -a1 and -a2 defining the bond vector must be specified\n");
    }
    d.a1--; d.a2--;
    
    /* Make sure that -a1 and -a2 exist in the trajectory */
    if (d.a1<0 || d.a1>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a1 is %d, which is outside of the range 0 to %d\n",d.a1+1,top->atoms.nr+1);
    }
    if (d.a2<0 || d.a2>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a2 is %d, which is outside of the range 0 to %d\n",d.a2+1,top->atoms.nr+1);
    }

    /* Parse through the frames */
    gmx_ana_do(trj, 0, &analyze_frame, &d);
    
    /* Now we analyze the water information */
    int maxhb = 5;
    maxhb++;
    std::vector<int> frames_with_lQ (maxhb,0);
    std::vector<int> frames_with_sigma (maxhb,0);
    std::vector<int> frames_with_pi (maxhb,0);
    for (int i=0; i<d.framen; i++) {
        /* Print off information about each frame */
        fprintf(d.fp,"%10i %10i %10i %10i %10i\n",i,d.water[i].size(), d.frame_lQ_nhb[i], d.frame_sigma_nhb[i], d.frame_pi_nhb[i]);
        /* Find information about the number of frames with different numbers of Hbonds */
        if (d.doLog)
        {
            for (int j=0; j<maxhb; j++)
            {
                if (d.frame_lQ_nhb[i] == j)
                {
                    frames_with_lQ[j] ++;
                }
                if (d.frame_sigma_nhb[i] == j)
                {
                    frames_with_sigma[j] ++;
                }
                if (d.frame_pi_nhb[i] == j)
                {
                    frames_with_pi[j] ++;
                }
            }
        }
    }
    ffclose(d.fp);
    /* Write information about the number of frames with different numbers of Hbonds */
    if (d.doLog)
    {
        for (int i=0;i<maxhb;i++)
        {
            fprintf(d.fpa,"%10i %10i %10i %10i\n",i,frames_with_lQ[i],frames_with_sigma[i],frames_with_pi[i]);
        }
        ffclose(d.fpa);
    }
    /* Get information about persistant water and hydrogen bonds */
    if (d.doPersistent) {
        std::vector<int> persistant_water (d.framen,0);
        std::vector<int> persistant_hb (d.framen,0);
        for (int i=0; i<d.framen; i++) {
            for (int j=0 ; j<(int) d.water[i].size() ; j++) {
                bool in_previous_frame = false;
                if (i>0) {
                    for (int k=0; k<d.water[i-1].size(); k++) {
                        if (d.water[i][j].resid == d.water[i-1][k].resid) {
                            in_previous_frame = true;
                            break;
                        }
                    }
                }
                if (! in_previous_frame) {
                    int k = i+1;
                    int npersist = 0;
                    int nhb = 0;
                    bool persistant = true;
                    while (persistant && k < d.framen) {
                        persistant = false;
                        for (int l=0 ; l<int(d.water[k].size()) ; l++) {
                            if (d.water[i][j].resid == d.water[k][l].resid) {
                                persistant = true;
                                npersist++;
                                if (d.water[i][j].is_hb && d.water[k][l].is_hb) {
                                    nhb++;
                                }
                                break;
                            }
                        }
                        k++;
                    }
                    persistant_water[npersist]++;
                    persistant_hb[nhb]++;
                }
            }
        }
        std::vector<int> total_water (d.framen,0);
        std::vector<int> total_hb (d.framen,0);
        for (int i=0; i<d.framen; i++) {
            for (int j=i; j<d.framen; j++) {
                total_water[i] += persistant_water[j];
                total_hb[i] += persistant_hb[j];
            }
        }
        /* Write to xvg file */
        for (int i=0; i<d.framen;i++){
            fprintf(d.fpp,"%10i %10i %10i %10i %10i\n",i,persistant_water[i],persistant_hb[i],total_water[i], total_hb[i]);
        }
        ffclose(d.fpp);
    }
}



/*! \brief
 * The main function.
 *
 * In Gromacs, most analysis programs are implemented such that the \p main
 * function is only a wrapper for a \p gmx_something function that does all
 * the work, and that convention is also followed here. */
int main(int argc, char *argv[])
{
    gmx_nitrile_hbond(argc, argv);
    return 0;
}
