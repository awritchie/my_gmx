#ifndef gmx_nitrile_hbond_hpp
#define gmx_nitrile_hbond_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>


#ifdef __cplusplus
extern "C"
#endif
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/trajana.h>

#ifndef MPI
#define MPI (4.*atan(1.))
#endif

#define RAD2DEG(X) (X*180./MPI)
#define DEG2RAD(X) (X*MPI/180.)

struct t_mol
{
    int resid;
    int is_hb;
    float nh;
    float cnh;
    float nho;
};

typedef std::vector<t_mol> t_water;

/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    gmx_ana_selection_t                 *refsel;
    FILE                                *fp;
    FILE                                *fpp;
    FILE                                *fpa;
    FILE                                *fpg;
    int                                 framen;
    int                                 a1, a2;
    int                                 nhb;
    int                                 nwater;
    std::vector<int>                    p_water; // numer of persistent waters per water molecule
    std::vector<int>                    p_hbonds; // number of persisten hbonds per water molecule
    std::vector<t_water>                water;
    std::vector<int>                    frame_sigma_nhb;
    std::vector<int>                    frame_pi_nhb;
    std::vector<int>                    frame_lQ_nhb;
    float                               cr;
    float                               r;
    float                               cnh;
    float                               nho;
    bool                                bVerbose;
    bool                                doPersistent, doLog, doGeo;
} t_analysisdata;

int parse_water(const float *a1, const float *a2, const float *ow, const float *h1, const float *h2, const float &r, t_mol &mol);

void vsub( const float *a1, const float *a2, const int &size, float *d );

float vlen_( const float *v );

float vlen(const float *a1, const float *a2 );

float vangle(const float *a1, const float *a2, const float *a3);

float ndot( const float *a1, const float *a2 );

#endif