#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <cstring>
// gromacs c headers
#ifdef __cplusplus
extern "C"  
#endif
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
#include <gromacs/tpxio.h>

#ifndef M_PI
#define M_PI atan(1)
#endif

struct t_hbdata
{
    FILE   *fp;
    FILE   *fpp;
    const char   *fpa;
    std::vector<int> nps; // Number of persistent waters per water molecule
    std::vector<int> nhbs; // Number of persistent h-bonds per water molecule
    int nhb;
};

struct t_mol
{
    int resid;
    int rstatus;
};

struct parms
{
    real r;
    real cnh;
    real nho;
    real rmax;
};

typedef std::vector<t_mol> t_water;

typedef std::vector< std::vector<double> > Matrix;

int gmx_CNCHOH(int argc, char *argv[]);

t_water analyze_frame(t_trxframe *fr, t_topology *top, parms cutoffs);

int parse_water(std::vector<double> CD, std::vector<double> NE, std::vector<double> OW, std::vector<double> HW1, std::vector<double> HW2, parms cutoffs);

int is_Hbond(std::vector<double> CD, std::vector<double> NE, std::vector<double> OW, std::vector<double> HW, parms cutoffs);

double vecangle(std::vector<double> r1, std::vector<double> r2, std::vector<double> r3);

double veclen(std::vector<double> r1, std::vector<double> r2);

std::vector<double> vecsub(std::vector<double> r1, std::vector<double> r2);

