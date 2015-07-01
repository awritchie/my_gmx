#include <iostream>
#include <vector>
#include <cmath>
// my cpp blas and lapack headers
#include "cpp_blas.h"
#include "cpp_lapack.h"
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
#define M_PI = atan(1.)*4.
#endif

#ifndef RAD2DEG
#define RAD2DEG 180./M_PI
#endif

float det3x3(float * A);

typedef struct
{
    FILE   *fp;
    std::vector<float> time;
    std::vector<float> rotation;
    std::vector<float> distance;
    std::vector<float> rmsd;
    std::vector<float> x_rotation;
    std::vector<float> y_rotation;
    std::vector<float> z_rotation;
    std::vector<float> x_rotation_axis;
    std::vector<float> y_rotation_axis;
    std::vector<float> z_rotation_axis;
} t_tiltdata;

typedef std::vector< std::vector<float> > Matrix;

int gmx_tilt(int argc, char *argv[]);
void readCoords(int n_atoms, atom_id ind_atoms[], t_trxframe *fr, t_topology *top, 
                std::vector<float> &coords, gmx_bool bVerbose);
void readCoords(int n_atoms, atom_id ind_atoms[], rvec *x, t_topology *top, 
                std::vector<float> &coords, gmx_bool bVerbose);
void displacement( std::vector<float> reference, std::vector<float> frame, t_tiltdata &data, gmx_bool bVerbose );
void average_coordinate(std::vector<float> coords, std::vector<float> &xyz);
void kabsch_alignment( std::vector<float> ref, std::vector<float> tar, t_tiltdata &data, gmx_bool bVerbose);
