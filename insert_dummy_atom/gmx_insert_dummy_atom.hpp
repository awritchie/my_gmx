#include <iostream>
#include <vector>
#include <cmath>
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

int gmx_insert_dummy_atom(int argc, char *argv[]);
