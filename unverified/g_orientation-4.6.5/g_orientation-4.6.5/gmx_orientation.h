//
//  gmx_orientation.h
//  g_orientation-4.6.5
//
//  Created by Andrew Ritchie on 2/28/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#ifndef g_orientation_4_6_5_gmx_orientation_h
#define g_orientation_4_6_5_gmx_orientation_h

#include <stdio.h>
#include <config.h>
#include <gromacs/statutil.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/vec.h>
#include <gromacs/copyrite.h>
#include <gromacs/tpxio.h>
#include <gromacs/index.h>
#include <math.h>
#include "eigensolver.h"
#include "eigio.h"

#endif

typedef struct
{
    int a;
    float b;
} my_type;

int gmx_orientation(int argc, const char * argv[]);
void cross( rvec a, rvec b, rvec result );
void projection( rvec a, rvec b, rvec result );
double angle_periodicity( double angle_in );
double my_acos( double ip );
int parse_args();