//
//  my_structs.h
//  gmx2pqr-4.6.5
//
//  Created by Andrew Ritchie on 2/12/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#ifndef gmx2pqr_my_structs_h
#define gmx2pqr_my_structs_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#endif


#ifndef MPI
#define MPI 4*atan(1)
#endif

#define eps0 8.8541878e-12 //Vacuum permittivity C^2 s^2 /(kg m^3)
#define ec 1.6021773e-19 // Coulombs/e-
#define NA 6.0221367e23 // Avogadro's number, molecules/mol
#define kb 1.3806581e-23 // J/K



typedef struct gmx2amb gmx2amb;
typedef struct gmx2amb_dat gmx2amb_dat;
typedef struct hbond_array hbond_array;

struct gmx2amb {
    char resname[20];
    char gmxname[20];
    char ambername[20];
    float charge, radius;
};

struct gmx2amb_dat {
    gmx2amb *dat;
    int n;
};