/*
 * $Id: template.c,v 1.5 2008/05/29 08:36:53 hess Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: template.c,v 1.5 2008/05/29 08:36:53 hess Exp $";

#include <math.h>
#include "statutil.h"
#include "vec.h"
#include "copyrite.h"
#include "index.h"
#include "smalloc.h"
#include "tpxio.h"
#include "gmx_fatal.h"
#include "xvgr.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute the dot product of a bond vector (atoms given in ",
    "index.ndx) and a bond vector (option -vec). If a projection ",
    "vector (-proj) is given, then projection is subtracted from the ",
    "bond vector before the dot product is computed."
  };
  
  static rvec refvector={0,0,0};
  static rvec projvector={0,0,0}; // default value
  rvec bondvector;
  int a1=-1, a2=-1;
  float cosine, angle;
  double dotproduct; // dot product of two vectors (for proj sub)
  rvec projection; // projection onto a vector (for proj sub)
  FILE *outputfile;

  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-vec", TRUE, etRVEC, {refvector}, "Reference vector for computing the dot product" },
    { "-proj", TRUE, etRVEC, {projvector}, "Subtract the projection onto this vector before computing the dot product" },
    { "-a1", TRUE, etINT, {&a1}, "Atom number 1 (must be set)" },
    { "-a2", TRUE, etINT, {&a2}, "Atom number 2 (must be set)" }
  };
 
  int i;
  t_trxframe fr;
  int        status;
  int        flags = TRX_READ_X;

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },    
    { efXVG, NULL, "dot", ffOPTWR }
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* Check for xtc and initialize first frame */
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  /* check for required parameters */
  if( a1<1 || a2<1 ) 
    gmx_fatal(FARGS, "Atom numbers a1 and a2 defining the bond vector must be specified\n" );

  /* internally, atom numbering starts from 0 */
  a1--;
  a2--;
 
  if( refvector[0]==0 && refvector[1]==0 && refvector[2]==0 )
    gmx_fatal(FARGS, "Trying to compute dot product with zero vector. Please set -vec.\n" );

  /* output file */
  outputfile=ffopen( opt2fn("-o",NFILE,fnm),"w" );
  if(outputfile==NULL) {
    gmx_fatal(FARGS,"Failed to open output file '%s'\n", opt2fn("-o",NFILE,fnm) );
  }

  /* some setup */
  printf( "Computing dot product of bond vector with " );
  for(i=0; i<3; i++ ) printf( "% 5.5f ", refvector[i] );
  printf("\n");
  if( norm(refvector) != 1 ){
    printf( "Rescaling reference vector, " );
    unitv( refvector, refvector );
    printf( "new vector: " );
    for(i=0; i<3; i++) printf("% 5.5f ", refvector[i] );
    printf("\n");
  }

  if( projvector[0] != 0 || projvector[1] != 0 || projvector[2] != 0 ){
    printf( "Subtracting projection onto " );
    for (i=0; i<3; i++ ) printf( "% 5.5f ", projvector[i] );
    printf( "\n" );
    
    if( norm(projvector) != 1 ){
      printf( "Rescaling projection vector, " );
      unitv( projvector, projvector );
      printf( "new vector: " );
      for(i=0; i<3; i++) printf("% 5.5f ", projvector[i] );
      printf("\n");
    }
  }
    

  /* loop over frames */
  do {
    /* coordinates are available in the vector fr.x */
    /*
    printf( "atom 1: " );
    for(i=0; i<3; i++) printf("% 5.5f ", fr.x[a1][i] );
    printf( "\natom 2: " );
    for(i=0; i<3; i++) printf("% 5.5f ", fr.x[a2][i] );
    printf("\n");
    */

    rvec_sub(fr.x[a2],fr.x[a1],bondvector);

    /* project. If projvector is the zero vector then this leaves
       bondvector unchanged */
    dotproduct=iprod(bondvector,projvector);
    svmul(dotproduct, projvector, projection );
    rvec_sub(bondvector, projection, bondvector); 

    /*
    printf( "vector " );
    for(i=0; i<3; i++) printf("% 5.5f ", bondvector[i] );
    printf("\n");
    */
    unitv( bondvector, bondvector );

    /*
    printf( "vector: " );
    for(i=0; i<3; i++) printf("% 5.5f ", bondvector[i] );
    printf("\n");
    */
    /* answer to display is the dot product of 'bondvector' with 'refvector' */
    cosine=iprod( bondvector, refvector ); // /norm(bondvector);
    if(cosine>=1) { 
      angle=0.0;
    } else if(cosine<=-1) {
      angle=180.0;
    } else {
      angle=acos(cosine)*180/M_PI;
    }

    /* done with frame */
    fprintf( outputfile, "%-10.4f\t% 1.4f\t% 1.4f\n", fr.time, cosine, angle );
  } while(read_next_frame(status,&fr));

  fclose( outputfile );
  thanx(stderr);
  return 0;
}

