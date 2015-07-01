#include "cpp_lapack.h"

/* self writtein determinant of a 3x3 matrix */
float det3x3(float *A)
{
    return A[0]*(A[4]*A[8]-A[7]*A[5])-A[3]*(A[1]*A[8]-A[7]*A[2])+A[6]*(A[1]*A[5]-A[4]*A[2]);
}

/* lapack single precision SVD function */
void sgesvd(char jobu, char jobvt, int m, int n, float *a,
            int lda, float *s, float *u, int ldu, float *vt, int ldvt,
            float *work, int lwork, int info)
{
    return sgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);
}