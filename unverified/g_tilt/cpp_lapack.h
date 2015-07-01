
/* self writtein determinant of a 3x3 matrix */
float det3x3(float *A);

/* SGESVD prototype */
extern "C" {
    void sgesvd_( char* jobu, char* jobvt, int* m, int* n, float* a,
                 int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt,
                 float* work, int* lwork, int* info );
}
/* lapack single precision SVD function */
void sgesvd(char jobu, char jobvt, int m, int n, float *a,
            int lda, float *s, float *u, int ldu, float *vt, int ldvt,
            float *work, int lwork, int info);