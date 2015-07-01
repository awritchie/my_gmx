
/* SGEMV prototype */
extern "C" {
    void sgemv_(const char *trans, const int *m, const int *n, const float *alpha, const float *a, const int *lda,
                const float *x, const int *incx, const float *beta, const float *y, const int *incy);
}
/* blas single precision matrix-vector multiplication prototype */
void sgemv(char trans, int m, int n, float alpha, float *a, int lda,
           float *x, int incx, float beta, float *y, int incy);

/* DDOT prototype */
extern "C" {
    double ddot_( int *n, double *x, int *incx, double *y, int *incy );
}
/* blas double precision dot product */
double ddot(int n, double *x, int incx, double *y, int incy);

/* SDOT prototype */
extern "C" {
    double sdot_( int *n, float *x, int *incx, float *y, int *incy );
}
/* blas single precision dot product */
double sdot(int n, float *x, int incx, float *y, int incy);

/* SGEMM prototype */
extern "C" {
    void sgemm_(char* transa, char* transb, int* m, int* n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
}
/* blas single prevision matrix matrix multiplication */
void sgemm(char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc);