/* blas_shim4.c -- route TROVE's dsyevr_/dgemm_ (LP64, gfortran) to the multithreaded ILP64 OpenBLAS
 * shipped with numpy (libscipy_openblas64_*.so, symbols scipy_*_64_).
 *
 * The OpenBLAS library is NOT linked: it is dlopen'ed with RTLD_LOCAL|RTLD_DEEPBIND, because it exports a
 * few unprefixed ILP64 LAPACK auxiliaries (dlarf1f_, droundup_lwork_, ...) that would otherwise capture the
 * calls made by netlib LAPACK 3.12 (e.g. dsyev -> dorgtr -> dlarf1f) and crash.  Everything else in TROVE
 * keeps using netlib LAPACK/BLAS.  Calls made from inside OpenMP parallel regions, or when the library
 * cannot be loaded, go to netlib (dlsym(RTLD_NEXT)).
 * Library path: $TROVE_OPENBLAS64, default below.  Threads: $OPENBLAS_NUM_THREADS. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
typedef int32_t i32; typedef int64_t i64;

#ifndef OB64_DEFAULT
#define OB64_DEFAULT "libscipy_openblas64_.so"  /* set TROVE_OPENBLAS64 to the full path of numpy's scipy-openblas64 library */
#endif

typedef void (*dsyevr32_t)(const char*, const char*, const char*, const i32*, double*, const i32*,
  const double*, const double*, const i32*, const i32*, const double*, i32*, double*, double*, const i32*,
  i32*, double*, const i32*, i32*, const i32*, i32*, size_t, size_t, size_t);
typedef void (*dgemm32_t)(const char*, const char*, const i32*, const i32*, const i32*, const double*,
  const double*, const i32*, const double*, const i32*, const double*, double*, const i32*, size_t, size_t);
typedef void (*dsyevr64_t)(const char*, const char*, const char*, const i64*, double*, const i64*,
  const double*, const double*, const i64*, const i64*, const double*, i64*, double*, double*, const i64*,
  i64*, double*, const i64*, i64*, const i64*, i64*, size_t, size_t, size_t);
typedef void (*dgemm64_t)(const char*, const char*, const i64*, const i64*, const i64*, const double*,
  const double*, const i64*, const double*, const i64*, const double*, double*, const i64*, size_t, size_t);

static dsyevr32_t netlib_dsyevr = NULL;
static dgemm32_t netlib_dgemm = NULL;
static dsyevr64_t ob_dsyevr = NULL;
static dgemm64_t ob_dgemm = NULL;

__attribute__((constructor)) static void shim_init(void) {
  netlib_dsyevr = (dsyevr32_t)dlsym(RTLD_NEXT, "dsyevr_");
  netlib_dgemm = (dgemm32_t)dlsym(RTLD_NEXT, "dgemm_");
  const char *path = getenv("TROVE_OPENBLAS64");
  if (!path || !*path) path = OB64_DEFAULT;
  void *h = dlopen(path, RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
  if (h) {
    ob_dsyevr = (dsyevr64_t)dlsym(h, "scipy_dsyevr_64_");
    ob_dgemm = (dgemm64_t)dlsym(h, "scipy_dgemm_64_");
  }
  if (!h || !ob_dsyevr || !ob_dgemm) {
    fprintf(stderr, "blas_shim4: OpenBLAS64 not available (%s); netlib dsyevr/dgemm will be used\n",
            h ? "missing symbols" : dlerror());
    ob_dsyevr = NULL; ob_dgemm = NULL;
  }
}

void dsyevr_(const char* jobz, const char* range, const char* uplo, const i32* n, double* a, const i32* lda,
  const double* vl, const double* vu, const i32* il, const i32* iu, const double* abstol, i32* m,
  double* w, double* z, const i32* ldz, i32* isuppz, double* work, const i32* lwork, i32* iwork,
  const i32* liwork, i32* info, size_t l1, size_t l2, size_t l3)
{
  if (!ob_dsyevr || omp_in_parallel()) {
    netlib_dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
                  iwork, liwork, info, l1, l2, l3);
    return;
  }
  i64 N = *n, LDA = *lda, IL = *il, IU = *iu, M = 0, LDZ = *ldz, LWORK = *lwork, LIWORK = *liwork, INFO = 0;
  i64 nsup = 2 * (N > 0 ? N : 1);
  i64 liw = (LIWORK > 0) ? LIWORK : 1;
  i64 *isup = (i64*)calloc((size_t)nsup, sizeof(i64));
  i64 *iw = (i64*)calloc((size_t)liw, sizeof(i64));
  ob_dsyevr(jobz, range, uplo, &N, a, &LDA, vl, vu, &IL, &IU, abstol, &M, w, z, &LDZ,
            isup, work, &LWORK, iw, &LIWORK, &INFO, 1, 1, 1);
  *info = (i32)INFO;
  if (LWORK == -1 || LIWORK == -1) {
    iwork[0] = (i32)iw[0];               /* workspace query: optimal liwork */
  } else {
    *m = (i32)M;                         /* isuppz is not used by TROVE: not copied back */
  }
  free(isup); free(iw);
}

void dgemm_(const char* ta, const char* tb, const i32* m, const i32* n, const i32* k, const double* alpha,
  const double* a, const i32* lda, const double* b, const i32* ldb, const double* beta, double* c,
  const i32* ldc, size_t lta, size_t ltb)
{
  if (!ob_dgemm || omp_in_parallel()) {
    netlib_dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, lta, ltb);
    return;
  }
  i64 M = *m, N = *n, K = *k, LDA = *lda, LDB = *ldb, LDC = *ldc;
  ob_dgemm(ta, tb, &M, &N, &K, alpha, a, &LDA, b, &LDB, beta, c, &LDC, 1, 1);
}
