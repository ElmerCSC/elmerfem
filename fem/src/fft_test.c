/*
 * fft_test.c — correctness + OpenMP thread-safety test for fft.c
 *
 * Compile:
 *   gcc -O2 -fopenmp -I$ELMER_BUILD/fem \
 *       fft_test.c $ELMER_SRC/fem/src/fft.c \
 *       -lm -o fft_test
 * Run:
 *   OMP_NUM_THREADS=8 ./fft_test
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

typedef struct { double Real; double Imag; } COMPLEX;

void cfftf(int N, COMPLEX *T, COMPLEX *F);
void cfftb(int N, COMPLEX *F, COMPLEX *T);
void rfftf(int N, double *T, COMPLEX *F);
void rfftb(int N, COMPLEX *F, double *T);

static int  pass = 0, fail = 0;
static void check(const char *name, double err, double tol)
{
    if (err <= tol) {
        printf("  PASS  %-40s  err=%.2e\n", name, err);
        pass++;
    } else {
        printf("  FAIL  %-40s  err=%.2e  (tol %.2e)\n", name, err, tol);
        fail++;
    }
}

/* max absolute error between two complex arrays */
static double cmaxerr(COMPLEX *a, COMPLEX *b, int N)
{
    double e = 0.0;
    for (int k = 0; k < N; k++) {
        double dr = a[k].Real - b[k].Real;
        double di = a[k].Imag - b[k].Imag;
        double d  = sqrt(dr*dr + di*di);
        if (d > e) e = d;
    }
    return e;
}

/* max absolute error between two real arrays */
static double rmaxerr(double *a, double *b, int N)
{
    double e = 0.0;
    for (int k = 0; k < N; k++) {
        double d = fabs(a[k] - b[k]);
        if (d > e) e = d;
    }
    return e;
}

/* ---- serial correctness tests ------------------------------------------ */

static void test_impulse(int N)
{
    /* delta at t=0  ->  F[n] = 1 for all n */
    COMPLEX *T = calloc(N, sizeof(COMPLEX));
    COMPLEX *F = calloc(N, sizeof(COMPLEX));
    T[0].Real = 1.0;
    cfftf(N, T, F);
    double e = 0.0;
    for (int k = 0; k < N; k++) {
        double dr = F[k].Real - 1.0, di = F[k].Imag;
        double d  = sqrt(dr*dr + di*di);
        if (d > e) e = d;
    }
    char name[64]; snprintf(name, sizeof(name), "cfftf impulse N=%d", N);
    check(name, e, 1e-10);
    free(T); free(F);
}

static void test_roundtrip_cfft(int N)
{
    /* forward then inverse should recover original (scaled by N) */
    double pi = acos(-1.0);
    COMPLEX *in  = malloc(N * sizeof(COMPLEX));
    COMPLEX *F   = malloc(N * sizeof(COMPLEX));
    COMPLEX *out = malloc(N * sizeof(COMPLEX));
    for (int k = 0; k < N; k++) {
        in[k].Real = cos(2.0*pi*3*k/N) + 0.5*cos(2.0*pi*7*k/N);
        in[k].Imag = 0.0;
    }
    cfftf(N, in, F);
    cfftb(N, F, out);
    /* cfftb returns N * original */
    COMPLEX *scaled = malloc(N * sizeof(COMPLEX));
    for (int k = 0; k < N; k++) {
        scaled[k].Real = in[k].Real * N;
        scaled[k].Imag = in[k].Imag * N;
    }
    double e = cmaxerr(out, scaled, N);
    char name[64]; snprintf(name, sizeof(name), "cfftf/cfftb roundtrip N=%d", N);
    check(name, e, 1e-8 * N);
    free(in); free(F); free(out); free(scaled);
}

static void test_sine_peak(int N, int freq)
{
    /* real cosine at bin 'freq'  ->  F[freq].Real ≈ N/2, F[N-freq].Real ≈ N/2 */
    double pi = acos(-1.0);
    COMPLEX *T = calloc(N, sizeof(COMPLEX));
    COMPLEX *F = calloc(N, sizeof(COMPLEX));
    for (int k = 0; k < N; k++)
        T[k].Real = cos(2.0*pi*freq*k/N);
    cfftf(N, T, F);
    double expected = N / 2.0;
    double e = fabs(F[freq].Real - expected);
    e = fmax(e, fabs(F[N - freq].Real - expected));
    e = fmax(e, fabs(F[freq].Imag));
    char name[64]; snprintf(name, sizeof(name), "cfftf cosine peak N=%d f=%d", N, freq);
    check(name, e, 1e-8 * N);
    free(T); free(F);
}

static void test_roundtrip_rfft(int N)
{
    double pi = acos(-1.0);
    double *in  = malloc(N * sizeof(double));
    COMPLEX *F  = malloc((N/2 + 1) * sizeof(COMPLEX));
    double *out = malloc(N * sizeof(double));
    for (int k = 0; k < N; k++)
        in[k] = cos(2.0*pi*5*k/N) - 0.3*cos(2.0*pi*11*k/N);
    rfftf(N, in, F);
    rfftb(N, F, out);
    /* rfftb returns N * original */
    double *scaled = malloc(N * sizeof(double));
    for (int k = 0; k < N; k++) scaled[k] = in[k] * N;
    double e = rmaxerr(out, scaled, N);
    char name[64]; snprintf(name, sizeof(name), "rfftf/rfftb roundtrip N=%d", N);
    check(name, e, 1e-8 * N);
    free(in); free(F); free(out); free(scaled);
}

/* ---- OpenMP thread-safety test ----------------------------------------- */

static void test_openmp_parallel(int N, int nreps)
{
    /* each thread independently computes the same FFT and checks against serial */
    double pi = acos(-1.0);

    /* compute serial reference */
    COMPLEX *ref_in = malloc(N * sizeof(COMPLEX));
    COMPLEX *ref_F  = malloc(N * sizeof(COMPLEX));
    for (int k = 0; k < N; k++) {
        ref_in[k].Real = cos(2.0*pi*5*k/N) + sin(2.0*pi*13*k/N);
        ref_in[k].Imag = cos(2.0*pi*2*k/N);
    }
    cfftf(N, ref_in, ref_F);

    int nerrors = 0;

    #pragma omp parallel for reduction(+:nerrors) schedule(dynamic)
    for (int t = 0; t < nreps; t++) {
        COMPLEX *in = malloc(N * sizeof(COMPLEX));
        COMPLEX *F  = malloc(N * sizeof(COMPLEX));
        for (int k = 0; k < N; k++) {
            in[k].Real = cos(2.0*pi*5*k/N) + sin(2.0*pi*13*k/N);
            in[k].Imag = cos(2.0*pi*2*k/N);
        }
        cfftf(N, in, F);
        double e = cmaxerr(F, ref_F, N);
        if (e > 1e-8 * N) nerrors++;
        free(in); free(F);
    }

    char name[80];
    snprintf(name, sizeof(name), "OpenMP %d reps cfftf N=%d", nreps, N);
    check(name, (double)nerrors, 0.0);

    free(ref_in); free(ref_F);
}

static void test_openmp_rfft_parallel(int N, int nreps)
{
    double pi = acos(-1.0);
    double *ref_in = malloc(N * sizeof(double));
    COMPLEX *ref_F = malloc((N/2+1) * sizeof(COMPLEX));
    for (int k = 0; k < N; k++)
        ref_in[k] = sin(2.0*pi*7*k/N) - cos(2.0*pi*17*k/N);
    rfftf(N, ref_in, ref_F);

    int nerrors = 0;

    #pragma omp parallel for reduction(+:nerrors) schedule(dynamic)
    for (int t = 0; t < nreps; t++) {
        double  *in = malloc(N * sizeof(double));
        COMPLEX *F  = malloc((N/2+1) * sizeof(COMPLEX));
        for (int k = 0; k < N; k++)
            in[k] = sin(2.0*pi*7*k/N) - cos(2.0*pi*17*k/N);
        rfftf(N, in, F);
        double e = cmaxerr(F, ref_F, N/2+1);
        if (e > 1e-8 * N) nerrors++;
        free(in); free(F);
    }

    char name[80];
    snprintf(name, sizeof(name), "OpenMP %d reps rfftf N=%d", nreps, N);
    check(name, (double)nerrors, 0.0);

    free(ref_in); free(ref_F);
}

/* ---- main --------------------------------------------------------------- */

int main(void)
{
    printf("fft_test  (OMP_NUM_THREADS=%d)\n\n", omp_get_max_threads());

    printf("--- serial correctness ---\n");
    test_impulse(16);
    test_impulse(256);
    test_impulse(1024);
    test_sine_peak(256, 7);
    test_sine_peak(1024, 31);
    test_roundtrip_cfft(64);
    test_roundtrip_cfft(1024);
    test_roundtrip_rfft(64);
    test_roundtrip_rfft(1024);

    printf("\n--- OpenMP thread safety ---\n");
    test_openmp_parallel(1024, 1000);
    test_openmp_rfft_parallel(1024, 1000);
    test_openmp_parallel(4096, 200);
    test_openmp_rfft_parallel(4096, 200);

    printf("\n%d passed, %d failed\n", pass, fail);
    return fail ? 1 : 0;
}
