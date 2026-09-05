/*
 * fft_test.c — correctness + OpenMP thread-safety test for fft.c
 *
 * Compile:
 *   gcc -O2 -Wall -DUSE_ISO_C_BINDINGS -fopenmp -I$ELMER_BUILD/fem/src \
 *       fft_test.c $ELMER_SRC/fem/src/fft.c \
 *       -lm -o fft_test
 * Run:
 *   OMP_NUM_THREADS=8 ./fft_test
 *
 * Note the include path is $ELMER_BUILD/fem/src, not $ELMER_BUILD/fem: fft.c
 * includes "../config.h", which is resolved relative to the -I directory, so
 * pointing at fem looks for $ELMER_BUILD/config.h and fails.
 *
 * -DUSE_ISO_C_BINDINGS is needed as well; the real build passes it on the
 * command line, and without it config.h leaves FC_FUNC undefined.
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

/* independent brute-force DFT reference, deliberately not sharing any code
   with fft.c:  R[n] = sum(T[k]*exp(-i*2*pi*n*k/N)) */
static void refdft(int N, COMPLEX *T, COMPLEX *R)
{
    double pi = acos(-1.0);
    for (int n = 0; n < N; n++) {
        long double sr = 0.0L, si = 0.0L;
        for (int k = 0; k < N; k++) {
            long double a = -2.0L * (long double)pi * (long double)n * k / N;
            long double c = cosl(a), s = sinl(a);
            sr += (long double)T[k].Real * c - (long double)T[k].Imag * s;
            si += (long double)T[k].Real * s + (long double)T[k].Imag * c;
        }
        R[n].Real = (double)sr;
        R[n].Imag = (double)si;
    }
}

/* ---- serial correctness tests ------------------------------------------ */

/* cfftf against the brute-force reference, for any length */
static void test_vs_reference(int N)
{
    double pi = acos(-1.0);
    COMPLEX *T = malloc(N * sizeof(COMPLEX));
    COMPLEX *F = malloc(N * sizeof(COMPLEX));
    COMPLEX *R = malloc(N * sizeof(COMPLEX));
    for (int k = 0; k < N; k++) {
        T[k].Real = cos(2.0*pi*3*k/N) + 0.25*k/N;
        T[k].Imag = sin(2.0*pi*5*k/N) - 0.75;
    }
    refdft(N, T, R);
    cfftf(N, T, F);
    double e = cmaxerr(F, R, N);
    char name[64]; snprintf(name, sizeof(name), "cfftf vs reference N=%d", N);
    check(name, e, 1e-9 * N);
    free(T); free(F); free(R);
}

/* same, computed in place (F aliases T) */
static void test_in_place(int N)
{
    double pi = acos(-1.0);
    COMPLEX *T = malloc(N * sizeof(COMPLEX));
    COMPLEX *R = malloc(N * sizeof(COMPLEX));
    for (int k = 0; k < N; k++) {
        T[k].Real = cos(2.0*pi*2*k/N);
        T[k].Imag = 0.5 - sin(2.0*pi*k/N);
    }
    refdft(N, T, R);
    cfftf(N, T, T);
    double e = cmaxerr(T, R, N);
    char name[64]; snprintf(name, sizeof(name), "cfftf in place N=%d", N);
    check(name, e, 1e-9 * N);
    free(T); free(R);
}


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

/* rfftb against a brute-force inverse of the Hermitian spectrum F[0..N/2].
   tested directly rather than only through an rfftf/rfftb roundtrip, so that a
   pair of compensating errors in the two cannot hide. */
static void test_rfftb_vs_reference(int N)
{
    int M = N/2;
    double pi = acos(-1.0);
    COMPLEX *F = malloc((M + 1) * sizeof(COMPLEX));
    double  *T = malloc((N + 2) * sizeof(double));
    double  *R = malloc(N * sizeof(double));

    /* dense spectrum, nothing zero except the imaginary parts at the two ends
       (which must vanish for the result to be real) */
    for (int k = 0; k <= M; k++) {
        F[k].Real = cos(1.7*k) + 0.3*k/M;
        F[k].Imag = sin(2.3*k) - 0.4;
    }
    F[0].Imag = 0.0;
    F[M].Imag = 0.0;

    for (int j = 0; j < N; j++) {
        long double s = 0.0L;
        for (int n = 0; n < N; n++) {
            long double re = (n <= M) ?  F[n].Real :  F[N-n].Real;
            long double im = (n <= M) ?  F[n].Imag : -F[N-n].Imag;
            long double a  = 2.0L * (long double)pi * j * n / N;
            s += re*cosl(a) - im*sinl(a);
        }
        R[j] = (double)s;
    }

    rfftb(N, F, T);

    double e = rmaxerr(T, R, N);
    char name[64]; snprintf(name, sizeof(name), "rfftb vs reference N=%d", N);
    check(name, e, 1e-9 * N);
    free(F); free(T); free(R);
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

    printf("\n--- arbitrary lengths (vs brute-force DFT) ---\n");
    /* trivial, small, prime, prime power, highly composite, mixed */
    static const int lens[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 17,
                                24, 25, 31, 36, 60, 64, 100, 121, 128, 210,
                                243, 256, 360, 512, 1000, 1009, 1024, 1331 };
    for (unsigned i = 0; i < sizeof(lens)/sizeof(lens[0]); i++)
        test_vs_reference(lens[i]);

    printf("\n--- arbitrary lengths, in place ---\n");
    for (unsigned i = 0; i < sizeof(lens)/sizeof(lens[0]); i++)
        test_in_place(lens[i]);

    printf("\n--- arbitrary lengths, cfftf/cfftb roundtrip ---\n");
    for (unsigned i = 0; i < sizeof(lens)/sizeof(lens[0]); i++)
        test_roundtrip_cfft(lens[i]);

    printf("\n--- arbitrary lengths, impulse & cosine peak ---\n");
    test_impulse(3);   test_impulse(12);  test_impulse(1009);
    test_sine_peak(6, 1); test_sine_peak(100, 7); test_sine_peak(1009, 31);

    printf("\n--- real transforms, non power of two ---\n");
    /* rfftf and rfftb need an even length, nothing more. the N%4==2 entries
       exercise the odd half length (N/2) path in rfftb, where the boundary
       term the even case handles separately does not exist. */
    static const int rlens[] = { 2, 4, 6, 8, 10, 12, 14, 18, 20, 24, 30, 36,
                                 50, 66, 100, 126, 250, 360, 1000, 1002 };
    for (unsigned i = 0; i < sizeof(rlens)/sizeof(rlens[0]); i++)
        test_roundtrip_rfft(rlens[i]);

    printf("\n--- rfftb vs brute-force Hermitian inverse ---\n");
    for (unsigned i = 0; i < sizeof(rlens)/sizeof(rlens[0]); i++)
        test_rfftb_vs_reference(rlens[i]);

    printf("\n--- OpenMP thread safety ---\n");
    test_openmp_parallel(1024, 1000);
    test_openmp_rfft_parallel(1024, 1000);
    test_openmp_parallel(4096, 200);
    test_openmp_rfft_parallel(4096, 200);

    printf("\n%d passed, %d failed\n", pass, fail);
    return fail ? 1 : 0;
}
