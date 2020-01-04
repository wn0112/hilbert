#include "stdlib.h"
#include "fftw3.h"
#include "math.h"

void hilbert_r2c(int size, double in[], fftw_complex out[]);
void hilbert_r2r(int size, double in[], double out[]);

const int N = 16;
int main()
{
    // create input array
    double x[N];
    for(int i = 0; i < N; i++)
        x[i] = i + 1;

    fftw_complex *hil_out_c;
    double *hil_out_r;

    hil_out_c = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    hil_out_r = (double *) malloc(sizeof(double) * N);

    hilbert_r2c(N, x, hil_out_c);
    hilbert_r2r(N, x, hil_out_r);

    printf("===================\n");
    printf("Complex output:\n\n");
    for(int i = 0; i < N; i++)
    {
        printf("%lf, %lf\n", hil_out_c[i][0], hil_out_c[i][1]);
    }
    printf("===================\n");
    printf("Real output:\n\n");
    for(int i = 0; i < N; i++)
    {
        printf("%lf\n", hil_out_r[i]);
    }
    fftw_free(hil_out_c);
    free(hil_out_r);
    return 0;
}

/* 
 * ********************************************
 * NOTE: 
 * fftw3 is required.
 * scipy is a python package to process digital signal, having a function hilbert(), this function has same output with scipy
 * ********************************************
 * :param: size     : size of in and out
 * :param: in       : double array
 * :param: out      : complex array
 */
void hilbert_r2c(int size, double in[], fftw_complex out[])
{
    int middle_pos;
    double *h;
    fftw_plan fft_plan;
    fftw_plan ifft_plan;
    fftw_complex *fft_out;

    h = (double *)calloc(size, sizeof(double));
    fft_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);

    // create h array
    h[0] = 1;
    if(size % 2 == 0)
    {
        h[size / 2] = 1;
        middle_pos = size / 2;
        for(int i = 1; i < middle_pos; i++)
            h[i] = 2;
    } else {
        middle_pos = (size + 1) / 2;
        for(int i = 1; i < middle_pos; i++)
            h[i] = 2;
    }

    // fft
    fft_plan = fftw_plan_dft_r2c_1d(size, in, fft_out, FFTW_ESTIMATE);
    fftw_execute(fft_plan);

    // multiplied by h array
    for(int i = 0; i < size; i++)
    {
        fft_out[i][0] *= h[i];
        fft_out[i][1] *= h[i];
    }

    // ifft
    ifft_plan = fftw_plan_dft_1d(size, fft_out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft_plan);

    // correct
    for(int i = 0; i < size; i++)
    {
        out[i][0] /= size;
        out[i][1] /= size;
    }

    // free memory
    fftw_destroy_plan(fft_plan);
    fftw_destroy_plan(ifft_plan);
    fftw_free(fft_out);
    free(h);
}

/*
 * :param: size     : size of in and out
 * :param: in       : double array
 * :param: out      : double array
 */
void hilbert_r2r(int size, double in[], double out[])
{
    int middle_pos;
    double *h;
    fftw_plan fft_plan;
    fftw_plan ifft_plan;
    fftw_complex *fft_out;
    fftw_complex *hil_out;

    h = (double *)calloc(size, sizeof(double));
    fft_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);
    hil_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size);

    // create h array
    h[0] = 1;
    if(size % 2 == 0)
    {
        h[size / 2] = 1;
        middle_pos = size / 2;
        for(int i = 1; i < middle_pos; i++)
            h[i] = 2;
    } else {
        middle_pos = (size + 1) / 2;
        for(int i = 1; i < middle_pos; i++)
            h[i] = 2;
    }

    // fft
    fft_plan = fftw_plan_dft_r2c_1d(size, in, fft_out, FFTW_ESTIMATE);
    fftw_execute(fft_plan);

    // multiplied by h array
    for(int i = 0; i < size; i++)
    {
        fft_out[i][0] *= h[i];
        fft_out[i][1] *= h[i];
    }

    // ifft
    ifft_plan = fftw_plan_dft_1d(size, fft_out, hil_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft_plan);

    // correct
    for(int i = 0; i < size; i++)
    {
        hil_out[i][0] /= size;
        hil_out[i][1] /= size;
        out[i] = sqrt(pow(hil_out[i][0], 2) + pow(hil_out[i][1], 2));
    }

    // free memory
    fftw_destroy_plan(fft_plan);
    fftw_destroy_plan(ifft_plan);
    fftw_free(fft_out);
    fftw_free(hil_out);
    free(h);
}