#include "mex.h"
#include "vary_dynamics.h"
#include "lie_theory.h"

#include <math.h>
#include <gsl/gsl_matrix.h>

void printVector(const gsl_vector *v) {
    int size = v->size;

    printf("Vector:\n");

    int i;
    for (i=0;i<size;i++) {
        printf("%f\n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void printMatrix(const gsl_matrix *m) {
    int height = m->size1;
    int width = m->size2;

    printf("Matrix: \n");

    int i,j;
    for (i=0;i<height;i++) {
        for (j=0;j<width-1;j++) {
            printf("%f ", gsl_matrix_get(m,i,j));
        }
        printf("%f\n", gsl_matrix_get(m,i,width-1));
    }
    printf("\n");
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*
    Interface to the stepdynamics function for the variaitonal lie integrator of the beam

    for now just passing g, xi ,eta, mu, and lambda back and forth
    leaving the properties in the step dynamics function

    matlab is transposed

    */

    int i,j;

    int n; //assuming everything is the right size
    const mwSize *dim;
    dim = mxGetDimensions(prhs[0]);
    n = (int)dim[1];
    double *g_in;
    g_in = mxGetPr(prhs[0]);
    gsl_matrix *g = gsl_matrix_alloc(n,12);
    for (i=0; i<n; i++) {
        for (j=0; j<12; j++) {
            gsl_matrix_set(g, i, j, g_in[i*12+j]);
        }
    }



    double *xi_in;
    double *eta_in;
    double *mu_in;
    double *lambda_in;
    xi_in = mxGetPr(prhs[1]);
    eta_in = mxGetPr(prhs[2]);
    mu_in = mxGetPr(prhs[3]);
    lambda_in = mxGetPr(prhs[4]);
    gsl_matrix *xi = gsl_matrix_alloc(n,6);
    gsl_matrix *eta = gsl_matrix_alloc(n,6);
    gsl_matrix *mu = gsl_matrix_alloc(n,6);
    gsl_matrix *lambda = gsl_matrix_alloc(n,6);
    for (i=0; i<n; i++) {
        for (j=0; j<6; j++) {
            gsl_matrix_set(xi, i, j, xi_in[i*6+j]);
            gsl_matrix_set(eta, i, j, eta_in[i*6+j]);
            gsl_matrix_set(mu, i, j, mu_in[i*6+j]);
            gsl_matrix_set(lambda, i, j, lambda_in[i*6+j]);
        }
    }

    double *q_in;
    q_in = mxGetPr(prhs[5]);
    const mwSize *q_dim;
    q_dim = mxGetDimensions(prhs[5]);
    int N_act = (int)(q_dim[1]*q_dim[0]);
    gsl_vector *q = gsl_vector_alloc(N_act);
    for (i=0; i<N_act; i++) {
        gsl_vector_set(q, i, q_in[i]);
    }



    gsl_matrix *g_next = gsl_matrix_alloc(n,12);
    gsl_matrix *xi_next = gsl_matrix_alloc(n,6);
    gsl_matrix *eta_next = gsl_matrix_alloc(n,6);
    gsl_matrix *mu_next = gsl_matrix_alloc(n,6);
    gsl_matrix *lambda_next = gsl_matrix_alloc(n,6);

    double dt = 0.1*(10e-2/(n-1))/pow(100e3/1000,0.5);
    //double dt = 1e-5;

    stepDynamics(n, dt, q, g, xi, eta, mu, lambda, g_next, xi_next, eta_next, mu_next, lambda_next);


    //put it all back in matlab
    double *g_out_vec, *xi_out_vec, *eta_out_vec, *mu_out_vec, *lambda_out_vec;
    mxArray *g_out, *xi_out, *eta_out, *mu_out, *lambda_out;
    g_out = plhs[0] = mxCreateDoubleMatrix(12,n,mxREAL);
    g_out_vec = mxGetPr(g_out);
    xi_out = plhs[1] = mxCreateDoubleMatrix(6,n,mxREAL);
    xi_out_vec = mxGetPr(xi_out);
    eta_out = plhs[2] = mxCreateDoubleMatrix(6,n,mxREAL);
    eta_out_vec = mxGetPr(eta_out);
    mu_out = plhs[3] = mxCreateDoubleMatrix(6,n,mxREAL);
    mu_out_vec = mxGetPr(mu_out);
    lambda_out = plhs[4] = mxCreateDoubleMatrix(6,n,mxREAL);
    lambda_out_vec = mxGetPr(lambda_out);


    for (i=0;i<n;i++) {
        for (j=0; j<12; j++) {
            g_out_vec[i*12+j] = gsl_matrix_get(g_next,i,j);
        }
        for (j=0;j<6;j++) {
            xi_out_vec[i*6+j] = gsl_matrix_get(xi_next, i, j);
            eta_out_vec[i*6+j] = gsl_matrix_get(eta_next, i, j);
            mu_out_vec[i*6+j] = gsl_matrix_get(mu_next, i, j);
            lambda_out_vec[i*6+j] = gsl_matrix_get(lambda_next, i, j);
        }
    }

    gsl_vector_free(q);

    gsl_matrix_free(g);
    gsl_matrix_free(xi);
    gsl_matrix_free(eta);
    gsl_matrix_free(mu);
    gsl_matrix_free(lambda);
    gsl_matrix_free(g_next);
    gsl_matrix_free(xi_next);
    gsl_matrix_free(eta_next);
    gsl_matrix_free(mu_next);
    gsl_matrix_free(lambda_next);

}
