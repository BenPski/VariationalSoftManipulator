#ifndef VARY_DYNAMICS_H
#define VARY_DYNAMICS_H

#include <gsl/gsl_matrix.h>

void stepDynamics(int n, double dt, const gsl_vector *q, const gsl_matrix *g, const gsl_matrix *xi, const gsl_matrix *eta, const gsl_matrix *mu, const gsl_matrix *lambda, gsl_matrix *g_next, gsl_matrix *xi_next, gsl_matrix *eta_next, gsl_matrix *mu_next, gsl_matrix *lambda_next);

#endif
