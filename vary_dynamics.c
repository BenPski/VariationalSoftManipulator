/*

The variational dynamics integrator for the rod

The basic algorithm is to compute the updates for g, xi, eta, mu, and lambda (config, strain, velocity, momentum, and stress) at each step
    g, xi, mu, and lambda are explicit and eta is implicit which can be solved quickly (usually two steps of newton raphson)

    it is not necessary to compute mu and lambda, but it reduces the number of repeated computations


Can use various different diffeomorpisms for this and an approximation to the differential exponential
    for now using exp rather than cay for testing
    use only the first approximation to the differential (I + ad_x/2)

*/

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "lie_theory.h"


void printVectorHelp(const gsl_vector *v) {
    int size = v->size;

    printf("Vector:\n");

    int i;
    for (i=0;i<size;i++) {
        printf("%f\n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void printMatrixHelp(const gsl_matrix *m) {
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

/*
struct IntegratorParams {
    int n; //discretizations
    double dt; //time step
    double ds; //space step

    gsl_matrix *g;
    gsl_matrix *xi;
    gsl_matrix *eta;
    gsl_matrix *mu;
    gsl_matrix *lambda;
}
*/
void stepDynamics(int n, double dt, const gsl_vector *q, const gsl_matrix *g, const gsl_matrix *xi, const gsl_matrix *eta, const gsl_matrix *mu, const gsl_matrix *lambda, gsl_matrix *g_next, gsl_matrix *xi_next, gsl_matrix *eta_next, gsl_matrix *mu_next, gsl_matrix *lambda_next) {
    //for now a very simple implementation
    //no tip forces
    //only forces are gravity and viscosity

    //all the variables and shit
    double n1, n2, n3, n4, n5, n6;
    double M1, M2, M3, M4, M5, M6;
    double L, ds;
    double E, G, A, I, Jz, D, rho, vis;
    int s;

    int N_act = q->size;

    //utility vectors and matrices
    gsl_vector *F = gsl_vector_alloc(6);
    gsl_vector *dxi = gsl_vector_alloc(6);
    gsl_vector *g_row1 = gsl_vector_alloc(12);
    gsl_vector *eta_guess = gsl_vector_alloc(6);
    gsl_vector *eta_temp = gsl_vector_alloc(6);
    gsl_vector *eta_update = gsl_vector_alloc(6);
    gsl_vector *W_bar = gsl_vector_alloc(6);
    gsl_matrix *dtau = gsl_matrix_alloc(6,6);
    gsl_matrix *J = gsl_matrix_alloc(6,6);
    gsl_matrix *g_rowM = gsl_matrix_alloc(4,4);
    gsl_matrix *g_row1M = gsl_matrix_alloc(4,4);
    gsl_permutation *p = gsl_permutation_alloc(6);

    gsl_vector *nu_der = gsl_vector_alloc(3);
    gsl_vector *omega_der = gsl_vector_alloc(3);
    gsl_vector *r = gsl_vector_alloc(3);
    gsl_vector *pa_der = gsl_vector_alloc(3);
    gsl_vector *pa_dder = gsl_vector_alloc(3);
    gsl_vector *ta_der = gsl_vector_alloc(3);
    gsl_matrix *omega_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *omega_der_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *r_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *pa_der_skew = gsl_matrix_alloc(3,3);

    gsl_vector *W_ext = gsl_vector_calloc(6);



    //info
    gsl_vector *xi_ref = gsl_vector_alloc(6);
    gsl_matrix *M = gsl_matrix_calloc(6,6);
    gsl_matrix *K = gsl_matrix_calloc(6,6);
    gsl_matrix *K_inv = gsl_matrix_calloc(6,6);
    gsl_matrix *V = gsl_matrix_calloc(6,6);
    gsl_vector *grav = gsl_vector_alloc(3);

    //temporary vectors and matrices
    gsl_vector *temp1 = gsl_vector_alloc(6);
    gsl_vector *temp2 = gsl_vector_alloc(6);
    gsl_vector *temp3 = gsl_vector_alloc(6);
    gsl_vector *temp4 = gsl_vector_alloc(6);
    gsl_vector *temp5 = gsl_vector_alloc(6);
    gsl_vector *temp6 = gsl_vector_alloc(6);
    gsl_vector *temp7 = gsl_vector_alloc(6);
    gsl_vector *temp8 = gsl_vector_alloc(6);
    gsl_vector *temp9 = gsl_vector_alloc(6);
    gsl_vector *temp10 = gsl_vector_alloc(6);
    gsl_vector *temp11 = gsl_vector_alloc(3);
    gsl_vector *temp12 = gsl_vector_alloc(6);
    gsl_vector *temp13 = gsl_vector_alloc(6);
    gsl_vector *temp14 = gsl_vector_alloc(3);
    gsl_vector *temp15 = gsl_vector_alloc(3);
    gsl_vector *temp16 = gsl_vector_alloc(3);
    gsl_vector *temp17 = gsl_vector_alloc(3);
    gsl_vector *temp18 = gsl_vector_alloc(3);
    gsl_vector *temp19 = gsl_vector_alloc(3);

    gsl_matrix *tempM1 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM2 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM3 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM4 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM5 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM6 = gsl_matrix_alloc(6,6);
    gsl_matrix *tempM7 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM8 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM9 = gsl_matrix_alloc(6,6);
    gsl_matrix *tempM10 = gsl_matrix_alloc(4,4);
    gsl_matrix *tempM11 = gsl_matrix_alloc(6,6);
    gsl_matrix *tempM12 = gsl_matrix_alloc(6,6);
    gsl_matrix *g_backM = gsl_matrix_alloc(4,4);


    int i, j, k;
    for (i=0; i<5; i++) {
        gsl_vector_set(xi_ref, i, 0);
    }
    gsl_vector_set(xi_ref, 5, 1);

    gsl_vector_set(grav,0,0);
    gsl_vector_set(grav,1,0);
    gsl_vector_set(grav,2,0);


    L = 10e-2;
    ds = L/(n-1);
    E = 100e3;
    G = E/3;
    D = 1e-2;
    A = M_PI*D*D/4;
    I = M_PI/64*pow(D,4);
    Jz = 2*I;
    rho = 1000;
    vis = 300;

    M1 = rho*I;
    M2 = rho*I;
    M3 = rho*Jz;
    M4 = rho*A;
    M5 = rho*A;
    M6 = rho*A;
    gsl_matrix_set(M, 0, 0, M1);
    gsl_matrix_set(M, 1, 1, M2);
    gsl_matrix_set(M, 2, 2, M3);
    gsl_matrix_set(M, 3, 3, M4);
    gsl_matrix_set(M, 4, 4, M5);
    gsl_matrix_set(M, 5, 5, M6);

    gsl_matrix_set(K, 0, 0, E*I);
    gsl_matrix_set(K, 1, 1, E*I);
    gsl_matrix_set(K, 2, 2, G*Jz);
    gsl_matrix_set(K, 3, 3, G*A);
    gsl_matrix_set(K, 4, 4, G*A);
    gsl_matrix_set(K, 5, 5, E*A);
    //just since it is easy
    gsl_matrix_set(K_inv, 0, 0, 1/(E*I));
    gsl_matrix_set(K_inv, 1, 1, 1/(E*I));
    gsl_matrix_set(K_inv, 2, 2, 1/(G*Jz));
    gsl_matrix_set(K_inv, 3, 3, 1/(G*A));
    gsl_matrix_set(K_inv, 4, 4, 1/(G*A));
    gsl_matrix_set(K_inv, 5, 5, 1/(E*A));

    gsl_matrix_set(V, 0, 0, vis*3*I);
    gsl_matrix_set(V, 1, 1, vis*3*I);
    gsl_matrix_set(V, 2, 2, vis*Jz);
    gsl_matrix_set(V, 3, 3, vis*A);
    gsl_matrix_set(V, 4, 4, vis*A);
    gsl_matrix_set(V, 5, 5, vis*3*A);

    // the tip and the base are slightly different from the internal points so they have to be done separately
    //also there is a lot of room to reuse and overwrite computations in the internal nodes
    //  for now just make new memory wherever necessary and reuse as much as possible, later can reduce extra memory

    //boundary conditions
    //the external tip load, W_ext = [skew(r)*[0;0;F];[0;0;F]]
    for (i=0; i<N_act; i++) {
        gsl_vector_set(r, 0, D/4*cos(i*2*M_PI/N_act));
        gsl_vector_set(r, 1, D/4*sin(i*2*M_PI/N_act));
        gsl_vector_set(r, 2, 0);
        skew(r,r_skew);
        gsl_vector_set(temp19, 0, 0);
        gsl_vector_set(temp19, 1, 0);
        gsl_vector_set(temp19, 2, gsl_vector_get(q,i));

        gsl_blas_dgemv(CblasNoTrans, 1.0, r_skew, temp19, 0.0, temp18);
        for (j=0; j<3; j++) {
            gsl_vector_set(W_ext, j, gsl_vector_get(W_ext, j) + gsl_vector_get(temp18, j));
            gsl_vector_set(W_ext, j+3, gsl_vector_get(W_ext, j+3) + gsl_vector_get(temp19, j));
        }

    }
    gsl_blas_dgemv(CblasNoTrans, 1.0, K_inv, W_ext, 0.0, temp1);
    adjointSE3(temp1, dtau);
    gsl_matrix_scale(dtau, -ds/2);
    for (i=0; i<6; i++) {
        gsl_matrix_set(dtau, i, i, 1+gsl_matrix_get(dtau,i,i));
    }
    gsl_blas_dgemv(CblasNoTrans, 1.0, K, temp1, 0.0, temp2);
    gsl_vector_add(temp1, xi_ref);
    gsl_blas_dgemv(CblasTrans, 1.0, dtau, temp2, 0.0, temp3);
    for (i=0; i<6; i++) {
        gsl_matrix_set(xi_next, n-1, i, gsl_vector_get(temp1, i));
        gsl_matrix_set(lambda_next, n-1, i, gsl_vector_get(temp3, i));

        gsl_matrix_set(eta_next, 0, i, 0);
        gsl_matrix_set(mu_next, 0, i, 0);
    }
    gsl_matrix_set(g_next, 0, 0, 1);
    gsl_matrix_set(g_next, 0, 1, 0);
    gsl_matrix_set(g_next, 0, 2, 0);
    gsl_matrix_set(g_next, 0, 3, 0);
    gsl_matrix_set(g_next, 0, 4, 1);
    gsl_matrix_set(g_next, 0, 5, 0);
    gsl_matrix_set(g_next, 0, 6, 0);
    gsl_matrix_set(g_next, 0, 7, 0);
    gsl_matrix_set(g_next, 0, 8, 1);
    gsl_matrix_set(g_next, 0, 9, 0);
    gsl_matrix_set(g_next, 0, 10, 0);
    gsl_matrix_set(g_next, 0, 11, 0);

    //the base for xi and lambda
    gsl_vector_const_view eta_curr = gsl_matrix_const_row(eta,0);
    gsl_vector_const_view xi_curr = gsl_matrix_const_row(xi,0);
    gsl_vector_const_view eta_ahead = gsl_matrix_const_row(eta,1);
    for (j=0; j<6; j++) {
        gsl_vector_set(temp1, j, -dt*gsl_vector_get(&eta_curr.vector,j)); // -dt*eta(:,i)
        gsl_vector_set(temp2, j, ds*gsl_vector_get(&xi_curr.vector,j)); // ds*xi(:,i)
        gsl_vector_set(temp3, j, dt*gsl_vector_get(&eta_ahead.vector,j)); //dt*eta(:,i+1)
    }

    expSE3(temp1, tempM1); //exp(-dt*eta(:,i))
    expSE3(temp2, tempM2); //exp(ds*xi(:,i))
    expSE3(temp3, tempM3); //exp(dt*eta(:,i+1))
    //caySE3(temp1, tempM1); //exp(-dt*eta(:,i))
    //caySE3(temp2, tempM2); //exp(ds*xi(:,i))
    //caySE3(temp3, tempM3); //exp(dt*eta(:,i+1))



    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempM1, tempM2, 0.0, tempM4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempM4, tempM3, 0.0, tempM5);


    logSE3(tempM5, temp4); //ds*xi_next
    //icaySE3(tempM5, temp4); //ds*xi_next
    for (j=0; j<6; j++) {
        gsl_matrix_set(xi_next, 0, j, gsl_vector_get(temp4, j)/ds); //xi_next
    }

    //lambda_next(:,i) = dtauInv(ds*xi_next(:,i))'*K*(xi_next(:,i)-xi_ref), dtauinv = I-ad(x)/2
    adjointSE3(temp4, dtau);
    gsl_matrix_scale(dtau,-0.5);
    for (j=0; j<6; j++) {
        gsl_matrix_set(dtau, j, j, 1+gsl_matrix_get(dtau, j, j));
        gsl_vector_set(dxi, j, gsl_vector_get(temp4,j)/ds - gsl_vector_get(xi_ref,j));
    }
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, dtau, K, 0.0, tempM6);
    gsl_blas_dgemv(CblasNoTrans, 1.0, tempM6, dxi, 0.0, temp5);
    for (j=0;j<6;j++) {
        gsl_matrix_set(lambda_next, 0, j, gsl_vector_get(temp5, j)); //lambda_next
    }

    //have to do one step of g outside the loop
    for (j=0; j<6; j++) {
        gsl_vector_set(temp4, j, ds*gsl_matrix_get(xi_next, 0, j));
    }
    expSE3(temp4, tempM7); //exp(ds*xi_next(:,i))
    //caySE3(temp4, tempM7); //exp(ds*xi_next(:,i))
    gsl_vector_const_view g_row = gsl_matrix_const_row(g_next, 0);
    unflattenConfig(&g_row.vector, g_rowM);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g_rowM, tempM7, 0.0, g_row1M);
    flattenConfig(g_row1M,g_row1);
    for (j=0; j<12; j++) {
        gsl_matrix_set(g_next, 1, j, gsl_vector_get(g_row1, j));
    }

    //the internal nodes
    //since the tip node is exactly the same as the internal nodes other than xi and lambda
    //just make the tip a special case for xi and lambda
    for (i=1; i<n; i++) {

        gsl_vector_const_view eta_i = gsl_matrix_const_row(eta, i);
        if (i<(n-1)) {

            // xi_next(:,i) = logSE3(expSE3(-dt*eta(:,i))*expSE3(ds*xi(:,i))*expSE3(dt*eta(:,i+1)))/ds
            gsl_vector_const_view eta_i1 = gsl_matrix_const_row(eta, i+1);
            gsl_vector_const_view xi_i = gsl_matrix_const_row(xi, i);

            for (j=0; j<6; j++) {
                gsl_vector_set(temp1, j, -dt*gsl_vector_get(&eta_i.vector,j)); // -dt*eta(:,i)
                gsl_vector_set(temp2, j, ds*gsl_vector_get(&xi_i.vector,j)); // ds*xi(:,i)
                gsl_vector_set(temp3, j, dt*gsl_vector_get(&eta_i1.vector,j)); //dt*eta(:,i+1)
            }

            expSE3(temp1, tempM1); //exp(-dt*eta(:,i))
            expSE3(temp2, tempM2); //exp(ds*xi(:,i))
            expSE3(temp3, tempM3); //exp(dt*eta(:,i+1))
            //caySE3(temp1, tempM1); //exp(-dt*eta(:,i))
            //caySE3(temp2, tempM2); //exp(ds*xi(:,i))
            //caySE3(temp3, tempM3); //exp(dt*eta(:,i+1))

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempM1, tempM2, 0.0, tempM4);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempM4, tempM3, 0.0, tempM5);

            logSE3(tempM5, temp4); //xi_next*ds
            //icaySE3(tempM5, temp4); //xi_next*ds
            for (j=0; j<6; j++) {
                gsl_matrix_set(xi_next, i, j, gsl_vector_get(temp4, j)/ds); //xi_next
            }

            //lambda_next(:,i) = dtauInv(ds*xi_next(:,i))'*K*(xi_next(:,i)-xi_ref), dtauinv = I-ad(x)/2
            adjointSE3(temp4, dtau);
            gsl_matrix_scale(dtau,-0.5);
            for (j=0; j<6; j++) {
                gsl_matrix_set(dtau, j, j, 1+gsl_matrix_get(dtau, j, j));
                gsl_vector_set(dxi, j, gsl_vector_get(temp4,j)/ds - gsl_vector_get(xi_ref,j));
            }
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, dtau, K, 0.0, tempM6);
            gsl_blas_dgemv(CblasNoTrans, 1.0, tempM6, dxi, 0.0, temp5);
            for (j=0;j<6;j++) {
                gsl_matrix_set(lambda_next, i, j, gsl_vector_get(temp5, j));
            }

            //g_next(:,i) = unflatten(g_next(:,i-1)*exp(ds*xi_next(:,i-1)))
            //technically one step ahead, but no big deal
            expSE3(temp4, tempM7); //exp(ds*xi_next(:,i))
            //caySE3(temp4, tempM7);
            gsl_vector_const_view g_row = gsl_matrix_const_row(g_next, i);
            unflattenConfig(&g_row.vector, g_rowM);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g_rowM, tempM7, 0.0, g_row1M);
            flattenConfig(g_row1M,g_row1);
            for (j=0; j<12; j++) {
                gsl_matrix_set(g_next, i+1, j, gsl_vector_get(g_row1, j));
            }
        } else {
            //need to just pull out xi_next and lambda_next for usage in eta and mu
            //only lambda is used
            for (j=0; j<6; j++) {
                gsl_vector_set(temp5, j, gsl_matrix_get(lambda_next, n-1, j));
            }
        }

        //mu_next(:,i) = Ad(tau(dt*eta(:,i)))'*mu(:,i) + dt/ds*(lambda_next(:,i) - Ad(tau(ds*xi_next(:,i-1))'*lambda_next(:,i-1)))
        gsl_vector_const_view mu_i = gsl_matrix_const_row(mu, i);
        gsl_vector_const_view xi_next_i1 = gsl_matrix_const_row(xi_next, i-1);
        gsl_vector_const_view lambda_next_i1 = gsl_matrix_const_row(lambda_next, i-1);

        for (j=0; j<6; j++) {
            //technically these are operations already done in the previous node step, could optimize away
            gsl_vector_set(temp6, j, gsl_vector_get(&eta_i.vector, j)*dt);
            gsl_vector_set(temp7, j, gsl_vector_get(&xi_next_i1.vector, j)*ds);
        }
        expSE3(temp6, tempM8);
        //caySE3(temp6, tempM8);
        AdjointSE3(tempM8, tempM9);
        expSE3(temp7, tempM10);
        //caySE3(temp7, tempM10);
        AdjointSE3(tempM10, tempM11);

        gsl_blas_dgemv(CblasTrans, 1.0, tempM9, &mu_i.vector, 0.0, temp8);
        gsl_blas_dgemv(CblasTrans, 1.0, tempM11, &lambda_next_i1.vector, 0.0, temp9);

        //external forces (back one time step and spatial step for the simplest implementation)
        for (j=0;j<6;j++) {
            gsl_vector_set(W_bar, j, 0);
        }
        //gravity W_bar = [0;0;0;R'g]*rho*A;
        gsl_vector_const_view g_back = gsl_matrix_const_row(g, i-1);
        unflattenConfig(&g_back.vector, g_backM);
        gsl_matrix_const_view R = gsl_matrix_const_submatrix(g_backM, 0, 0, 3, 3);
        gsl_blas_dgemv(CblasTrans, rho*A, &R.matrix, grav, 0.0, temp11);
        for (j=0; j<3; j++) {
            gsl_vector_set(W_bar, j+3, gsl_vector_get(temp11,j));
        }

        //viscosity W_bar = V*(xi_next(:,i-1)-xi(:,i-1))/dt (simple approximation of xi_dot)
        gsl_vector_const_view xi_i1 = gsl_matrix_const_row(xi, i-1);
        for (j=0; j<6; j++) {
            gsl_vector_set(temp12, j, (gsl_vector_get(&xi_next_i1.vector, j) - gsl_vector_get(&xi_i1.vector, j))/dt);
        }
        gsl_blas_dgemv(CblasNoTrans, -1.0, V, temp12, 0.0, temp13);
        gsl_vector_add(W_bar, temp13);


        //external cable forces W_bar = sum(-q*[skew(r)*R^T*ta_der;R^T*ta_der])
        // ta_der = -skew(pa_der)^2*pa_dder/norm(pa_der)^3
        // pa_der = R*(nu+skew(omega)*r)
        // pa_dder = R*(skew(omega)*(nu+skew(omega)*r)+nu_der+skew(omega_der)*r)
        // nu = nu(i-1,j-1)
        // omega = omega(i-1,j-1)
        // nu_der = (nu(i-1,j)-nu(i-1,j-1))/ds
        // omega_der = (omega(i-1,j)-omega(i-1,j-1))/ds;
        gsl_vector_const_view nu = gsl_vector_const_subvector(&xi_i1.vector, 3,3);
        gsl_vector_const_view omega = gsl_vector_const_subvector(&xi_i1.vector, 0,3);

        gsl_vector_const_view xi_i = gsl_matrix_const_row(xi, i);
        gsl_vector_const_view nu_up = gsl_vector_const_subvector(&xi_i.vector, 3,3);
        gsl_vector_const_view omega_up = gsl_vector_const_subvector(&xi_i.vector, 0,3);
        for (j=0; j<3; j++) {
            gsl_vector_set(nu_der, j, (gsl_vector_get(&nu_up.vector, j) - gsl_vector_get(&nu.vector, j))/ds);
            gsl_vector_set(omega_der, j, (gsl_vector_get(&omega_up.vector, j) - gsl_vector_get(&omega.vector, j))/ds);
        }
        skew(&omega.vector,omega_skew);
        skew(omega_der,omega_der_skew);

        //values that change between actuator
        for (j=0;j<N_act;j++) {
            //load r, for now just assuming the layout is circular and half the radius
            gsl_vector_set(r, 0, D/4*cos(j*2*M_PI/N_act));
            gsl_vector_set(r, 1, D/4*sin(j*2*M_PI/N_act));
            gsl_vector_set(r, 2, 0);
            skew(r,r_skew);

            //pa_der
            gsl_blas_dgemv(CblasNoTrans, 1.0, omega_skew, r, 0.0, temp14);
            gsl_vector_add(temp14,&nu.vector);
            gsl_blas_dgemv(CblasNoTrans, 1.0, &R.matrix, temp14, 0.0, pa_der);
            skew(pa_der, pa_der_skew);

            //pa_dder
            gsl_blas_dgemv(CblasNoTrans, 1.0, omega_skew, temp14, 0.0, temp15);
            gsl_blas_dgemv(CblasNoTrans, 1.0, omega_der_skew, r, 0.0, temp16);
            gsl_vector_add(temp16,temp15);
            gsl_vector_add(temp16,nu_der);
            gsl_blas_dgemv(CblasNoTrans, 1.0, &R.matrix, temp16, 0.0, pa_dder);

            //ta_der
            gsl_blas_dgemv(CblasNoTrans, 1.0, pa_der_skew, pa_dder, 0.0, temp17);
            gsl_blas_dgemv(CblasNoTrans, -1/pow(gsl_blas_dnrm2(pa_der),3), pa_der_skew, temp17, 0.0, ta_der);

            //[skew(r)*R^T*ta_der;R^T*ta_der]
            gsl_blas_dgemv(CblasTrans, -gsl_vector_get(q, j) , &R.matrix, ta_der, 0.0, temp18);
            gsl_blas_dgemv(CblasNoTrans, 1.0, r_skew, temp18, 0.0, temp19);

            //add onto W_bar
            for (k=0; k<3; k++) {
                gsl_vector_set(W_bar,k, gsl_vector_get(W_bar,k) + gsl_vector_get(temp19, k));
                gsl_vector_set(W_bar,k+3, gsl_vector_get(W_bar,k+3) + gsl_vector_get(temp18, k));
            }

        }


        for (j=0; j<6; j++) {
            gsl_vector_set(temp10, j, gsl_vector_get(temp8, j) + dt/ds*(gsl_vector_get(temp5, j) - gsl_vector_get(temp9, j)) + dt*gsl_vector_get(W_bar, j)); //mu_next
            gsl_matrix_set(mu_next, i, j, gsl_vector_get(temp10, j));
        }


        //eta_next(:,i) = solve dtauInv(dt*eta_next(:,i))'*M*eta_next(:,i) - mu_next(:,i)
        //using the simplification that dtauinv = I - ad(x)/2 can do a bit extra since that gives an evaluation to the system and the jacobian pretty easily so can manually just do newton Iteration
        // x_next = x - J(x)\F(x)
        // in matlab testing it usually only 2 steps so could just hard code the 2 steps and assume that its good?

        //setup the guessed value for eta_next
        for (j=0; j<6; j++) {
            gsl_vector_set(eta_guess, j, gsl_vector_get(&eta_i.vector, j));
        }
        for (j=0; j<10; j++) { //solves quickly so very few iterations
            // first compute F(eta) = mu - (I-ad(dt*eta)/2)'*M*eta
            for (k=0; k<6; k++) {
                gsl_vector_set(eta_temp, k, gsl_vector_get(eta_guess,k)*dt);
                gsl_vector_set(F, k, gsl_vector_get(temp10,k)); //initialize with mu_next value
            }
            adjointSE3(eta_temp, dtau);
            gsl_matrix_scale(dtau,-0.5);
            for (k=0; k<6; k++) {
                gsl_matrix_set(dtau, k, k, 1+gsl_matrix_get(dtau, k, k));
            }
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1.0, dtau, M, 0.0, tempM12);
            gsl_blas_dgemv(CblasNoTrans, 1.0, tempM12, eta_guess, 1.0, F);

            //check if F is small
            if (gsl_blas_dnrm2(F) < 1e-10) { //small enough to quit iterating
                break;
            }

            //do the iteration
            //first compute the jacobian, but for simplicity just doing it manually
            /*
            should be negated, whoops
            J =
            [                          M1, (M3*dt*n3)/2 - (M2*dt*n3)/2, (M3*dt*n2)/2 - (M2*dt*n2)/2,                           0, (M6*dt*n6)/2 - (M5*dt*n6)/2, (M6*dt*n5)/2 - (M5*dt*n5)/2]
            [ (M1*dt*n3)/2 - (M3*dt*n3)/2,                          M2, (M1*dt*n1)/2 - (M3*dt*n1)/2, (M4*dt*n6)/2 - (M6*dt*n6)/2,                           0, (M4*dt*n4)/2 - (M6*dt*n4)/2]
            [ (M2*dt*n2)/2 - (M1*dt*n2)/2, (M2*dt*n1)/2 - (M1*dt*n1)/2,                          M3, (M5*dt*n5)/2 - (M4*dt*n5)/2, (M5*dt*n4)/2 - (M4*dt*n4)/2,                           0]
            [                           0,                (M6*dt*n6)/2,               -(M5*dt*n5)/2,                          M4,               -(M5*dt*n3)/2,                (M6*dt*n2)/2]
            [               -(M6*dt*n6)/2,                           0,                (M4*dt*n4)/2,                (M4*dt*n3)/2,                          M5,               -(M6*dt*n1)/2]
            [                (M5*dt*n5)/2,               -(M4*dt*n4)/2,                           0,               -(M4*dt*n2)/2,                (M5*dt*n1)/2,                          M6]

            */
            n1 = gsl_vector_get(eta_guess, 0);
            n2 = gsl_vector_get(eta_guess, 1);
            n3 = gsl_vector_get(eta_guess, 2);
            n4 = gsl_vector_get(eta_guess, 3);
            n5 = gsl_vector_get(eta_guess, 4);
            n6 = gsl_vector_get(eta_guess, 5);

            gsl_matrix_set(J, 0, 0, -M1);
            gsl_matrix_set(J, 0, 1, (M2*dt*n3)/2 - (M3*dt*n3)/2);
            gsl_matrix_set(J, 0, 2, (M2*dt*n2)/2 - (M3*dt*n2)/2);
            gsl_matrix_set(J, 0, 3, 0);
            gsl_matrix_set(J, 0, 4, (M5*dt*n6)/2 - (M6*dt*n6)/2);
            gsl_matrix_set(J, 0, 5, (M5*dt*n5)/2 - (M6*dt*n5)/2);
            gsl_matrix_set(J, 1, 0, (M3*dt*n3)/2 - (M1*dt*n3)/2);
            gsl_matrix_set(J, 1, 1, -M2);
            gsl_matrix_set(J, 1, 2, (M3*dt*n1)/2 - (M1*dt*n1)/2);
            gsl_matrix_set(J, 1, 3, (M6*dt*n6)/2 - (M4*dt*n6)/2);
            gsl_matrix_set(J, 1, 4, 0);
            gsl_matrix_set(J, 1, 5, (M6*dt*n4)/2 - (M4*dt*n4)/2);
            gsl_matrix_set(J, 2, 0, (M1*dt*n2)/2 - (M2*dt*n2)/2);
            gsl_matrix_set(J, 2, 1, (M1*dt*n1)/2 - (M2*dt*n1)/2);
            gsl_matrix_set(J, 2, 2, -M3);
            gsl_matrix_set(J, 2, 3, (M4*dt*n5)/2 - (M5*dt*n5)/2);
            gsl_matrix_set(J, 2, 4, (M4*dt*n4)/2 - (M5*dt*n4)/2);
            gsl_matrix_set(J, 2, 5, 0);
            gsl_matrix_set(J, 3, 0, 0);
            gsl_matrix_set(J, 3, 1, -(M6*dt*n6)/2);
            gsl_matrix_set(J, 3, 2, (M5*dt*n5)/2);
            gsl_matrix_set(J, 3, 3, -M4);
            gsl_matrix_set(J, 3, 4, (M5*dt*n3)/2);
            gsl_matrix_set(J, 3, 5, -(M6*dt*n2)/2);
            gsl_matrix_set(J, 4, 0, (M6*dt*n6)/2);
            gsl_matrix_set(J, 4, 1, 0);
            gsl_matrix_set(J, 4, 2, -(M4*dt*n4)/2);
            gsl_matrix_set(J, 4, 3, -(M4*dt*n3)/2);
            gsl_matrix_set(J, 4, 4, -M5);
            gsl_matrix_set(J, 4, 5, (M6*dt*n1)/2);
            gsl_matrix_set(J, 5, 0, -(M5*dt*n5)/2);
            gsl_matrix_set(J, 5, 1, (M4*dt*n4)/2);
            gsl_matrix_set(J, 5, 2, 0);
            gsl_matrix_set(J, 5, 3, (M4*dt*n2)/2);
            gsl_matrix_set(J, 5, 4, -(M5*dt*n1)/2);
            gsl_matrix_set(J, 5, 5, -M6);

            //don't have to invert J but solve the system J*x_update = F
            gsl_linalg_LU_decomp(J, p, &s);
            gsl_linalg_LU_solve(J, p, F, eta_update);

            //update the guess for eta
            for (k=0; k<6; k++) {
                gsl_vector_set(eta_guess, k, gsl_vector_get(eta_guess,k) - gsl_vector_get(eta_update,k));
            }

        }

        //set eta_next
        for (j=0; j<6; j++) {
            gsl_matrix_set(eta_next, i, j, gsl_vector_get(eta_guess, j));
        }

    }

    //free all the matrices and vectors
    gsl_vector_free(F);
    gsl_vector_free(dxi);
    gsl_vector_free(g_row1);
    gsl_vector_free(eta_guess);
    gsl_vector_free(eta_temp);
    gsl_vector_free(eta_update);
    gsl_vector_free(W_bar);
    gsl_matrix_free(dtau);
    gsl_matrix_free(J);
    gsl_matrix_free(g_rowM);
    gsl_matrix_free(g_row1M);
    gsl_matrix_free(g_backM);
    gsl_permutation_free(p);

    gsl_vector_free(nu_der);
    gsl_vector_free(omega_der);
    gsl_vector_free(r);
    gsl_vector_free(pa_der);
    gsl_vector_free(pa_dder);
    gsl_vector_free(ta_der);
    gsl_matrix_free(omega_skew);
    gsl_matrix_free(omega_der_skew);
    gsl_matrix_free(r_skew);
    gsl_matrix_free(pa_der_skew);

    gsl_vector_free(W_ext);

    gsl_vector_free(xi_ref);
    gsl_vector_free(grav);
    gsl_matrix_free(M);
    gsl_matrix_free(K);
    gsl_matrix_free(K_inv);
    gsl_matrix_free(V);

    gsl_vector_free(temp1);
    gsl_vector_free(temp2);
    gsl_vector_free(temp3);
    gsl_vector_free(temp4);
    gsl_vector_free(temp5);
    gsl_vector_free(temp6);
    gsl_vector_free(temp7);
    gsl_vector_free(temp8);
    gsl_vector_free(temp9);
    gsl_vector_free(temp10);
    gsl_vector_free(temp11);
    gsl_vector_free(temp12);
    gsl_vector_free(temp13);
    gsl_vector_free(temp14);
    gsl_vector_free(temp15);
    gsl_vector_free(temp16);
    gsl_vector_free(temp17);
    gsl_vector_free(temp18);
    gsl_vector_free(temp19);

    gsl_matrix_free(tempM1);
    gsl_matrix_free(tempM2);
    gsl_matrix_free(tempM3);
    gsl_matrix_free(tempM4);
    gsl_matrix_free(tempM5);
    gsl_matrix_free(tempM6);
    gsl_matrix_free(tempM7);
    gsl_matrix_free(tempM8);
    gsl_matrix_free(tempM9);
    gsl_matrix_free(tempM10);
    gsl_matrix_free(tempM11);
    gsl_matrix_free(tempM12);


}
