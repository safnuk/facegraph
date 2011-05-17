/* cricci.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "csimpson.h"
#include "cfile_io.h"
#include "cmesh.h"
#include "cricci.h"

void run_ricci_flow(mesh *m)
{
        ricci_solver r;
        ricci_config rc;
        initialize_ricci_solver(&r, m, &rc);
        calc_flat_metric(&r);
        deallocate_ricci_solver(&r);
}

void initialize_ricci_solver(ricci_solver *r, mesh *m, ricci_config *rc)
{
        // set default configuration settings
        rc->verbose = 2; // set to 2 for updates every iteration, 0 for no output
        rc->ds = 0.1; // step length for numerical integration
        rc->relative_error = 1.0e-6;
        rc->absolute_error = 1.0e-12;
        rc->wolfe_c1= 1.0e-4;
        rc->wolfe_c2 = 0.9;
        rc->max_iterations = 20;
        r->m = m;
        r->rc = rc;
        r->iteration = 0;
        r->status = RUNNING;
        r->s = (double *)malloc(m->ranks[0] * sizeof(double));
        r->K = (double *)malloc(m->ranks[0] * sizeof(double));
        r->step = (double *)malloc(m->ranks[0] * sizeof(double));
        if (!r->s || !r->K || !r->step) {
                printf("Memory allocation error in initialize_ricci_solver.\n");
                exit(1);
        }
        calc_initial_variables(r);
}

void deallocate_ricci_solver(ricci_solver *r)
{
        free(r->s);
        free(r->K);
        free(r->step);
}

void calc_flat_metric(ricci_solver *r)
{
        while (convergence_test(r) == RUNNING) {
                update_hessian(r);
                calc_next_step(r);
                calc_line_search(r);
        }
}


/* Tests if ricci flow has converged, and updates status
 * of ricci_solver.
 */
ricci_state convergence_test(ricci_solver *r)
{
        r->iteration++;
        if (r->iteration > r->rc->max_iterations) {
                r->status = TOO_MANY_ITERATIONS; 
        }
        r->K_norm = sup_norm(r->K, r->m->ranks[0]);
        if (r->K_norm < r->rc->relative_error * r->init_K_norm + r->rc->absolute_error) {
                r->status = CONVERGENT;
        }
        print_ricci_status(r);
        return r->status;
}

void update_hessian(ricci_solver *r)
{
}

void calc_next_step(ricci_solver *r)
{
}

void calc_line_search(ricci_solver *r)
{
}

/* Calculates the product H*x and returns the result in y,
 * where H is the hessian [dK_i / du_j] (symmetric, positive definite).
 *
 * n is the size (number of vertices) and instance is a pointer
 * to a ricci_solver struct.
 *
 * Assumes that calc_hessian(m) has been called for the mesh
 * being used.
 */
int calc_hessian_product(double *x, double *y, 
    void *instance)
{
        int i, j, k, l;
        mesh *m;
        ricci_solver *rs;
        vertex *v, *v1;
        triangle *t;


        rs = (ricci_solver *)instance;
        m = rs->m;
        for (i=0; i < m->ranks[0]; i++) {
                y[i] = 0;
                v = &(m->vertices[i]);
                for (j=0; j < v->degree - v->boundary; j++) {
                        t = (triangle *)(v->incident_triangles[j]);
                        for (k=0; k<3; k++) {
                                v1 = (vertex *)(t->vertices[k]);
                                l = v1->index;
                                y[i] -= x[l] * *(v->dtheta_du[j][k]);
                        }
                }
        }
}

void calc_initial_variables(ricci_solver *r)
{
        int i;
        for (i=0; i < r->m->ranks[0]; i++) {
                r->s[i] = r->m->vertices[i].s;
                r->step[i] = 0;
        }
        calc_edge_lengths(r->m);
        calc_curvatures(r->m, r->K);
        r->init_K_norm = sup_norm(r->K, r->m->ranks[0]);
        r->K_norm = r->init_K_norm;
        r->step_norm = 0;
        r->step_scale = 1;
}

int check_hessian_symmetry(ricci_solver *r)
{
        mesh *m = r->m;
        double *x, *y;
        int i,j;
        double h_ij, h_ji;
        x = (double *) malloc(m->ranks[0] * sizeof(double));
        y = (double *) malloc(m->ranks[0] * m->ranks[0] * sizeof(double));
        if (!x || !y) {
                printf("Memory allocation error in check_hessian_symmetry.\n");
                exit(1);
        }
        for (i=0; i < m->ranks[0]; i++) {
                x[i] = 0;
        }
        j=0;
        for (i=0; i < m->ranks[0]; i++) {
                x[i] = 1;
                calc_hessian_product(x, &y[j * m->ranks[0]], (void *)r);
                j++;
                x[i]= 0;
        }
        for (i=0; i < m->ranks[0]; i++) {
                for (j=i+1; j < m->ranks[0]; j++) {
                        h_ij = y[i * m->ranks[0] + j];
                        h_ji = y[j * m->ranks[0] + i];
                        if (abs(h_ij -  h_ji) > .00000001) {
                                free(x);
                                free(y);
                                return 0;
                        }
                }
        }
        free(x);
        free(y);
        return 1;
}

void print_ricci_status(ricci_solver *r)
{
        if ((r->status != RUNNING) && (r->rc->verbose > 0)) {
                printf("Ricci flow terminated after %i iterations with status code %i, ", r->iteration-1, 
                                r->status);
                printf("||K|| = %f, Step size = %f\n", r->K_norm, r->step_norm);
        }
        if ((r->status == RUNNING) && (r->rc->verbose == 2)) {
                printf("Iteration %i: ||K|| = %f, Step size = %f, Step scale = %f\n",
                     r->iteration-1, r->K_norm, r->step_norm, r->step_scale );
        }
}
/* Calculates the sup-norm of the n-dimensional vector v. 
 */
double sup_norm(double *v, int n)
{
        int i;
        double max = 0;
        for (i=0; i<n; i++) {
                if (abs(v[i]) > max) {
                        max = abs(v[i]);
                }
        }
        return max;
}
