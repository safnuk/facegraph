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
        initialize_ricci_solver(&r, m);
        calc_flat_metric(&r);
        deallocate_ricci_solver(&r);
}

void initialize_ricci_solver(ricci_solver *r, mesh *m)
{
        r->verbose = 2; // set to 2 for updates every iteration, 0 for no output
        r->ds = 0.1; // step length for numerical integration
        r->m = m;
        r->s = (double *)malloc(m->ranks[0] * sizeof(double));
        if (r->u == NULL) {
                printf("Memory allocation error.\n");
                exit(1);
        }
        calc_initial_variables(r);
}

void deallocate_ricci_solver(ricci_solver *r)
{
        free(r->s);
}

void calc_flat_metric(ricci_solver *r)
{
        int status_code;
        lbfgsfloatval_t f_of_u;
        ricci_state *s = &(r->current_state);
        status_code = lbfgs(r->m->ranks[0], r->u, &f_of_u,
                        ricci_evaluate, ricci_progress, r, &(r->param));
        if (r->verbose >= 1) {
                printf("Ricci solver finished %i iterations with status code %i.\n",
                                s->iterations, status_code);
                printf("    f(u) = %f, ||u|| = %f, ||K(u)|| = %f\n", 
                                f_of_u, s->u_norm, s->K_norm);
        }
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
int calc_hessian_product(double *x, double *y, int n,
    void *instance)
{
        int i, j, k, l;
        mesh *m;
        ricci_solver *rs;
        vertex *v, *v1;
        triangle *t;


        rs = (ricci_solver *)instance;
        m = rs->m;
        for (int i=0; i<n; i++) {
                y[i] = 0;
                v = &(m->vertices[i]);
                for (j=0; j < v->degree - v->boundary; j++) {
                        t = (triangle *)(v->incident_triangles[j])
                        for (k=0; k<3; k++) {
                                v1 = (vertex *)(t->vertices[k]);
                                l = v1->index;
                                y[i] -= x[l] * (v->dtheta_du[j][k])
                        }
                }
        }
}

void calc_initial_variables(ricci_solver *r)
{
        int i;
        for (i=0; i < r->m->ranks[0]; i++) {
                r->s[i] = r->m->vertices[i].s;
        }
}

