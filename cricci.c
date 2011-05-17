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
        if (r->s == NULL) {
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
        mesh *m = r->m;
        edge *e;
        int i;

        for (i=0; i < m->ranks[1]; i++) {
                e = &(m->edges[i]);
                calc_edge_length(e);
        }
        calc_inner_angles(m);
        calc_hessian(m);
        if (check_hessian_symmetry(r)) {
                printf("Hessian is symmetric.\n");
        }
        else {
                printf("Hessian is NOT symmetric.\n");
        }
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
        }
}

