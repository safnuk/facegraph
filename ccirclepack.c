// ccirclepack.c

#include <stdlib.h>
#include <stdio.h>
#include <lbfgs.h>
#include <math.h>
#include <list>
#include "cfile_io.h"
#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"
#include "geodesic.h"
#include "cmesh.h"
#include "ccirclepack.h"


/* Note: lbfgsfloatval_t is a double type.
 *
 */

double calc_circlepack_metric(mesh *m)
{
        optimizer opt;
        initialize_optimizer(&opt,  m);
        calc_optimum_radii(&opt);
        record_output_in_mesh(&opt);
        deallocate_optimizer(&opt);
}

void initialize_optimizer(optimizer *opt, mesh *m)
{
        opt->verbose = 1;  // set to 2 for verbose status outputs, 0 for no output
        opt->m = m;
        opt->ambient_lengths = lbfgs_malloc(m->ranks[1]);
        opt->r = lbfgs_malloc(m->ranks[0]);
        if ((opt->ambient_lengths == NULL) || (opt->r == NULL)) { 
                printf("Memory allocation error");
                exit(1);
        }
        calc_ambient_edge_lengths(opt);
        calc_initial_radii(opt);
        // use default values for optimization routine
        lbfgs_parameter_init(&(opt->param));
}

void deallocate_optimizer(optimizer *opt)
{
        lbfgs_free(opt->ambient_lengths);
        lbfgs_free(opt->r);
}

/* Calculates the lengths of all edges in the mesh,
 * as inherited from the coordinates of the vertices
 * in 3-space. i.e., edge lengths inherited from the
 * ambient metric.
 */
void calc_ambient_edge_lengths(optimizer *opt)
{
        int i;
        point *v1, *v2;
        mesh *m = opt->m;
        lbfgsfloatval_t *edge_lengths = opt->ambient_lengths;

        for (i=0; i < m->ranks[1]; i++)
        {
                v1 = get_coordinate(m, 
                                (vertex *)((m->edges[i]).vertices[0]));
                v2 = get_coordinate(m, 
                                (vertex *)((m->edges[i]).vertices[1]));
                edge_lengths[i] = calc_distance(v1, v2);
        }
}

/* Sets initial radii to be 70% of the average edge length.
 * I.e., for each vertex, 70% of the average length of
 * each incident edge is used for the inital radius.
 */
void calc_initial_radii(optimizer *opt)
{
        double radius;
        int i, j, ei;
        vertex *v;
        edge *e;
        for (i=0; i < opt->m->ranks[0]; i++) {
                radius = 0;
                v = &(opt->m->vertices[i]);
                for (j=0; j < v->degree; j++) {
                        ei = ((edge *)(v->incident_edges[j]))->index;
                        radius += opt->ambient_lengths[ei];
                }
                if (v->degree != 0) {
                        radius = (radius * 0.70) / v->degree;
                }
                opt->r[i] = radius;
        }
}

/* TODO: Implement verbose == 1 routine (display end status)
 */
void calc_optimum_radii(optimizer *opt)
{
        int status_code;
        lbfgsfloatval_t f_of_r;
        state *s = &(opt->current_state);
        if (opt->verbose >= 2) {
                status_code = lbfgs(opt->m->ranks[0], opt->r, &f_of_r,
                                evaluate, progress, opt, &(opt->param));
        } else if (opt->verbose == 1) {
                status_code = lbfgs(opt->m->ranks[0], opt->r, &f_of_r,
                                evaluate, quiet_progress, opt, &(opt->param));
                printf("Optimizer finished %i iterations with status code %i.\n", s->iterations, status_code);
                printf("   Error = %f  ||r|| = %f  ||Df(r)|| = %f Step = %f\n", s->error, s->r_norm,
                                s->g_norm, s->step);
        } else {
                status_code = lbfgs(opt->m->ranks[0], opt->r, &f_of_r,
                                evaluate, NULL, opt, &(opt->param));
        }
}

void record_output_in_mesh(optimizer *opt)
{
        copy_radii(opt);
        calc_circle_angles(opt);
}

/* Copy optimized radii from optimizer struct into
 * mesh struct, to save for future use.
 */
void copy_radii(optimizer *opt)
{
        int i;
        vertex *v;
        for(i=0; i < opt->m->ranks[0]; i++) {
                v = &(opt->m->vertices[i]);
                v->s = (exp(opt->r[i]) - 1) / (exp(opt->r[i]) + 1);
                v->s0 = v->s;
        }
}

/* Calculate the angles of intersection between circles
 * from incident vertices. i.e. there is an intersection
 * at every edge, and the angle is chosen to make
 * the circle packing metric as close as possible to 
 * the ambient metric.
 * 
 * This routine also calculates the induced edge length
 * from the circle packing metric, and stores it in the
 * edge struct (along with the intersection angle).
 *
 * Technically, what is stores is the cosine of the angle
 * and the hyperbolic cosine of the length, as these are the
 * natural parameters used for Ricci flow.
 */
void calc_circle_angles(optimizer *opt)
{
        int i, vi0, vi1;
        double r0, r1, l_bar, tmp;
        vertex *v0, *v1;
        edge *e;

        for (i=0; i < opt->m->ranks[1]; i++) {
                e = &(opt->m->edges[i]);
                v0 = (vertex *)(e->vertices[0]);
                v1 = (vertex *)(e->vertices[1]);
                vi0 = v0->index;
                vi1 = v1->index;
                r0 = opt->r[vi0];
                r1 = opt->r[vi1];
                l_bar = opt->ambient_lengths[i];
                if (l_bar >= r0 + r1) { // circles are tangent
                        e->cosh_length = cosh(r0 + r1);
                        e->cos_angle = -1;
                } else if (cosh(l_bar) <= cosh(r0) * cosh(r1)) { // maximum allowed overlap
                        e->cosh_length = cosh(r0) * cosh(r1);
                        e->cos_angle = 0; // angle = pi / 2
                } else { // arbitrary obtuse angle (no error)
                        e->cosh_length = cosh(l_bar);
                        tmp = sinh(r0) * sinh(r1);
                        if (tmp != 0.0) {
                                e->cos_angle = 
                                        ((cosh(r0) * cosh(r1) - cosh(l_bar)) / tmp);
                        } else {
                                e->cos_angle = 0;
                        }
                }
        }
}

/* Callback function used by lbfg routine to
 * evaluate the function being minimized
 * and its gradient at the point x.
 * isinstance is used to store parameters used by the function,
 * g is the array of gradient outputs, and n is the number of variables.
 *
 * In this case, the function being minimized is the sum of squares
 * of residuals between the ambient edge lengths and the lengths
 * induced by the circle packing metric (with vertex radii being the 
 * variables). 
 *
 * instance should be a pointer to the optimizer struct.
 */
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *r,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
        optimizer *opt = (optimizer *)instance;
        mesh *m = opt->m;
        lbfgsfloatval_t total_error = 0;
        double r0, r1, l_bar, l_actual;
        double cosh_r0, cosh_r1;
        double error, deriv;
        int i, vi0, vi1;

        for (i=0; i < n; i++) {
                g[i] = 0;
        }
        for (i=0; i < m->ranks[1]; i++) {
                vi0 = ((vertex *)((m->edges[i]).vertices[0]))->index;
                vi1 =((vertex *)(((m->edges)[i]).vertices[1]))->index;
                r0 = r[vi0];
                r1 = r[vi1];
                cosh_r0 = cosh(r0);
                cosh_r1 = cosh(r1);
                l_bar = opt->ambient_lengths[i];
                if (l_bar > r0 + r1) {
                        error = r0 + r1 - l_bar;
                        total_error += error * error;
                        deriv = 2 * (r0  + r1 - l_bar);
                        g[vi0] += deriv;
                        g[vi1] += deriv;
                } else if (cosh(l_bar) < cosh_r0 * cosh_r1) {
                        l_actual = acosh(cosh_r0 * cosh_r1);
                        error = l_actual - l_bar;
                        total_error += error * error;
                        deriv = 2 * error / sinh(l_actual);
                        g[vi0] += deriv * sinh(r0) * cosh_r1;
                        g[vi1] += deriv * sinh(r1) * cosh_r0;
                }
        }
        return total_error;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *r,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t f_of_r,
    const lbfgsfloatval_t r_norm,
    const lbfgsfloatval_t g_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
        printf("Iteration %d:\n", k);
        printf("  f(r) = %f, ||r|| = %f, ||Df(r)|| = %f, step = %f\n",
                        f_of_r, r_norm, g_norm, step);
        return 0;
}

static int quiet_progress(
    void *instance,
    const lbfgsfloatval_t *r,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t f_of_r,
    const lbfgsfloatval_t r_norm,
    const lbfgsfloatval_t g_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
        state *update = &(((optimizer*)instance)->current_state);
        update->error = f_of_r;
        update->r_norm = r_norm;
        update->g_norm = g_norm;
        update->step = step;
        update->iterations = k;

        return 0;
}
