/* cricci.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "csimpson.h"
#include "cfile_io.h"
#include "coord_double.h"
#include "comprow_double.h"
#include "mvblasd.h"
#include "diagpre_double.h"              // Preconditioners
#include "icpre_double.h"
#include "ilupre_double.h"
#include "cg.h"
#include "cmesh.h"
#include "cconj_grad.h"
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
        rc->verbose = 3; // set to 2 for updates every iteration, 0 for no output
        rc->integration_precision = 30; // number of subdivisions for simpson's integration
        rc->relative_error = 1.0e-11;
        rc->absolute_error = 1.0e-11;
        rc->wolfe_c1= 1.0e-4;
        rc->wolfe_c2 = 0.9;
        rc->strong_wolfe = 0;  // 0 for regular Wolfe conditions, 1 for strong
        rc->max_iterations = 30;
        rc->max_line_steps = 25;
        rc->cg_max_iterations = 200; // 0 for no max number of iterations
        rc->cg_precon = 0;  // 1 for Jacobi preconditioning, 0 for no precon.
        rc->cg_tolerance = 1.0e-6;
        r->m = m;
        r->rc = rc;
        r->iteration = 0;
        r->status = RUNNING;
        r->s.newsize(m->ranks[0]);
        r->s_next.newsize(m->ranks[0]);
        r->K.newsize(m->ranks[0]);
        r->K_next.newsize(m->ranks[0]);
        r->step.newsize(m->ranks[0]);
        calc_initial_variables(r);
}

void deallocate_ricci_solver(ricci_solver *r)
{
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
        r->K_norm = norm(r->K);
        r->step_norm = norm(r->step);
        if (r->K_norm < r->rc->relative_error * r->init_K_norm + r->rc->absolute_error) {
                r->status = CONVERGENT;
        }
        print_ricci_status(r);
        return r->status;
}

void update_hessian(ricci_solver *r)
{
        calc_hessian(r->m);
}

/* Use conjugate gradient method to solve the equation
 *      H step = K
 * This gives the direction to change the radii parameters which
 * most quickly approaches the minimium of f (where curvature is 0).
 */
void calc_next_step(ricci_solver *r)
{
        int result;
        int max_iterations = r->rc->cg_max_iterations;
        double tolerance = r->rc->cg_tolerance;
        r->step = 0;
        DiagPreconditioner_double D(r->m->hessian);

        result = CG(r->m->hessian, r->step, r->K, D, max_iterations, tolerance);

        if (r->rc->verbose > 2) {
                printf("CG flag = %i, iterations = %i, tolerance = %e\n", result,
                                max_iterations, tolerance);
        }
        /*
        if (r->rc->cg_precon) {
                pccg_solve(&calc_hessian_product, r->step, r->K, 
                        r->rc->cg_tolerance,
                        r->m->ranks[0], r->rc->cg_max_iterations, 
                        r->rc->verbose - 2,
                        (void *)r);
        } else {
                cg_solve(&calc_hessian_product, r->step, r->K, 
                        r->rc->cg_tolerance,
                        r->m->ranks[0], r->rc->cg_max_iterations, 
                        r->rc->verbose - 2,
                        (void *)r);
        }
        */
}

/* Perform a line search in the direction found by calc_next_step.
 * Primarily, this is to ensure
 * (a) That the step keeps the s-variables in the proper bounds
 *     (between 0 and 1).
 *
 * (b) Wolfe conditions are satisfied (to ensure convergence of
 *     the algorithm).
 *
 * The step direction is scaled by factors of (1/2)^i until both
 * (a) and (b) are met. If this does not happen within the max
 * number of iterations allowed, the status of ricci_solver is
 * updated accordingly.
 */
void calc_line_search(ricci_solver *r)
{
        int i=0;
        int wolfe_conditions_verified=0;
        r->step_scale = 2;
        r->K_supnorm = sup_norm(r->K);
        while (!wolfe_conditions_verified && (i < r->rc->max_line_steps)) {
                r->step_scale *= 0.5; 
                calc_next_s(r);
                if (vector_in_bounds(r->s_next, 1.0e-14, 1-1.0e-14)) {
                        update_s_and_edge_lengths(r->m, r->s_next);
                        calc_curvatures(r->m, r->K_next);
                        wolfe_conditions_verified = test_wolfe_conditions(r);
                }
                i++;
        }
        if (!wolfe_conditions_verified) {
                r->status = LINE_SEARCH_FAILED;
                return;
        }
        // Make the next-step quantities the current values
        r->K = r->K_next;
        r->s = r->s_next;
}

/* Calculates the new value for the radius parameters,
 * by the formula
 *      s_next = s * exp((1/2)^j step)
 *
 * This is because the gradient descent variables are
 * u, with s = exp(u), hence
 *      u_next = u + lambda * step
 * translates into the above formula.
 * The formula comes from solving the equation
 *  -H * step = K for Newton's method.
 */
void calc_next_s(ricci_solver *r)
{
        int i;
        for (i=0; i < r->m->ranks[0]; i++) {
                r->s_next[i] = r->s[i] * exp(r->step_scale * r->step[i]);
        }
}

/* Checks if Wolfe conditions are satisfied for the line search
 * being performed. Note that signs are opposite of the usual
 * because the step vector is multiplied by -1.
 *
 * Returns 1 if Wolfe conditions are satisfied, 0 otherwise.
 */
int test_wolfe_conditions(ricci_solver *r)
{
        double K_next_supnorm = sup_norm(r->K_next);
        if (K_next_supnorm >= r->K_supnorm) {
                return 0;
        }
        /*
        if (r->rc->strong_wolfe) {
                if (fabs(dot_product(r->step, r->K_next, n)) >
                                r->rc->wolfe_c2 * fabs(dot_product(r->step, r->K, n))) {
                        return 0;
                }
        } else {
                if (dot_product(r->step, r->K_next, n) >
                        r->rc->wolfe_c2 * dot_product(r->step, r->K, n)) {
                        return 0;
                }
        }
        */
        return 1;
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
        r->init_K_norm = norm(r->K);
        r->K_norm = r->init_K_norm;
        r->step_norm = 0;
        r->step_scale = 1;
}


void print_ricci_status(ricci_solver *r)
{
        if ((r->status != RUNNING) && (r->rc->verbose > 0)) {
                printf("Ricci flow terminated after %i iterations with status code %i, ", r->iteration-1, 
                                r->status);
                printf("||K|| = %e, K_max = %e, Step size = %f\n", r->K_norm, 
                                sup_norm(r->K), r->step_norm);
        }
        if ((r->status == RUNNING) && (r->rc->verbose >= 2)) {
                printf("Iteration %i: ||K|| = %e, K_max = %e, Step size = %f, Step scale = %f\n",
                     r->iteration-1, r->K_norm, 
                     sup_norm(r->K), r->step_norm, r->step_scale);
        }
}
/* Calculates the sup-norm of the n-dimensional vector v. 
 */
double sup_norm(const MV_Vector_double &x)
{
  double temp = 0;
  for (int i=0; i<x.size(); i++) {
    if (temp < fabs(x(i)))
      temp = fabs(x(i));
  }
  return temp;
}


/* Returns true if every entry of the vector satisfies
 *              min < v[i] < max.
 *  Returns false otherwise.
 */
int vector_in_bounds(const MV_Vector_double &v, double min, double max)
{
        int i;
        for (i=0; i<v.size(); i++) {
                if (v[i] <= min || v[i] >= max) {
                        return 0;
                }
        }
        return 1;
}


