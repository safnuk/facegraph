/* cricci.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lbfgs.h>
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
        r->u = lbfgs_malloc(m->ranks[0]);
        if (r->u == NULL) {
                printf("Memory allocation error.\n");
                exit(1);
        }
        calc_initial_variables(r);
        // use default values for optimization routine
        lbfgs_parameter_init(&(r->param));
        // r->param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_MORETHUENTE;

}

void deallocate_ricci_solver(ricci_solver *r)
{
        lbfgs_free(r->u);
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

void calc_initial_variables(ricci_solver *r)
{
        int i;
        double s;
        for (i=0; i < r->m->ranks[0]; i++) {
                s = r->m->vertices[i].s;
                r->u[i] = log(s);
        }
}

static lbfgsfloatval_t ricci_evaluate(
    void *instance,
    const lbfgsfloatval_t *u,
    lbfgsfloatval_t *K,
    const int n,
    const lbfgsfloatval_t step
    )
{
        ricci_solver *r = (ricci_solver *)instance;
        mesh *m = r->m;
        double K_avg;
        lbfgsfloatval_t f;

        printf("Range for u: %f <= u <= %f; ", min(u,n), max(u,n));
        f = update_f_and_s(m, u, r->ds);
        printf("f = %.4e; ", f);
        K_avg = calc_curvatures(m, K);
        printf("%f <= K <= %f; Avg K = %f\n", min(K,n), max(K,n), K_avg);
        return f;
}

static int ricci_progress(
    void *instance,
    const lbfgsfloatval_t *u,
    const lbfgsfloatval_t *K,
    const lbfgsfloatval_t f_of_u,
    const lbfgsfloatval_t u_norm,
    const lbfgsfloatval_t K_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
        ricci_state *update = &(((ricci_solver *)instance)->current_state);
        update->f_of_u = f_of_u;
        update->u_norm = u_norm;
        update->K_norm = K_norm;
        update->step = step;
        update->iterations = k;

        if (((ricci_solver *)instance)->verbose == 2) {
                printf("Iteration %i: ", k);
                printf("f(u) = %f, ||u|| = %f, ||K(u)|| = %f, step = %f\n", 
                                f_of_u, u_norm, K_norm, step);
        }
        return 0;
}
