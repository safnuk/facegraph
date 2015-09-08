// spherefit.c

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

void fit_points_to_sphere(point* points, int n) 
{
  const gsl_multifit_fdfsolver_type* T;
  gsl_multifit_fdfsolver* s;
  sphere_data sd;
  const size_t p = 4; // number of parameters to fit

  gsl_matrix* covariance = gsl_matrix_alloc(p,p);
  initialize_sphere_data(&sd, points, n);

  gsl_multifit_function_fdf f;
  double x_init[4] = {1.0, 1.0, 1.0, 1.0};
  gsl_vector_view x = gsl_vector_view_array(x_init, p);

  f.f = &sphere_f;
  f.df = &sphere_df;
  f.fdf = &sphere_fdf;
  f.n = n;
  f.p = p;
  f.params = &sd;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);
  gsl_multifit_fdfsolver_set(s, &f, &x.vector);

  iterate_sphere_solver(s);

  print_final_status(s);

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covariance);
}

void initialize sphere_data(sphere_data* sd, point* p, int n)
{
  sd->n = n;
  sd->points = p;
}

void iterate_sphere_solver(s)
{
  int iterations = 0;
  int status;
  print_sphere_state(iterations, s);
  do {
    iterations++;
    status = gsl_multifit_fdfsolver_iterate(s);
    printf("status = %s\n", gsl_strerror(status));
    print_sphere_state(iterations, s);
    if (status)
      break;
    status = gsl_mutlifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE && iterations < 500);
}

void print_final_sphere_status(s)
{
}

void print_sphere_state(s)
{
}
