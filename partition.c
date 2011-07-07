// partition.c

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <gsl/gsl_multimin>

#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"

#include "cfile_io.h"
#include "geodesic.h"
#include "cmesh.h"
#include "graph.h"
#include "partition.h"

const double kStepSize = 0.1;
const double kTol = 1e-4;
const double kErrTol = 1e-6;
const int kMaxIterations = 100;

double calc_optimal_partition_offset(vertex* v)
{
  gsl_multimin_function_fdf my_func;
  my_func.n = 2;
  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  vertex_params p[3];
  calc_vertex_params(v, p);
  my_func.params = (void*)p;

  gsl_multimin_fdfminimizer_type const* T = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, 2);
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, 0.0);
  gsl_vector_set(x, 1, 0.0);

  gsl_multimin_fdminimizer_set(s, &my_func, &x, kStepSize, kTol);

  int iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status) break;
    status = gsl_multimin_test_gradient(s->gradient, kErrTol);
    if (status == GSL_SUCCESS) {
      printf("Minimum found at:\n");
    }
    printf("%5d %.5f %.5f %10.5f\n", iter, 
        gsl_vector_get(s->x, 0),
        gsl_vector_get(s->x, 1),
        s->f);
  } while (status == GSL_CONTINUE && iter < kMaxIterations);

  double offset = calc_offset(p[0].length, gsl_vector_get(s->x, 0),
      gsl_vector_get(s->x, 1));

  gsk_multimin_fdminimizer_free(s);
  gsl_vector_free(x);
  return offset;
}

double calc_offset(double l, double alpha, double d)
{
  double a = normalize_angle(alpha);
  double sign = 1.0;
  if (a > M_PI) {
    sign = -1;
  }
  double l_new = asinh(sinh(l) * cosh(d) - cosh(l) * sinh(d) * cos(a));
  double offset = acosh((sinh(l_new) * sinh(l) + cosh(d)) / (cosh(l_new) * cosh(l)));
  return sign * offset;
}

double my_f(const gsl_vector* v, void* params)
{
  double d[3];
  double alphas[3];
  vertex_params* p = (vertex_params*) params;
  double alpha = gsl_vector_get(v, 0);
  double delta = gsl_vector_get(v, 1);
  calc_new_lengths(alpha, delta, p, d, alphas);
  double result=0;
  for (int i=0; i<3; ++i) {
    result += (d[i] - d[(i+1)%3]) * (d[i] - d[(i+1)%3]);
  }
  return result;
}

void my_df(const gsl_vector* v, void* params, gsl_vector* df)
{
  double d[3];
  double alphas[3];
  double dl_dalpha[3];
  double dl_ddelta[3];
  vertex_params* p = (vertex_params*) params;
  double alpha = gsl_vector_get(v, 0);
  double delta = gsl_vector_get(v, 1);
  calc_new_lengths(alpha, delta, p, d, alphas);
  calc_dl_dalpha(alphas, delta, p, d, dl_dalpha);
  calc_dl_ddelta(alphas, delta, p, d, dl_ddelta);
  double df_da = 0;
  double df_dd = 0;
  for (int i=0; i<3; ++i) {
    double c = 4 * d[i] - 2*d[(i+1)%3] - 2*d[(i+2)%3];
    df_da += c * dl_dalpha[i];
    df_dd += c * dl_ddelta[i];
  }
  gsl_vector_set(df, 0, df_da);
  gsl_vector_set(df, 1, df_dd);

}

void my_fdf(const gsl_vector* v, void* params, double* f gsl_vector* df)
{
  (*f) = my_f(v, params);
  my_df(v, params, df);
}

void calc_new_lengths(double alpha, double delta, 
    vertex_params const* p, double* d, double* alphas)
{
  alphas[0] = alpha; 
  alphas[1] = alpha + p[1].angle + p[2].angle;
  alphas[2] = alpha + p[2].angle;
  for (int i=0; i<3; ++i) {
    d[i] = calc_new_length(alphas[i], delta, p[i].length);
  }
}

double calc_new_length(double alpha, double delta, double length)
{
  double result = sinh(length) * cosh(delta);
  result -= cosh(length) * sinh(delta) * cos(alpha);
  return asinh(result);
}

void calc_dl_dalpha(double const* alphas, double delta, 
    vertex_params const* p, double const *d, double* dl_dalpha)
{
  for (int i=0; i<3; ++i) {
    dl_dalpha[i] = cosh(p[i].length) * sinh(delta) * sin(alphas[i]) / cosh(d[i]);
  }
}

void calc_dl_ddelta(double const* alphas, double delta, 
    vertex_params const* p, double const* d, double* dl_ddelta)
{
  for (int i=0; i<3; ++i) {
    dl_ddelta[i] = (sinh(p[i].length) * sinh(delta) - cosh(p[i].length) * cosh(delta)
        / cos(alphas[i]);
  }
}
