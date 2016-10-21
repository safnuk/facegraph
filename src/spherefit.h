// spherefit.h
#ifndef SphereFit
#define SphereFit

struct sphere_data {
  size_t n;
  point* points;
};

void fit_points_to_sphere();
void initialize sphere_data(sphere_data* sd, points* p, int n);
void iterate_sphere_solver(s);
void print_final_sphere_status(s);
void print_sphere_state(s);

#endif
