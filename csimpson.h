// csimpson.h
#ifndef Simpson
#define Simpson

double integrate(double (*f)(double, void *), 
    double a, double b, double dx, void *instance);
int calculate_number_of_intervals(double a, double b, double *dx);

#endif
