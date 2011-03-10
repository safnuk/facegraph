/* csimpson.c
 */

#include <stdlib.h>
#include <stdio.h>
#include "csimpson.h"

/*
double g(double x, void *instance)
{
        return 4*x*x*x*x + 2/(x*x);
}

int main()
{
        double result;
        double a, b, h;
        a=2; b=5; h=0.1;
        result = integrate(&g, a, b, h, NULL);
        printf("Int(g, %f, %f) = %f  (dx=%f)\n", a, b, result, h);
}
*/

/* Use Simpson's rule to estimate the integral of f over
 * the interval [a, b]. The subdivision used is an even number
 * of subintervals which have width <= dx.
 *
 * The instance variable points to user defined data for use
 * by the callback function.
 */
double integrate(double (*f)(double, void *), 
                double a, double b, double dx, void *instance)
{
        int i;
        double h = dx;
        int n = calculate_number_of_intervals(a, b, &h);
        double x = a; 
        double sum = f(x, instance);
        

        for (i=2; i<n; i+=2) {
                x = a + i*h;
                sum += 2 * f(x, instance);
        }
        for (i=1; i<n; i+=2) {
                x = a + i * h;
                sum += 4 * f(x, instance);
        }
        sum += f(b, instance);
        return sum * h / 3;
}

/* Returns the smallest integer n which is even and satisfies
 *      |b-a| / n  <  dx
 * Also resets dx to be the correct step size ( = |b-a|/n )
 */
int calculate_number_of_intervals(double a, double b, double *dx)
{
        double width = a < b ? (b-a) : (a-b); // width = |b-a|
        int n = (int)(width / (*dx) + 1);
        n =  n % 2 ? n+1 : n; // add 1 to n if it is odd
        *dx = width / n;
        return n;
}
