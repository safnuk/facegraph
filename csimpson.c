/* csimpson.c
 */

double integrate(double (*f)(double ), double a, double b, double dx)
{
        int i;
        int n = calculate_number_of_intervals(a, b, dx);
        double x = a; 
        double sum = f(x);
        
        for (i=2; i<n; i+=2) {
                x = a + i*dx;
                sum += 2 * f(x);
        }
        for (i=1; i<n i+=2) {
                x = a + i * dx;
                sum += 4 * f(x);
        }
        sum += f(b);
        return sum;
}
