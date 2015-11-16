/* File:    ex5.9_omp_trap3_schedule.c
 * Purpose: Estimate definite integral (or area under curve) using the
 *          trapezoidal rule.  This version uses a parallel for directive.
 *          It also uses a schedule(runtime) clause to modify the parallel
 *          for directive.  This is a modified version of omp_trap3.c.
 *
 * Compile: gcc -g -Wall -fopenmp -o ex5.9_omp_trap3_schedule 
 *                ex5.9_omp_trap3_schedule.c
 * Usage:   ./ex5.9_omp_trap3_schedule <number of threads>
 *
 * Input:   a, b, n
 * Output:  Which thread executes which iteration of the Trap for loop
 *          and an estimate of the integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Notes:   
 *   1.  The function f(x) is hardwired.
 *   2.  In this version, it's not necessary for n to be
 *       evenly divisible by thread_count.
 *   3.  Assumes n >= 2
 *
 * IPP:  Section 5.5 (pp. 224 and ff.), Exercise 5.9
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);
double f(double x);    /* Function we're integrating */
double Trap(double a, double b, int n, int thread_count);
void Print_iters(int iters[], int n);

int main(int argc, char* argv[]) {
   double  global_result = 0.0;  /* Store result in global_result */
   double  a, b;                 /* Left and right endpoints      */
   int     n;                    /* Total number of trapezoids    */
   int     thread_count;

   if (argc != 2) Usage(argv[0]);
   thread_count = strtol(argv[1], NULL, 10);
   printf("Enter a, b, and n\n");
   scanf("%lf %lf %d", &a, &b, &n);

   global_result = Trap(a, b, n, thread_count);

   printf("With n = %d trapezoids, our estimate\n", n);
   printf("of the integral from %f to %f = %.14e\n",
      a, b, global_result);
   return 0;
}  /* main */

/*--------------------------------------------------------------------
 * Function:    Usage
 * Purpose:     Print command line for function and terminate
 * In arg:      prog_name
 */
void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
   exit(0);
}  /* Usage */

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input arg:   x
 * Return val:  f(x)
 */
double f(double x) {
   double return_val;

   return_val = x*x;
   return return_val;
}  /* f */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Use trapezoidal rule to estimate definite integral
 * Input args:  
 *    a: left endpoint
 *    b: right endpoint
 *    n: number of trapezoids
 * Return val:
 *    approx:  estimate of integral from a to b of f(x)
 */
double Trap(double a, double b, int n, int thread_count) {
   double  h, approx;
   int  i;
   int* iters = (int *) malloc(n*sizeof(int));
   
   h = (b-a)/n; 
   approx = (f(a) + f(b))/2.0; 
#  pragma omp parallel for num_threads(thread_count) \
      reduction(+: approx) schedule(runtime)
   for (i = 1; i <= n-1; i++) {
      approx += f(a + i*h);
      iters[i] = omp_get_thread_num();
   }
   approx = h*approx; 

   Print_iters(iters, n);
   
   free(iters);
   return approx;
}  /* Trap */


/*------------------------------------------------------------------
 * Function:  Print_iters
 * Purpose:   Print the iters array generated by the Trap function
 * In args:   iters:  iters[i] = thread that executed iteration i
 *            n:      order of the iters array
 */
void Print_iters(int iters[], int n) {
   int i, curr_thread;

   printf("iters = ");
   for (i = 0; i < n; i++)
      printf("%d ", iters[i]);
   printf("\n\n");

   printf("Thread\tIterations\n");

   printf("  %d  \t  0 -- ", iters[0]);
   curr_thread = iters[0];

   for (i = 1; i < n-1; i++) 
      if (curr_thread != iters[i]) {
         printf("%d\n", i-1);
         printf("  %d  \t  %d -- ", iters[i], i);
         curr_thread = iters[i];
      }

   if (curr_thread == iters[n-1]) {
      printf("%d\n", n-1);
   } else {
      printf("%d\n", n-2);
      printf("  %d  \t  %d -- %d\n", iters[n-1], n-1, n-1);
   }
   printf("\n");
}  /* Print_iters */
