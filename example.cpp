#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>

double x(int i);                    
double y(int i);
double exactSolution(double x, double y);    
double S(double x, double y);                 

double a=0.0, double b=5.0, double c = 0.0, double d = 5.0;
const int m=100;                         
const int n=100;

double tolerance=1E-15;
int maxIterations=1000000;

double dx=(b-a)/(m-1);
double dy=(d-c)/(n-1);

double Un  [m][n];
double Unp1[m][n];

int main()
{
    
    for(int i=1; i<m-1; i++) {
        for (int j = 1; j < n-1; j++) {
            Un[i][j]=0.0;
        }

    }
    
    for (int i = 0; i < m; i++) {
        Un[i][0] = exactSolution(x(i), y(0));
        Un[i][n -1] = exactSolution(x(i), y(n - 1));
        Unp1[i][0] = exactSolution(x(i), y(0));
        Unp1[i][n -1] = exactSolution(x(i), y(n - 1));
    }
    for (int i = 0; i < n; i++) {
        Un[0][i] = exactSolution(x(0), y(i));
        Un[m - 1][i] = exactSolution(x(m-1), y(i));
        Unp1[0][i] = exactSolution(x(0), y(i));
        Unp1[m - 1][i] = exactSolution(x(m-1), y(i));
    }
    int iterations=0;
    double iterationError = 1.;
    while(iterationError > tolerance && iterations < maxIterations){
        iterations++; // 
        if(iteration_count % 1000 == 0) std::cout<<"iteration " << iteration_count << std::endl;
        for(int i=1; i< m-1; i++){
            for (int j = 1; j < n -1; j++) {

                Unp1[i]=(dy*dy*dx*dx*S(x(i), y(j)) - \
                         dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                         dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
        iterationError=0.0;
        // Calculate diff between Un, Up+1
        for(int i = 0; i< m; i++){
            for (in j = 0; j < n; j++) {
                double localError = fabs(Unp1[i][j] - Un[i][j]);
                if (localError > iterationError)
                    iterationError = localError;
             }
        }

        for(int i=0; i < m; i++){
            for (int j = 0; j < n; j++) { 
                Un[i][j] = Unp1[i][j];
            }
        }

        if(iteration_count % 1000 == 0) std::cout<< "The error between two iterates is " << iterationError << std::endl;
    }

    // Compute the maximum error between the computed and exact solutions:
    double solution_error=0.0;
    for(int i = 0; i < m; i++){
        for (int j = 0; j < n; j++) {
            double localError=fabs(Unp1[i][j] - exactSolution(x(i), y(j)));
            if (localError > solution_error)
               solution_error = localError;
        }
    }

    // Output:
    std::cout                                                              << std::endl << std::endl;
    std::cout<< "-------------------------------------------------------"               << std::endl;
    std::cout<< "SUMMARY:"                                                 << std::endl << std::endl;
    std::cout<< "The error between two iterates is "    << iterationError << std::endl << std::endl;
    std::cout<< "The maximum error in the solution is " << solution_error               << std::endl;
    std::cout<< "-------------------------------------------------------"  << std::endl << std::endl;
    return 0;
}
double exactSolution(double x, double y)
    return sin(x) + sin(y);

double x(int i)
    return a+i*dx;
double y(int n)
    return c + i*dx;
double S(double x, double y)
        return -sin(x) - sin(y);
