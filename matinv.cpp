//////////////////////////////////////////////////////////////////////////////////////
//
// http://www.graphviz.org/pub/graphviz/development/doxygen/html/matinv_8c_source.html
//
//////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>


//#include "global_extern.h"

#define N_NEW(n,t)       (t*)calloc((n),sizeof(t))

extern int lu_decompose(int n);
//extern void matinv(int n,double B[][128]);
extern void lu_solve(double *x, double *b, int n);

static double *scales;
static double **lu;
static int *ps;

//double B[25][25];
//double Binv[25][25] = {0};
double A[128][128]={0};
double Ainv[128][128]={0};

void free_array(double **rv)
{
        if (rv) {
                free(rv[0]);
                free(rv);
        }
}

double **new_array(int m, int n, double ival)
{
        double **rv;
        double *mem;
        int i, j;
        rv = N_NEW(m, double *);
        mem = N_NEW(m * n, double);
        for (i = 0; i < m; i++) {
                rv[i] = mem;
                mem = mem + n;
                for (j = 0; j < n; j++)
                        rv[i][j] = ival;
        }
        return rv;
}

/*void matinv()
{
        int i,j,k;
        char c,buff[100];
        FILE *fp;

        if((fp  = fopen("dgps.csv", "r" )) ==NULL ){
               printf( "nofile\n" );exit(1);
        }

        for(j=1; j<=25; j++){
                for(k=1; k<=25; k++){
                        i = 0;
                        while((c = fgetc(fp)) != ',' && c != '\n' && c != EOF)
                        {
                                buff[i++] = c;
                        }
                        buff[i] = '\0';

                        B[j-1][k-1] = atof(buff);
                }//k
        }//j=1 -> j=?

        matinv(25);
        return (0);
}
//*/

int matinv(int n,double B[])
{
	
        register int i, j;
        int k,l;
        double *b, temp;
//      double A[128][128]={0},Ainv[128][128]={0};

        for(k=0; k<n; k++){
                for(l=0; l<n; l++){
                        A[k][l] = B[k*n+l];
                }
        }
	
	
        /* Decompose matrix into L and U triangular matrices */

        if (lu_decompose(n) == 0)
                return (0); /* Singular */

        /* Invert matrix by solving n simultaneous equations n times */
        b = N_NEW(n, double);
        for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++)
                        b[j] = 0.0;
                b[i] = 1.0;
                lu_solve(Ainv[i], b, n); /* Into a row of Ainv: fix later */
        }
        free(b);

        /* Transpose matrix */
        for (i = 0; i < n; i++) {
                for (j = 0; j < i; j++) {
                        temp = Ainv[i][j];
                        Ainv[i][j] = Ainv[j][i];
                        Ainv[j][i] = temp;
                }
        }
		
        for(k=0; k<n; k++){
                for(l=0; l<n; l++){
                       B[k*n+l]=Ainv[k][l]; 
                }
        }
	        return (0);
}

 /* lu_decompose() decomposes the coefficient matrix A into upper and lower
* triangular matrices, the composite being the LU matrix.
*
* The arguments are:
*
* a - the (n x n) coefficient matrix
* n - the order of the matrix
*
* 1 is returned if the decomposition was successful,
* and 0 is returned if the coefficient matrix is singular.
*/

int lu_decompose(int n)
{
        register int i, j, k;
        int l,m;
        int pivotindex = 0;
        double pivot, biggest, mult, tempf;
        double a[128][128]={0};

        for(l=0; l<n; l++){
                for(m=0; m<n; m++){
                        a[l][m] = A[l][m];
                }
        }

        if (lu)
                free_array(lu);
        lu = new_array(n, n, 0.0);

        if (ps)
                free(ps);
        ps = N_NEW(n, int);

        if (scales)
                free(scales);
        scales = N_NEW(n, double);

        for (i = 0; i < n; i++) { /* For each row */
                /* Find the largest element in each row for row equilibration */
                biggest = 0.0;

                for (j = 0; j < n; j++)
                        if (biggest < (tempf = fabs(lu[i][j] = a[i][j])))
                                biggest = tempf;

                if (biggest != 0.0)
                        scales[i] = 1.0 / biggest;
                else {
                        scales[i] = 0.0;
                        return (0); /* Zero row: singular matrix */
                }
                ps[i] = i; /* Initialize pivot sequence */
        }

        for (k = 0; k < n - 1; k++) { /* For each column */
                /* Find the largest element in each column to pivot around */
                biggest = 0.0;

                for (i = k; i < n; i++) {
                        if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
                                biggest = tempf;
                                pivotindex = i;
                        }
                }

                if (biggest == 0.0)
                        return (0); /* Zero column: singular matrix */
                if (pivotindex != k) { /* Update pivot sequence */
                        j = ps[k];
                        ps[k] = ps[pivotindex];
                        ps[pivotindex] = j;
                }

                /* Pivot, eliminating an extra variable each time */
                pivot = lu[ps[k]][k];

                for (i = k + 1; i < n; i++) {
                        lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;

                        if (mult != 0.0) {
                                for (j = k + 1; j < n; j++)
                                        lu[ps[i]][j] -= mult * lu[ps[k]][j];
                        }
                }
        }

        if (lu[ps[n - 1]][n - 1] == 0.0)
                return (0); /* Singular matrix */

        return (1);
}

/* lu_solve() solves the linear equation (Ax = b) after the matrix A has
* been decomposed with lu_decompose() into the lower and upper triangular
* matrices L and U.
*
* The arguments are:
*
* x - the solution vector
* b - the constant vector
* n - the order of the equation
*/

void lu_solve(double *x, double *b, int n)
{
        register int i, j;
        double dot;
        /* Vector reduction using U triangular matrix */

        for (i = 0; i < n; i++) {
                dot = 0.0;
                for (j = 0; j < i; j++)
                        dot += lu[ps[i]][j] * x[j];
                x[i] = b[ps[i]] - dot;
        }

        /* Back substitution, in L triangular matrix */
        for (i = n - 1; i >= 0; i--) {
                dot = 0.0;
                for (j = i + 1; j < n; j++)
                        dot += lu[ps[i]][j] * x[j];
                x[i] = (x[i] - dot) / lu[ps[i]][i];
        }
}
