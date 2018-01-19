/*------------------------------------------------------------------------------
* rtkcmn.c : rtk common functions
*
*          Copyright (C) 2007 by T.TAKASU, All rights reserved.
*
* options : -DMKL      use intel MKL
*           -DNOLAPACK use embedded matrix routines
*
* version : $Revision: 1.2 $ $Date: 2007/01/23 16:28:05 $
* history : 2007/01/12 1.0 new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#ifdef MKL
#include <mkl_types.h>
#include <mkl_lapack.h>
#include <mkl_blas.h>
#endif

static const char rcsid[]="$Id: rtkcmn.c,v 1.2 2007/01/23 16:28:05 ttaka Exp $";

/* constants -----------------------------------------------------------------*/

#define MAXITR      10          /* max number of iteration for point pos */
#define MAXDTOE     10800.0     /* max time difference to ephemeris toe (sec) */
#define ERRA        1.5         /* measurement error factor a of code */
#define ERRB        1.0         /* measurement error factor b of code */
#define ERRI        0.3         /* ionospheric model error factor */

const static double gpst0[]={1980,1,6,0,0,0}; /* gpsweek 0 - day 1 */

const static double leaps[][7]={ /* leap seconds {y,m,d,h,m,s,utc-gpst,...} */
    {2006,1,1,0,0,0,-14},
    {1999,1,1,0,0,0,-13},
    {1997,7,1,0,0,0,-12},
    {1996,1,1,0,0,0,-11},
    {1994,7,1,0,0,0,-10},
    {1993,7,1,0,0,0, -9},
    {1992,7,1,0,0,0, -8},
    {1991,1,1,0,0,0, -7},
    {1990,1,1,0,0,0, -6},
    {1988,1,1,0,0,0, -5},
    {1985,7,1,0,0,0, -4},
    {1983,7,1,0,0,0, -3},
    {1982,7,1,0,0,0, -2},
    {1981,7,1,0,0,0, -1}
};
const double chisqr[45]={        /* chi-sqr(n) (alpha=0.001) */
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,31.3,32.9,34.5,36.1,37.7,
    39.3,40.8,42.3,43.8,45.3,46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,70.1,73.4,86.7,99.6,112.3,124.8,137.2
};
/* lapack/blas prototypes ----------------------------------------------------*/
#ifdef MKL
#define dgemm_      dgemm
#define dgetrf_     dgetrf
#define dgetri_     dgetri
#define dgetrs_     dgetrs
#else
extern void dgemm_(char *, char *, int *, int *, int *, double *, double *,
                   int *, double *, int *, double *, double *, int *);
extern void dgetrf_(int *, int *, double *, int *, int *, int *);
extern void dgetri_(int *, double *, int *, int *, double *, int *, int *);
extern void dgetrs_(char *, int *, int *, double *, int *, int *, double *,
                    int *, int *);
#endif
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
    double *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
		int temp=0;
    //    fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
    int *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
		int temp=0;
      //  fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *zeros(int n, int m)
{
    double *p;
    
#if NOCALLOC
    if ((p=mat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0.0;
#else
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(sizeof(double),n*m))) {
		int temp=0;
       // fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
#endif
    return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *eye(int n)
{
    double *p;
    int i;
    
    if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors 
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void cross3(const double *a, const double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args   : double *a        I   vector a (3 x 1)
*          double *b        O   normlized vector (3 x 1) || b || = 1
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int normv3(const double *a, double *b)
{
    double r;
    if ((r=norm(a,3))<=0.0) return 0;
    b[0]=a[0]/r;
    b[1]=a[1]/r;
    b[2]=a[2]/r;
    return 1;
}
/* copy matrix -----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
/* matrix routines -----------------------------------------------------------*/

#ifdef LAPACK /* with LAPACK/BLAS or MKL */

/* multiply matrix (wrapper of blas dgemm) -------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C)
* args   : char   *tr       I  transpose flags ("N":normal,"T":transpose)
*          int    n,k,m     I  size of (transposed) matrix A,B
*          double alpha     I  alpha
*          double *A,*B     I  (transposed) matrix A (n x m), B (m x k)
*          double beta      I  beta
*          double *C        IO matrix C (n x k)
* return : none
*-----------------------------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;
    
    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
           &ldb,&beta,C,&n);
}
/* inverse of matrix -----------------------------------------------------------
* inverse of matrix (A=A^-1)
* args   : double *A        IO  matrix (n x n)
*          int    n         I   size of matrix A
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double *work;
    int info,lwork=n*16,*ipiv=imat(n,1);
    
    work=mat(lwork,1);
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if (!info) dgetri_(&n,A,&n,ipiv,work,&lwork,&info);
    free(ipiv); free(work);
    return info;
}
/* solve linear equation -------------------------------------------------------
* solve linear equation (X=A\Y or X=A'\Y)
* args   : char   *tr       I   transpose flag ("N":normal,"T":transpose)
*          double *A        I   input matrix A (n x n)
*          double *Y        I   input matrix Y (n x m)
*          int    n,m       I   size of matrix A,Y
*          double *X        O   X=A\Y or X=A'\Y (n x m)
* return : status (0:ok,0>:error)
* notes  : matirix stored by column-major order (fortran convention)
*          X can be same as Y
*-----------------------------------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info,*ipiv=imat(n,1);
    
    matcpy(B,A,n,n);
    matcpy(X,Y,n,m);
    dgetrf_(&n,&n,B,&n,ipiv,&info);
    if (!info) dgetrs_((char *)tr,&n,&m,B,&n,ipiv,X,&n,&info);
    free(ipiv); free(B); 
    return info;
}

#else /* without LAPACK/BLAS or MKL */

/* multiply matrix -----------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);
    
    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;
    
    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
            }
            *d=-(*d); vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) {free(vv); return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;
    
    for (i=0;i<n;i++) {
        ip=indx[i]; s=b[ip]; b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
    }
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double d,*B;
    int i,j,*indx;
    
    indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
    if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0; A[j+j*n]=1.0;
        lubksb(B,n,indx,A+j*n);
    }
    free(indx); free(B);
    return 0;
}
/* solve linear equation -----------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info;
    
    matcpy(B,A,n,n);
    if (!(info=matinv(B,n))) matmul(tr[0]=='N'?"NN":"TN",n,m,n,1.0,B,Y,0.0,X);
    free(B);
    return info;
}
#endif
/* end of matrix routines ----------------------------------------------------*/

/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (x=(A*A')^-1*A*y)
* args   : double *A        I   transpose of (weighted) design matrix (n x m)
*          double *y        I   (weighted) measurements (m x 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *x        O   estmated parameters (n x 1)
*          double *Q        O   esimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : for weighted least square, replace A and y to A*w and w*y (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q)
{
    int info;
    double *Ay;
    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}
/* kalman filter ---------------------------------------------------------------
* kalman filter state update (K=P*H*(H'*P*H+R)^-1,xp=x+K*v,Pp=(I-K*H')*P)
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   updated states vector (n x 1)
*          double *Pp       O   updated covariance matrix of states (n x n)
* return : status (0:ok,0>:error)
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int filter(const double *x, const double *P, double *H, double *v,
                  const double *R, int n, int m, double *xp, double *Pp)
{
    int info;
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    memcpy(Q,R,sizeof(double)*m*m);
    memcpy(xp,x,sizeof(double)*n);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
    free(F); free(Q); free(K); free(I);
    return info;
}
/* print matrix ----------------------------------------------------------------
* print matrix to stdout
* args   : double *A        I   matrix (n x m)
*          int    n,m       I   number of rows and columns of A
*          int    p,q       I   total columns, columns under decimal point
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void matprint(const double *A, int n, int m, int p, int q)
{
    int i,j;
    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) printf("%*.*f ",p,q,A[i+j*n]);
        printf("\n");
    }
}
/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
*-----------------------------------------------------------------------------*/
extern double str2num(const char *s, int i, int n)
{
    char str[256],*p=str;
    double value;
    if (i<0||(int)strlen(s)<i||sizeof(str)-1<i) return 0.0;
    for (s+=i;*s&&--n>=0;s++) *p++=*s=='D'?'E':*s; *p='\0';
    return sscanf(str,"%lf",&value)==1?value:0.0;
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int str2time(const char *s, int i, int n, gtime_t *t)
{
    char str[256],*p=str;
    double ep[6];
    if (i<0||(int)strlen(s)<i||sizeof(str)-1<i) return -1;
    for (s+=i;*s&&--n>=0;) *p++=*s++; *p='\0';
    if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
        return -1;
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    *t=epoch2time(ep);
    return 0;
}
/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t epoch2time(const double *ep)
{
    struct tm t={0};
    gtime_t time;
    t.tm_year=(int)ep[0]-1900; t.tm_mon=(int)ep[1]-1; t.tm_mday=(int)ep[2];
    t.tm_hour=(int)ep[3]; t.tm_min=(int)ep[4]; t.tm_sec=(int)ep[5];
    time.time=mktime(&t); time.sec=ep[5]-t.tm_sec;
    return time;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
    struct tm *tt=localtime(&t.time);
    ep[0]=tt->tm_year+1900; ep[1]=tt->tm_mon+1; ep[2]=tt->tm_mday;
    ep[3]=tt->tm_hour; ep[4]=tt->tm_min; ep[5]=tt->tm_sec+t.sec;
}
/* gpstime to time -------------------------------------------------------------
* convert gps week and tow to gtime_t struct
* args   : int    week      I   gps week number
*          double sec       I   gpstime (time of week) (sec)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2time(int week, double sec)
{
    gtime_t t=epoch2time(gpst0);
    t.time+=86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}
/* time to gpstime -------------------------------------------------------------
* convert gtime_t struct to gps week and tow
* args   : gtime_t t        I   gtime_t struct
*          int    *week     O   gps week number
* return : gpstime (time of week) (sec)
*-----------------------------------------------------------------------------*/
extern double time2gpst(gtime_t t, int *week)
{
    gtime_t t0=epoch2time(gpst0);
    time_t sec=t.time-t0.time;
    *week=(int)(sec/(86400*7));
    return (double)(sec-(long)*week*86400*7)+t.sec;
}
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (sec)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
extern gtime_t timeadd(gtime_t t, double sec)
{
    double tt;
    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (sec)
*-----------------------------------------------------------------------------*/
extern double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}
/* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2utc(gtime_t t)
{
    int i;
    gtime_t tu;
    for (i=0;i<sizeof(leaps)/sizeof(*leaps);i++) {
        tu=timeadd(t,leaps[i][6]);
        if (timediff(tu,epoch2time(leaps[i]))>=0.0) return tu;
    }
    return t;
}
/* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
*-----------------------------------------------------------------------------*/
extern gtime_t utc2gpst(gtime_t t)
{
    int i;
    for (i=0;i<sizeof(leaps)/sizeof(*leaps);i++) {
        if (timediff(t,epoch2time(leaps[i]))>=0.0) return timeadd(t,-leaps[i][6]);
    }
    return t;
}
/* time to string --------------------------------------------------------------
* convert gtime_t struct to string
* args   : gtime_t t        I   gtime_t struct
*          char   *s        O   string ("yyyy/mm/dd hh:mm:ss.s")
*          int    n         I   columns of sec under decimal point
* return : none
*-----------------------------------------------------------------------------*/
extern void time2str(gtime_t t, char *s, int n)
{
    double ep[6];
    if (n>=0&&1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE*(2.0-FE),r2=r[0]*r[0]+r[1]*r[1],z,zk,v=RE,sinp;
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=atan2(r[1],r[0]);
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void pos2ecef(const double *pos, double *r)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    double e2=FE*(2.0-FE),v=RE/sqrt(1.0-e2*sinp*sinp);
    r[0]=(v+pos[2])*cosp*cosl;
    r[1]=(v+pos[2])*cosp*sinl;
    r[2]=(v*(1.0-e2)+pos[2])*sinp;
}
/* ecef to local coordinate transfromation matrix ----------------------------*/
static void tolocal(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}
/* transform coordinate between ecef and local (dir=0:to ecef,1:to local) ----*/
static void eceflocal(const double *pos, const double *a, double *b, int dir)
{
    double E[9];
    tolocal(pos,E);
    matmul(dir?"NN":"TN",3,1,3,1.0,E,a,0.0,b);
}
/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2enu(const double *pos, const double *r, double *e)
{
    eceflocal(pos,r,e,1);
}
/* transform local vector to ecef coordinate -----------------------------------
* transform local tangental coordinate vector to ecef
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *e        I   vector in local tangental coordinate {e,n,u}
*          double *r        O   vector in ecef coordinate {x,y,z}
* return : none
*-----------------------------------------------------------------------------*/
extern void enu2ecef(const double *pos, const double *e, double *r)
{
    eceflocal(pos,e,r,0);
}
/* transform covariance to local tangental coordinate --------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *P        I   covariance in ecef coordinate
*          double *Q        O   covariance in local tangental coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covenu(const double *pos, const double *P, double *Q)
{
    double E[9],EP[9];
    tolocal(pos,E);
    matmul("NN",3,3,3,1.0,E,P,0.0,EP);
    matmul("NT",3,3,3,1.0,EP,E,0.0,Q);
}
/* satellite ephemeris to satellite position/clock-bias ------------------------
* compute satellite position and clock-bias from satellite ephemeris
* args   : gtime_t t        I   time
*          eph_t *eph       I   satellite ephemeris and clock parameter
*          double pr        I   pseudorange measurement (m)
*          double *rs       O   satellilte positions (ecef) (m)
*          double *dts      O   satellilte clock-bias (sec)
* return : none
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t t, const eph_t *eph, double pr, double *rs,
                    double *dts)
{
    double tk,tc,tt,M,E,Ek,phi,sinE,cosE,sin2p,cos2p,r,i,x,y,OMG,sinO,cosO,cosi;
    
    tk=timediff(t,eph->toe); tc=timediff(t,eph->toc);
    tt=pr/CLIGHT+(eph->f0+eph->f1*tc+eph->f2*tc*tc);
    tk-=tt; tc-=tt; /* signal transmission time */
    M=eph->M0+(sqrt(MU/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    for (E=M,sinE=Ek=0.0;fabs(E-Ek)>1E-12;) {
        Ek=E; sinE=sin(Ek); E=M+eph->e*sinE;
    }
    cosE=cos(E);
    phi=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    sin2p=sin(2.0*phi); cos2p=cos(2.0*phi);
    phi+=eph->cus*sin2p+eph->cuc*cos2p;
    r=eph->A*(1.0-eph->e*cosE)+eph->crs*sin2p+eph->crc*cos2p;
    i=eph->i0+eph->cis*sin2p+eph->cic*cos2p+eph->idot*tk;
    OMG=eph->OMG0+(eph->OMGd-OMGE)*tk-OMGE*eph->toes;
    x=r*cos(phi); y=r*sin(phi);
    sinO=sin(OMG); cosO=cos(OMG); cosi=cos(i);
    rs[0]=x*cosO-y*cosi*sinO;
    rs[1]=x*sinO+y*cosi*cosO;
    rs[2]=y*sin(i);
    *dts=eph->f0+eph->f1*tc+eph->f2*tc*tc;
    *dts-=2.0*sqrt(MU*eph->A)*eph->e*sinE/CLIGHT/CLIGHT+eph->tgd;
}
/* satellite positions/clock-biases --------------------------------------------
* compute satellite positions and clock-biases from satellite ephemerides
* args   : obsd_t *obs      I   observation data records
*          int     n        I   number of observation data records
*          nav_t  *nav      I   navigation messages (sorted by transmssion time)
*          double *rs       O   satellilte positions (ecef) (m)
*          double *dts      O   satellilte clock-bias (sec)
* return : none
* notes  : rs[(0:2)+i*3],dts[i] = obs[i] satellite position/clock-bias
*          if no navigation messages, set 0.0 to rs[], dts[]
*-----------------------------------------------------------------------------*/
extern void satpos(const obsd_t *obs, int n, const nav_t *nav, double *rs,
                   double *dts)
{
    int i,j,k;
    double tt,tmin;
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        /* search ephemeris toe closest to obs time */
        for (j=0,k=-1,tmin=1E99;j<nav->n;j++) {
            if (nav->eph[j].sat!=obs[i].sat) continue;
            tt=fabs(timediff(nav->eph[j].toe,obs[i].time));
            if (nav->eph[j].svh!=0&&tt<MAXDTOE) {k=-1; break;}
            if (tt<tmin) {k=j; tmin=tt;}
        }
        if (k>=0) eph2pos(obs[i].time,nav->eph+k,obs[i].P[0],rs+i*3,dts+i);
        else rs[i*3]=rs[1+i*3]=rs[2+i*3]=dts[i]=0.0;
    }
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   receiver-to-satellilte unit vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
    int i,j;
    double r,rk,rot;
    
    if (norm(rs,3)<RE) return -1.0;
    e[2]=rs[2]-rr[2];
    for (i=0,r=2E7,rk=0.0;i<10;i++,rk=r) {
        rot=-OMGE*r/CLIGHT;
        e[0]=rs[0]-rs[1]*rot-rr[0]; /* e=Rz(-OMGE*r/C)*rs-rr */
        e[1]=rs[0]*rot+rs[1]-rr[1];
        r=norm(e,3);
        if (fabs(r-rk)<1E-4) {for (j=0;j<3;j++) e[j]/=r; return r;}
    }
    return -1.0;
}
/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     O   azimuth/elevation angle {az,el} (rad)
* return : none
*-----------------------------------------------------------------------------*/
extern void satazel(const double *pos, const double *e, double *azel)
{
    double enu[3];
    if (pos[2]<=-RE) {azel[0]=0.0; azel[1]=PI/2.0; return;}
    ecef2enu(pos,e,enu);
    azel[0]=atan2(enu[0],enu[1]);
    azel[1]=asin(enu[2]);
}
/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   gtime_t struct
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
{
    int week;
    double tt,f,psi,phi,lam,amp,per,x;
    
    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    psi=0.0137/(azel[1]/PI+0.11)-0.022;
    phi=pos[0]/PI+psi*cos(azel[0]);
    phi=phi<-0.416?-0.416:(phi>0.416?0.416:phi);
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);
    phi+=0.064*cos((lam-1.617)*PI);
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<0.0?0.0:amp; per=per<72000.0?72000.0:per;
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);
    tt=4.32E4*lam+time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */
    x=2.0*PI*(tt-50400.0)/per;
    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}
/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(const double *pos, const double *azel)
{
    double h,pres,temp,e,humi=50.0,tz;
    
    if (pos[2]<0||1E4<pos[2]||azel[1]<=0) return 0.0;
    h=pos[2]<0.0?0.0:pos[2];
    pres=1013.25*pow(1.0-2.2557E-5*h,5.257);
    temp=15.0-6.5E-3*h;
    e=6.108*pow(10.0,7.5*temp/(237.3+temp))*humi/100.0;
    tz=tan(PI/2.0-azel[1]);
    return 0.002277/sin(azel[1])*(pres+(1255.0/(temp+273.15)+0.05)*e-tz*tz);
}
/* pseudorange measurement error variance ------------------------------------*/
static double varerr(double el, double dion)
{
    double a2=ERRA*ERRA,b2=ERRB*ERRB,c2=ERRI*ERRI,sinel=sin(el);
    return a2+b2/sinel/sinel+c2*dion*dion;
}
/* pseudorange residuals -----------------------------------------------------*/
static int resduals(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *x, const double *ion,
                    double elmin, double *v, double *H, double *azel)
{
    int i,j,nv=0;
    double y,r,dion,sig,pos[3],e[3];
    
    ecef2pos(x,pos); /* receiver geodetic postion */ 
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        azel[i*2]=azel[1+i*2]=0.0;
        
        /* pseudorange/geometric distance/satellite azimuth/elevation angle */
        if ((y=obs[i].P[0])==0.0||(r=geodist(rs+i*3,x,e))<0.0) continue;
        satazel(pos,e,azel+i*2); if (azel[1+i*2]<elmin) continue;
        
        /* pseudorange model */
        dion=ionmodel(obs[i].time,ion,pos,azel+i*2);
        r+=x[3]-CLIGHT*dts[i]+dion+tropmodel(pos,azel+i*2);
        
        /* measurement error standard deviation (m) */
        sig=sqrt(varerr(azel[1+i*2],dion));
        
        /* design matrix (normalized) */
        for (j=0;j<4;j++) H[j+4*nv]=(j<3?-e[j]:1.0)/sig;
        
        /* residuals (normalized) */
        v[nv++]=(y-r)/sig; 
    }
    return nv;
}
/* single point positioning ----------------------------------------------------
* compute receiver position and clock-bias by single point positioning
* args   : obsd_t *obs      I   observation data records
*          int    n         I   number of observation data records
*          int    nf        I   number of frequency (1:L1,2:L1+L2)
*          double *rs       I   satellilte positions (ecef) (m)
*          double *dts      I   satellilte clock-biases (sec)
*          double *ion      I   broadcast ionosphere model parameters
*          double elmin     I   elevation cutoff angle (rad)
*          double *rr       O   receiver position (ecef) and clock-bias (m)
*          double *Qr       O   estimated position covarience (3 x 3)
*          double *dtr      O   estimated receiver clock-bias (sec)
*          double *azel     O   satellite azimuth/elevation {az,el,...} (rad)
* return : number of valid satellites(0>:error)
* notes  : rs[(0:2)+i*3],dts[i] = obs[i] satellite position/clock-bias
*          azel[(0:1)+i*2]      = obs[i] satellite azimuth/elevation angle
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *ion, double elmin, double *rr, double *Qr,
                  double *dtr, double *azel)
{
    int i,j,k,nv,info;
    double x[4],dx[4],Q[16],v[MAXOBS],H[4*MAXOBS];
    
    for (i=0;i<4;i++) x[i]=0.0; /* initial position/clock-bias */
    
    for (i=0;i<MAXITR;i++) {
        
        /* residuals */
        nv=resduals(obs,n,rs,dts,x,ion,elmin,v,H,azel);
        
        /* least square estimation */
        if ((info=lsq(H,v,4,nv,dx,Q))) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-4) {
            /* validation of solution */
            for (i=0;i<n;i++) if (azel[1+i*2]<0.0) return -1;
            if (nv>4&&dot(v,v,nv)>chisqr[nv-5]) return -1;
            
            for (j=0;j<3;j++) {
                rr[j]=x[j]; for (k=0;k<3;k++) Qr[j+k*3]=Q[j+k*4];
            }
            *dtr=x[3]/CLIGHT;
            return nv;
        }
    }
    return -1;
}
