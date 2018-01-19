/*------------------------------------------------------------------------------
* rtklib.h : rtk positioning common libraries
*
*          Copyright (C) 2007 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.8 $ $Date: 2007/03/20 08:44:46 $
* history : 2007/01/13 1.0 new
*           2007/03/20 1.1 add member rb of type sol_t
*-----------------------------------------------------------------------------*/
#ifndef RTKLIB_H
#define RTKLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef __cplusplus
extern "C" {
#endif

/* constants -----------------------------------------------------------------*/

#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define MU          3.986005E14         /* earth gravitational constant (ICD-GPS) */
#define OMGE        7.2921151467E-5     /* earth angular velocity (ICD-GPS) (rad/s) */
#define RE          6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE          (1.0/298.257223563) /* earth flattening (WGS84) */
#define FREQ1       1.57542E9           /* GPS L1 frequency (Hz) */
#define FREQ2       1.22760E9           /* GPS L2 frequency (Hz) */

#define NFREQ       2                   /* number of frequencies */
#define MAXSAT      32                  /* max satellite number (1-MAXSAT) */
#define MAXOBS      32                  /* max number of obs in an epoch */
#define MAXDT       0.05                /* max difference of epoch (sec) */

extern const double chisqr[];           /* chi-sqr(n) table (alpha=0.001) */

/* types ---------------------------------------------------------------------*/

typedef struct {        /* time struct */
    time_t time;        /* time (sec) expressed by standard time_t */
    double sec;         /* fraction of sec under 1 sec */
} gtime_t;

typedef struct {        /* observation data record */
    gtime_t time;       /* receiver sampling time */
    int sat,rcv;        /* satellite/receiver number */
    double L[NFREQ];    /* observation data carrier-phase (cycle) */
    double P[NFREQ];    /* observation data pseudorange (m) */
    short LLI[NFREQ];   /* loss of lock indicator */
} obsd_t;

typedef struct {        /* observation data */
    int n;              /* number of obervation data records */
    obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* satellite ephemeris and clock parameters */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva,svh;        /* sv accuracy, sv health */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot,crc,crs,cuc,cus,cic,cis,toes;
                        /* sv ephemeris parameters */
    double f0,f1,f2,tgd; /* sv clock parameters */
} eph_t;

typedef struct {        /* navigation messages */
    int n;              /* number of ephemeris and clock parameters */
    eph_t *eph;         /* satellite ephemeris and clock parameters */
    double ion[8];      /* iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double utc[4];      /* delta-utc parameters */
} nav_t;

typedef struct {        /* positioning solution */
    double rr[3];       /* estimated receiver position (ecef) (m) */
    double Qr[9];       /* estimated receiver positoin covariance */
    double rb[3];       /* reference station position (ecef) (m) */
} sol_t;

/* callback function types */
typedef int (* infunc_t)(obsd_t **, nav_t **);
typedef void (* outfunc_t)(gtime_t, sol_t, sol_t, int, int);

/* matrix and vector functions -----------------------------------------------*/
extern double *mat(int n, int m);
extern double *zeros(int n, int m);
extern double *eye(int n);
extern double dot(const double *a, const double *b, int n);
extern double norm(const double *a, int n);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern int matinv(double *A, int n);
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X);
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q);
extern int filter(const double *x, const double *P, double *H, double *v,
                  const double *R, int n, int m, double *xp, double *Pp);
extern void matprint(const double *A, int n, int m, int p, int q);

/* time and string functions -------------------------------------------------*/
extern double str2num(const char *s, int i, int n);
extern int str2time(const char *s, int i, int n, gtime_t *t);
extern void time2str(gtime_t t, char *str, int n);
extern gtime_t epoch2time(const double *ep);
extern void time2epoch(gtime_t t, double *ep);
extern gtime_t gpst2time(int week, double sec);
extern double time2gpst(gtime_t t, int *week);
extern gtime_t timeadd(gtime_t t, double sec);
extern double timediff(gtime_t t1, gtime_t t2);
extern gtime_t gpst2utc(gtime_t t);
extern gtime_t utc2gpst(gtime_t t);

/* coordinates functions -----------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);
extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);
extern void covenu(const double *pos, const double *P, double *Q);
extern double geoidh(const double *pos);

/* input/output functions ----------------------------------------------------*/
extern int readrnx(char **files, int n, obs_t *obs, nav_t *nav);

/* integer least-square estimation -------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s);

/* navigation functions ------------------------------------------------------*/
extern void eph2pos(gtime_t t, const eph_t *eph, double pr, double *rs,
                    double *dts);
extern void satpos(const obsd_t *obs, int n, const nav_t *nav, double *rs,
                   double *dts);
extern void satazel(const double *pos, const double *e, double *azel);
extern double geodist(const double *rs, const double *rr, double *e);
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel);
extern double tropmodel(const double *pos, const double *azel);
extern int pntpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *ion, double elmin, double *rr, double *Qr,
                  double *dtr, double *azel);
extern void rtkpos(infunc_t input, outfunc_t output, int mode, int nf,
                   const double *rb, double elmin, double thres, int opt);

#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */
