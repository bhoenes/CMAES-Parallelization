/* --------------------------------------------------------- */
/* --- File: cmaes.h ----------- Author: Nikolaus Hansen --- */
/* ---------------------- last modified: IX 2010         --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003-2010  Nikolaus Hansen. 
     e-mail: nikolaus.hansen (you know what) inria.fr
      
     License: see file cmaes.c
   
*/
#ifndef NH_cmaes_h /* only include ones */ 
#define NH_cmaes_h 

#include <time.h>

typedef struct 
/* random_t 
 * sets up a pseudo random number generator instance 
 */
{
  /* Variables for Uniform() */
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int *rgrand;
  
  /* Variables for Gauss() */
  short flgstored;
  double hold;
} random_t;

typedef struct 
/* timings_t 
 * time measurement, used to time eigendecomposition 
 */
{
  /* for outside use */
  double totaltime; /* zeroed by calling re-calling timings_start */
  double totaltotaltime;
  double tictoctime; 
  double lasttictoctime;
  
  /* local fields */
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted; 

  double lastdiff;
  double tictoczwischensumme;
} timings_t;

typedef struct 
/* readpara_t
 * collects all parameters, in particular those that are read from 
 * a file before to start. This should split in future? 
 */
{
  /* input parameter */
  int N; /* problem dimension, must stay constant */
  unsigned int seed; 
  double * xstart; 
  double * typicalX; 
  int typicalXcase;
  double * rgInitialStds;
  double * rgDiffMinChange; 

  /* termination parameters */
  double stopMaxFunEvals; 
  double facmaxeval;
  double stopMaxIter; 
  struct { int flg; double val; } stStopFitness; 
  double stopTolFun;
  double stopTolFunHist;
  double stopTolX;
  double stopTolUpXFactor;

  /* internal evolution strategy parameters */
  int lambda;          /* -> mu, <- N */
  int mu;              /* -> weights, (lambda) */
  double mucov, mueff; /* <- weights */
  double *weights;     /* <- mu, -> mueff, mucov, ccov */
  double damps;        /* <- cs, maxeval, lambda */
  double cs;           /* -> damps, <- N */
  double ccumcov;      /* <- N */
  double ccov;         /* <- mucov, <- N */
  double diagonalCov;  /* number of initial iterations */
  struct { int flgalways; double modulo; double maxtime; } updateCmode;
  double facupdateCmode;

  /* supplementary variables */

  char *weigkey; 
  char resumefile[99];
  char **rgsformat;
  void **rgpadr;
  char **rgskeyar;
  double ***rgp2adr;
  int n1para, n1outpara;
  int n2para;
} readpara_t;

typedef struct 
/* cmaes_t 
 * CMA-ES "object" 
 */
{
  char *version;
  readpara_t sp;
  random_t rand; /* random number generator */

  double sigma;  /* step size */

  double *rgxmean;  /* mean x vector, "parent" */
  double *rgxbestever; 
  double **rgrgx;   /* range of x-vectors, lambda offspring */
  int *index;       /* sorting index of sample pop. */
  double *arFuncValueHist;

  short flgIniphase; /* not really in use anymore */
  short flgStop; 

  double chiN; 
  double **C;  /* lower triangular matrix: i>=j for C[i][j] */
  double **B;  /* matrix with normalize eigenvectors in columns */
  double *rgD; /* axis lengths */

  double *rgpc;
  double *rgps;
  double *rgxold; 
  double *rgout; 
  double *rgBDz;   /* for B*D*z */
  double *rgdTmp;  /* temporary (random) vector used in different places */
  double *rgFuncValue; 
  double *publicFitness; /* returned by cmaes_init() */

  double gen; /* Generation number */
  double countevals;
  double state; /* 1 == sampled, 2 == not in use anymore, 3 == updated */

  double maxdiagC; /* repeatedly used for output */
  double mindiagC;
  double maxEW;
  double minEW;

  char sOutString[330]; /* 4x80 */

  short flgEigensysIsUptodate;
  short flgCheckEigen; /* control via signals.par */
  double genOfEigensysUpdate; 
  timings_t eigenTimings;
 
  double dMaxSignifKond; 				     
  double dLastMinEWgroesserNull;

  short flgresumedone; 

  time_t printtime; 
  time_t writetime; /* ideally should keep track for each output file */
  time_t firstwritetime;
  time_t firstprinttime; 

} cmaes_t; 

//Doug Hakkarinen
//Moving the methods for objects to here to be able to access them.

long   random_init(random_t *, long unsigned seed /* 0==clock */);
void   random_exit(random_t *);
double random_Gauss(random_t *); /* (0,1)-normally distributed */
double random_Uniform(random_t *);
long   random_Start(random_t *, long unsigned seed /* 0==1 */);

void   timings_init(timings_t *timing);
void   timings_start(timings_t *timing); /* fields totaltime and tictoctime */
double timings_update(timings_t *timing);
void   timings_tic(timings_t *timing);
double timings_toc(timings_t *timing);

void readpara_init (readpara_t *, int dim, int seed,  const double * xstart, 
                    const double * sigma, int lambda, const char * filename);
void readpara_exit(readpara_t *);
void readpara_ReadFromFile(readpara_t *, const char *szFileName);
void readpara_SupplementDefaults(readpara_t *);
void readpara_SetWeights(readpara_t *, const char * mode);
void readpara_WriteToFile(readpara_t *, const char *filenamedest, 
                          const char *parafilesource);

double const * cmaes_SetMean(cmaes_t *, const double *xmean);
double * cmaes_PerturbSolutionInto(cmaes_t *t, double *xout, 
                                   double const *xin, double eps);
void cmaes_WriteToFile(cmaes_t *, const char *key, const char *name);
void cmaes_WriteToFileAW(cmaes_t *t, const char *key, const char *name, 
                         char * append);
void cmaes_WriteToFilePtr(cmaes_t *, const char *key, FILE *fp);
void cmaes_ReadFromFilePtr(cmaes_t *, FILE *fp); 
void cmaes_FATAL(char const *s1, char const *s2, 
                 char const *s3, char const *s4);

/* ------------------- Locally visibly ----------------------- */

char * getTimeStr(void); 
void TestMinStdDevs( cmaes_t *);
/* static void WriteMaxErrorInfo( cmaes_t *); */

void Eigen( int N,  double **C, double *diag, double **Q, 
                   double *rgtmp);
int  Check_Eigen( int N,  double **C, double *diag, double **Q);
void QLalgo2 (int n, double *d, double *e, double **V); 
void Householder2(int n, double **V, double *d, double *e); 
void Adapt_C2(cmaes_t *t, int hsig);

void FATAL(char const *sz1, char const *s2, 
                  char const *s3, char const *s4);
void ERRORMESSAGE(char const *sz1, char const *s2, 
                         char const *s3, char const *s4);
void   Sorted_index( const double *rgFunVal, int *index, int n);
int    SignOfDiff( const void *d1, const void * d2);
double douSquare(double);
double rgdouMax( const double *rgd, int len);
double rgdouMin( const double *rgd, int len);
double douMax( double d1, double d2);
double douMin( double d1, double d2);
int    intMin( int i, int j);
int    MaxIdx( const double *rgd, int len);
int    MinIdx( const double *rgd, int len);
double myhypot(double a, double b);
double * new_double( int n);
void * new_void( int n, size_t size); 


#endif 
