/* --------------------------------------------------------- */
/*  This version splits the lambda among all the processors evenly, then finds the one that performed best, resets all the rest to have that mean, and start the next cycle.  As such, it enforces the parallelism at a cyclic level.
*/
/* --------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.h"
#include "forwardmodel.h"
#include "cmaes.h"
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "myconfig.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index, p, n) (((p)*((index)+1)-1)/(n))

double* reinit(cmaes_t * evo, int lambda, int numDipoles);
void resetSignals(cmaes_t *evo, int numDipoles);
int is_feasible(double * dipole, int N);


double fitfun(double const *x, int dim); 

/* the objective (fitness) function to be minimized */
double fitfun(double const *x, int N) {
  double * predictions;
  double rho;

  rho = 720;
 
  predictions = forwardModel(N,x,rho, 12.0, 16.0, 20.0);
  return objectiveFunction(32, predictions, observations);

}

/* the optimization loop */
int main(int argc, char **argv) {
  cmaes_t evo; /* an CMA-ES type struct or "object" */
  double *arFunvals,  *xfinal, *const*pop;
  int i,j;
  int numberDipoles;
  int id;  //Rank
  int p;  //Number processors
  double elapsed_time;//Time from beginning.
  double bestValue;
  int lambda;
  int maxLambda;
  int * sendCnts; //For MPI_Alltoallv for arFunVals
  int * sdispls;  //For MPI_Alltoallv for arFunVals
  int * recvCnts; //For MPI_Alltoallv for arFunVals
  int * rdispls; //For MPI_Alltoallv for arFunVals
  int * sendCntsPop; //For MPI_Alltoallv for pop
  int * sdisplsPop; //For MPI_Alltoallv for pop
  int * recvCntsPop; //For MPI_Alltoallv for pop
  int * rdisplsPop; //For MPI_Alltoallv for pop
  int canTerminate;
  int canTerminateBuffer;

  //Start MPI
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime(); //Set initial time.

  MPI_Comm_rank(MPI_COMM_WORLD, &id); //Set id
  MPI_Comm_size(MPI_COMM_WORLD, &p); //set p



  for (i=0;i<32;i++)
  {
    observations[i]/=1000.0;
  }

  //Set number of dipoles, either first argument or default value of 2.
  numberDipoles=2; 
  if (argc>=2)
  {
    numberDipoles=atoi(argv[1]);
  }

  //Set lambda based on entry, default of 40
  maxLambda=40;
  if (argc>=3)
  {
    maxLambda=atoi(argv[2]);
  }

  if (id==0)
  {
    printf("Dipoles:%d MaxLambda:%d\n",numberDipoles,maxLambda);
  }

  //Allocate lambda pieces to each processor, based on the size of maxLambda and the number of processors.
  lambda = BLOCK_SIZE(id,p,maxLambda);

  printf("Id:%d Lambda:%d\n",id,lambda);

  //Setup send and receive buffers for function evaluations and populations that resulted in those evaluation.
  sendCnts = malloc(p*sizeof(int));
  sdispls = malloc(p*sizeof(int));
  recvCnts = malloc(p*sizeof(int));
  rdispls = malloc(p*sizeof(int));
  sendCntsPop = malloc(p*sizeof(int));
  sdisplsPop = malloc(p*sizeof(int));
  recvCntsPop = malloc(p*sizeof(int));
  rdisplsPop = malloc(p*sizeof(int));

  for (i=0;i<p;i++)
  {

    sendCnts[i]=lambda;//Same for all others
    sdispls[i] = BLOCK_LOW(id,p,maxLambda);//Same for all others
    recvCnts[i] = BLOCK_SIZE(i,p,maxLambda);//Depends on which we receive from.
    rdispls[i] = BLOCK_LOW(i,p,maxLambda);

    sendCntsPop[i]=lambda*((numberDipoles*6+2));//Same for all others
    sdisplsPop[i] = BLOCK_LOW(id,p,maxLambda)*(numberDipoles*6+2);//Same for all others
    recvCntsPop[i] = BLOCK_SIZE(i,p,maxLambda)*(numberDipoles*6+2);//Depends on which we receive from.
    rdisplsPop[i] = BLOCK_LOW(i,p,maxLambda)*(numberDipoles*6+2);

  }

  for (i=0;i<p;i++)
  {

    printf("Id: %d recvCnts[%d]=%d\n",id,i,recvCnts[i]);
    printf("Id: %d rdispls[%d]=%d\n",id,i,rdispls[i]);

    printf("Id: %d recvCntsPop[%d]=%d\n",id,i,recvCntsPop[i]);
    printf("Id: %d rdisplsPop[%d]=%d\n",id,i,rdisplsPop[i]);

  }
  





  /* Initialize everything into the struct evo, 0 means default */
  //arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "initials.par"); 

//  printf("0\n");
  
  //The maxLambda parameter has been added so all of them will have enough space to store the results
  arFunvals = reinit(&evo, maxLambda, numberDipoles);

  //outputCMAES_t(evo,1);

//  printf("1\n");

  resetSignals(&evo, numberDipoles);  /* write header and initial values */

  //Reset the seed value based on processor (so they don't all come out the same!
  evo.sp.seed=evo.sp.seed*(id+1)/p;
  printf("proc: %d seed: %d\n",id,evo.sp.seed);


  //outputCMAES_t(evo,0);

//  printf("2\n");

//  printf("%s\n", cmaes_SayHello(&evo));
//  i=40;

//  for (i=32;i<40;i*=2)
//  { 

//    arFunvals = reinit(&evo, i);
    //outputCMAES_t(evo);


  evo.sp.lambda=lambda;
  canTerminate = (0==1);
  /* Iterate until stop criterion holds */
  while(!canTerminate)
    { 
      /* generate lambda new search points, sample population */
      pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

      /* Here you may resample each solution point pop[i] until it
	 becomes feasible, e.g. for box constraints (variable
	 boundaries). function is_feasible(...) needs to be
	 user-defined.  
	 Assumptions: the feasible domain is convex, the optimum is
	 not on (or very close to) the domain boundary, initialX is
	 feasible and initialStandardDeviations are sufficiently small
	 to prevent quasi-infinite looping.
      */
      /*for (i = 0; i < lambda; ++i) 
      {
          cmaes_ReSampleSingle(&evo, i); 
      }*/
      for (i = 0; i < lambda; ++i) 
      {
	   while (!is_feasible(evo.rgrgx[i],(int) cmaes_Get(&evo, "dim"))) 
	   {
             cmaes_ReSampleSingle(&evo, i); 
           }
      }

      for (i=0;i<lambda;i++)
      {
         for(j=0;j<(6*numberDipoles)+2;j++)
         {
	   evo.rgrgx[BLOCK_LOW(id,p,maxLambda)+i][j]=evo.rgrgx[i][j];
         }
      }
 
      /* evaluate the new search points using fitfun from above */ 
      for (i = BLOCK_LOW(id,p,maxLambda); i <= BLOCK_HIGH(id,p,maxLambda); ++i) {
	arFunvals[i] = fitfun(evo.rgrgx[i], (int) cmaes_Get(&evo, "dim"));
        //printf("ID:%d, arFunvals[%d]=%lf\n",id,i,arFunvals[i]);
      }

      



      //Now communicate the arFunvals around
      MPI_Alltoallv(arFunvals,sendCnts,sdispls,MPI_DOUBLE,arFunvals,recvCnts,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);


      //Now communicate the populations being looked at around
      MPI_Alltoallv(&evo.rgrgx[0][0],sendCntsPop,sdisplsPop,MPI_DOUBLE,&evo.rgrgx[0][0],recvCntsPop,rdisplsPop,MPI_DOUBLE,MPI_COMM_WORLD);


      /* update the search distribution used for cmaes_SampleDistribution() */
      cmaes_UpdateDistribution(&evo, arFunvals);  


      //Test for any that can terminate.
      canTerminate = cmaes_TestForTermination(&evo);
      if (canTerminate)
      {
	printf("id:%d can terminate for reason:%s\n",id,cmaes_TestForTermination(&evo));
      }
      MPI_Allreduce(&canTerminate,&canTerminateBuffer,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);//Get the max, if any are >0, then someone has terminated.
      canTerminate = canTerminateBuffer;//Reset so the loop will exit.

      /* read instructions for printing output or changing termination conditions */ 
//      cmaes_ReadSignals(&evo, "signals.par");   
//      fflush(stdout); /* useful in MinGW */
    }
//  printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */

//  cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */


  elapsed_time += MPI_Wtime();



  /* get best estimator for the optimum, xmean */
  xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
  bestValue = fitfun(xfinal, (int) cmaes_Get(&evo, "dim"));
  printf("Proccesor:%d has last mean of:%lf elapsedTime:%lf\n",id,bestValue,elapsed_time);
  for (i=0;i<6*numberDipoles;i++)
  {
    printf("(%d:%d:%lf)\n",id,i,xfinal[i]);
  }

//  cmaes_exit(&evo); /* release memory */ 
  /* do something with final solution and finally release memory */
  free(xfinal); 
  free(sendCnts);
  free(sdispls);
  free(recvCnts);
  free(rdispls);
  free(sendCntsPop);
  free(sdisplsPop);
  free(recvCntsPop);
  free(rdisplsPop);

  MPI_Finalize();

//}

  return 0;
}


double * reinit(cmaes_t * evo, int lambda, int numDipoles)
{
  int i, j, k, N;
  double dtest, trace;
  double *rgrgxdata;

/** Unknown if these are needed, most seem to be related to reading the init file (which we are removing) **/

  evo->sp.rgsformat = (char **) new_void(55, sizeof(char *));
  evo->sp.rgpadr = (void **) new_void(55, sizeof(void *)); 
  evo->sp.rgskeyar = (char **) new_void(11, sizeof(char *));
  evo->sp.rgp2adr = (double ***) new_void(11, sizeof(double **));
  evo->sp.weigkey = (char *)new_void(7, sizeof(char)); 

  /* All scalars:  */
  i = 0;
  evo->sp.rgsformat[i] = " N %d";        evo->sp.rgpadr[i++] = (void *) &evo->sp.N; 
  evo->sp.rgsformat[i] = " seed %d";    evo->sp.rgpadr[i++] = (void *) &evo->sp.seed;
  evo->sp.rgsformat[i] = " stopMaxFunEvals %lg"; evo->sp.rgpadr[i++] = (void *) &evo->sp.stopMaxFunEvals;
  evo->sp.rgsformat[i] = " stopMaxIter %lg"; evo->sp.rgpadr[i++] = (void *) &evo->sp.stopMaxIter;
  evo->sp.rgsformat[i] = " stopFitness %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.stStopFitness.val;
  evo->sp.rgsformat[i] = " stopTolFun %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.stopTolFun;
  evo->sp.rgsformat[i] = " stopTolFunHist %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.stopTolFunHist;
  evo->sp.rgsformat[i] = " stopTolX %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.stopTolX;
  evo->sp.rgsformat[i] = " stopTolUpXFactor %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.stopTolUpXFactor;
  evo->sp.rgsformat[i] = " lambda %d";      evo->sp.rgpadr[i++] = (void *) &evo->sp.lambda;
  evo->sp.rgsformat[i] = " mu %d";          evo->sp.rgpadr[i++] = (void *) &evo->sp.mu;
  evo->sp.rgsformat[i] = " weights %5s";    evo->sp.rgpadr[i++] = (void *) evo->sp.weigkey;
  evo->sp.rgsformat[i] = " fac*cs %lg";evo->sp.rgpadr[i++] = (void *) &evo->sp.cs;
  evo->sp.rgsformat[i] = " fac*damps %lg";   evo->sp.rgpadr[i++] = (void *) &evo->sp.damps;
  evo->sp.rgsformat[i] = " ccumcov %lg";    evo->sp.rgpadr[i++] = (void *) &evo->sp.ccumcov;
  evo->sp.rgsformat[i] = " mucov %lg";     evo->sp.rgpadr[i++] = (void *) &evo->sp.mucov;
  evo->sp.rgsformat[i] = " fac*ccov %lg";  evo->sp.rgpadr[i++]=(void *) &evo->sp.ccov;
  evo->sp.rgsformat[i] = " diagonalCovarianceMatrix %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.diagonalCov;
  evo->sp.rgsformat[i] = " updatecov %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.updateCmode.modulo;
  evo->sp.rgsformat[i] = " maxTimeFractionForEigendecompostion %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.updateCmode.maxtime;
  evo->sp.rgsformat[i] = " resume %59s";    evo->sp.rgpadr[i++] = (void *) evo->sp.resumefile;
  evo->sp.rgsformat[i] = " fac*maxFunEvals %lg";   evo->sp.rgpadr[i++] = (void *) &evo->sp.facmaxeval;
  evo->sp.rgsformat[i] = " fac*updatecov %lg"; evo->sp.rgpadr[i++]=(void *) &evo->sp.facupdateCmode;
  evo->sp.n1para = i; 
  evo->sp.n1outpara = i-2; /* disregard last parameters in WriteToFile() */

  /* arrays */
  i = 0;
  evo->sp.rgskeyar[i]  = " typicalX %d";   evo->sp.rgp2adr[i++] = &evo->sp.typicalX;
  evo->sp.rgskeyar[i]  = " initialX %d";   evo->sp.rgp2adr[i++] = &evo->sp.xstart;
  evo->sp.rgskeyar[i]  = " initialStandardDeviations %d"; evo->sp.rgp2adr[i++] = &evo->sp.rgInitialStds;
  evo->sp.rgskeyar[i]  = " diffMinChange %d"; evo->sp.rgp2adr[i++] = &evo->sp.rgDiffMinChange;
  evo->sp.n2para = i;  

/** End Unknown section **/

  //First reset the strategy parameters (sp)
  evo->sp.N=6*numDipoles;
  evo->sp.seed=0;
  evo->sp.xstart=NULL;
  evo->sp.typicalXcase=0;
  evo->sp.rgInitialStds = NULL;
  evo->sp.rgDiffMinChange = NULL; 
  evo->sp.stopMaxFunEvals = -1;
  evo->sp.stopMaxIter = -1;
  evo->sp.facmaxeval = 1; 
  evo->sp.stStopFitness.flg = 1;
  evo->sp.stStopFitness.val = 1e-5;
  evo->sp.stopTolFun = 1e-12; 
  evo->sp.stopTolFunHist = 1e-13; 
  evo->sp.stopTolX = 0; /* 1e-11*insigma would also be reasonable */ 
  evo->sp.stopTolUpXFactor = 1e3; 

  evo->sp.lambda=lambda;
  evo->sp.mu=-1;
  evo->sp.mucov = -1;
  evo->sp.weights = NULL;
  evo->sp.weigkey=malloc(5*sizeof(char));
  strcpy(evo->sp.weigkey, "log\0");

  evo->sp.cs = -1;
  evo->sp.ccumcov = -1;
  evo->sp.damps = -1;
  evo->sp.ccov = -1;

  evo->sp.diagonalCov = 0; /* default is 0, but this might change in future, see below */

  evo->sp.updateCmode.modulo = -1;  
  evo->sp.updateCmode.maxtime = -1;
  evo->sp.updateCmode.flgalways = 0;
  evo->sp.facupdateCmode = 1;

  strcpy(evo->sp.resumefile, "_no_");

//printf("reinit:0\n");

  readpara_SupplementDefaults(&(evo->sp));

//printf("reinit:1\n");
  //Now reset the state variables
  evo->version = "3.11.00.beta";
  evo->sp.seed = random_init( &evo->rand, (long unsigned int) evo->sp.seed);

  N = evo->sp.N; /* for convenience */


  //printf("reinit:2\n");

  evo->sp.rgInitialStds = new_double(N);
  evo->sp.xstart = new_double(N);
  for (k=0;k<numDipoles;k++)
  {
    evo->sp.rgInitialStds[k*6+0] = 3.000000;
    evo->sp.rgInitialStds[k*6+1] = 4.0;
    evo->sp.rgInitialStds[k*6+2] = 5.0;
    evo->sp.rgInitialStds[k*6+3] = 90.0;
    evo->sp.rgInitialStds[k*6+4] = 90.0;
    evo->sp.rgInitialStds[k*6+5] = 3.75;


    evo->sp.xstart[k*6+0] = 6.000000;
    evo->sp.xstart[k*6+1] = 8.0;
    evo->sp.xstart[k*6+2] = 15.0;
    evo->sp.xstart[k*6+3] = 90.0;
    evo->sp.xstart[k*6+4] = 90.0;
    evo->sp.xstart[k*6+5] = 0.5;
  }
  /* initialization  */
  for (i = 0, trace = 0.; i < N; ++i)
    trace += evo->sp.rgInitialStds[i]*evo->sp.rgInitialStds[i];
  evo->sigma = sqrt(trace/N); /* evo->sp.mueff/(0.2*evo->sp.mueff+sqrt(N)) * sqrt(trace/N); */

//printf("reinit:3\n");

  evo->chiN = sqrt((double) N) * (1. - 1./(4.*N) + 1./(21.*N*N));
  evo->flgEigensysIsUptodate = 1;
  evo->flgCheckEigen = 0; 
  evo->genOfEigensysUpdate = 0;
  timings_init(&evo->eigenTimings);
  evo->flgIniphase = 0; /* do not use iniphase, hsig does the job now */
  evo->flgresumedone = 0;
  evo->flgStop = 0;


  for (dtest = 1.; dtest && dtest < 1.1 * dtest; dtest *= 2.) 
    if (dtest == dtest + 1.)
      break;
  evo->dMaxSignifKond = dtest / 1000.; /* not sure whether this is really safe, 100 does not work well enough */

  evo->gen = 0;
  evo->countevals = 0;
  evo->state = 0;
  evo->dLastMinEWgroesserNull = 1.0;
  evo->printtime = evo->writetime = evo->firstwritetime = evo->firstprinttime = 0; 

  evo->rgpc = new_double(N);
  evo->rgps = new_double(N);
  evo->rgdTmp = new_double(N+1);
  evo->rgBDz = new_double(N);
  evo->rgxmean = new_double(N+2); evo->rgxmean[0] = N; ++evo->rgxmean;
  evo->rgxold = new_double(N+2); evo->rgxold[0] = N; ++evo->rgxold; 
  evo->rgxbestever = new_double(N+3); evo->rgxbestever[0] = N; ++evo->rgxbestever; 
  evo->rgout = new_double(N+2); evo->rgout[0] = N; ++evo->rgout;
  evo->rgD = new_double(N);
  evo->C = (double**)new_void(N, sizeof(double*));
  evo->B = (double**)new_void(N, sizeof(double*));
  evo->publicFitness = new_double(evo->sp.lambda); 
  evo->rgFuncValue = new_double(evo->sp.lambda+1); 
  evo->rgFuncValue[0]=evo->sp.lambda; ++evo->rgFuncValue;
  evo->arFuncValueHist = new_double(10+(int)ceil(3.*10.*N/evo->sp.lambda)+1);
  evo->arFuncValueHist[0] = (double)(10+(int)ceil(3.*10.*N/evo->sp.lambda));
  evo->arFuncValueHist++; 

//printf("reinit:4\n");

  for (i = 0; i < N; ++i) {
      evo->C[i] = new_double(i+1);
      evo->B[i] = new_double(N);
    }
  evo->index = (int *) new_void(evo->sp.lambda, sizeof(int));
  for (i = 0; i < evo->sp.lambda; ++i) 
    evo->index[i] = i; /* should not be necessary */
  evo->rgrgx = (double **)new_void(evo->sp.lambda, sizeof(double*));
  
  rgrgxdata = new_double(evo->sp.lambda*(N+2));
  for (i = 0; i < evo->sp.lambda; ++i) {
    evo->rgrgx[i] = &rgrgxdata[i*(N+2)];
    evo->rgrgx[i][0] = N; 
    evo->rgrgx[i]++;
  }
//  printf("Arggg:%d\n",evo->sp.lambda*(N+2));

//printf("reinit:5\n");
  /* Initialize newed space  */

  for (i = 0; i < N; ++i)
    for (j = 0; j < i; ++j)
       evo->C[i][j] = evo->B[i][j] = evo->B[j][i] = 0.;
        
  for (i = 0; i < N; ++i)
    {
      evo->B[i][i] = 1.;
      evo->C[i][i] = evo->rgD[i] = evo->sp.rgInitialStds[i] * sqrt(N / trace);
      evo->C[i][i] *= evo->C[i][i];
      evo->rgpc[i] = evo->rgps[i] = 0.;
    }

//printf("reinit:6\n");

  evo->minEW = rgdouMin(evo->rgD, N); evo->minEW = evo->minEW * evo->minEW;
  evo->maxEW = rgdouMax(evo->rgD, N); evo->maxEW = evo->maxEW * evo->maxEW; 

  evo->maxdiagC=evo->C[0][0]; for(i=1;i<N;++i) if(evo->maxdiagC<evo->C[i][i]) evo->maxdiagC=evo->C[i][i];
  evo->mindiagC=evo->C[0][0]; for(i=1;i<N;++i) if(evo->mindiagC>evo->C[i][i]) evo->mindiagC=evo->C[i][i];

  /* set xmean */
  for (i = 0; i < N; ++i)
  {
    evo->rgxmean[i] = evo->rgxold[i] = evo->sp.xstart[i]; 

  }



//printf("reinit:7\n");
  return evo->publicFitness;
}

void resetSignals(cmaes_t *evo, int numDipoles)
{
  int k;
  //Reset typical x -- all other mandatory settings are overwritten in reinit()

  evo->sp.typicalX = new_double(evo->sp.N);

  for (k=0;k<numDipoles;k++)
  {
    evo->sp.typicalX[k*6+0] = 6.000000;
    evo->sp.typicalX[k*6+1] = 8.0;
    evo->sp.typicalX[k*6+2] = 15.0;
    evo->sp.typicalX[k*6+3] = 90.0;
    evo->sp.typicalX[k*6+4] = 90.0;
    evo->sp.typicalX[k*6+5] = 0.5;
  }


  //Reset optional settings:

  evo->sp.stopMaxFunEvals=1e299;
  evo->sp.facmaxeval=1.0;
  evo->sp.stopMaxIter=1e299;
  evo->sp.stopTolFun=1e-12;
  evo->sp.stopTolFunHist=1e-13;
  evo->sp.stopTolX=1e-11;
  evo->sp.stopTolUpXFactor=1e3;

  evo->sp.updateCmode.maxtime=1;
  
}

//Doug Hakkarinen - added for bounds checking.
int is_feasible(double * dipole, int N)
{
  int i;
  for (i=0;i<N;i++)
  {
    if (dipole[i]<MINVALUES[i%6]||dipole[i]>MAXVALUES[i%6])
    {
      return (1==0);//Return false
    }
  }
  return (0==0);//Return true
}

