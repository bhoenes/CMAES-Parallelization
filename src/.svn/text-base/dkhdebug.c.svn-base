//This file provides implementations for several debug functions created by Doug Hakkarinen

#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.h"
#include "forwardmodel.h"
#include "cmaes.h"
#include <math.h>
#include <string.h>
#include "dkhdebug.h"

void outputCMAES_t(cmaes_t t, int subs)
{
  int n;
  n = t.sp.N;
  printf("Output for cmaes_t object:\n");
  printf("%s\n",t.version);
  if (subs==0)
  {
    outputReadPara(t.sp);
  }
  printf("%lf",t.sigma); 

  outputArray(t.rgxmean,n);  
  outputArray(t.rgxbestever,n); 

  printf("%lf\n",t.chiN);

  printf("%lf\n",t.gen);
  printf("%lf\n",t.countevals);
  printf("%lf\n",t.state);

  printf("%lf\n",t.maxdiagC);
  printf("%lf\n",t.mindiagC);
  printf("%lf\n",t.maxEW);
  printf("%lf\n",t.minEW);

  printf("End output for cmaes_t object.\n");

}

void outputArray(double * data, int n)
{
  int i;
  for (i=0;i<n;i++)
  {
    printf("%lf ",data[i]);
  }
  printf("\n");
}

void outputStringArray(char ** data, int n)
{
  int i;
  for (i=0;i<n;i++)
  {
    printf("%s ",data[i]);
  }
  printf("\n");
}

void outputReadPara(readpara_t sp)
{
   printf("Output of readpara_t\n");
   printf("Input Parameter\n");
   printf("N:%d\n",sp.N);
   printf("seed:%d\n",sp.seed);
   outputArray(sp.xstart,sp.N);
   outputArray(sp.typicalX,sp.N);
   printf("typicalXcase:%d\n",sp.typicalXcase);
   outputArray(sp.rgInitialStds,sp.N);
   //outputArray(sp.rgDiffMinChange,sp.N);
 
   printf("Termination Parameters\n");
   printf("stopMaxFunEvals:%lf\n",sp.stopMaxFunEvals);
   printf("facmaxeval:%lf\n",sp.facmaxeval);
   printf("stopMaxIter:%lf\n",sp.stopMaxIter);
   printf("stStopFitness.flg .val:%d %lf\n",sp.stStopFitness.flg, sp.stStopFitness.val);
   printf("stopTolFun:%lf\n",sp.stopTolFun);
   printf("stopTolFunHist:%lf\n",sp.stopTolFunHist);
   printf("stopTolX:%lf\n",sp.stopTolX);
   printf("stopTolUpXFactor:%lf\n",sp.stopTolUpXFactor);

   printf("Internal Evolutionary Strategy Parameters\n");
   printf("lambda:%d\n",sp.lambda);
   printf("mu:%d\n",sp.mu);
   printf("mucov:%lf\n",sp.mucov);
   printf("mueff:%lf\n",sp.mueff);
   outputArray(sp.weights,sp.mu);
   printf("damps:%lf\n",sp.damps);
   printf("cs:%lf\n",sp.cs);
   printf("ccumcov:%lf\n",sp.ccumcov);
   printf("ccov:%lf\n",sp.ccov);
   printf("diagonalCov:%lf\n",sp.diagonalCov);
   printf("updateCmode.flgalways .modulo .maxtime:%d %lf %lf\n",sp.updateCmode.flgalways,sp.updateCmode.modulo,sp.updateCmode.maxtime);  
   printf("facupdateCmode:%lf\n",sp.facupdateCmode);

   printf("weigkey:%s\n", sp.weigkey);
   printf("resumefile:%s\n",sp.resumefile);
   outputStringArray(sp.rgsformat,sp.n1outpara);
   printf("n1para n1outpara n2para:%d %d %d\n",sp.n1para,sp.n1outpara,sp.n2para);
   printf("End readpara_t\n");
}


void my_cmaes_exit(cmaes_t *t)
{
  //int exitpoint;
  int i, N = t->sp.N;
  
  t->state = -1; /* not really useful at the moment */
  //exitpoint=0;
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;

  free( t->rgpc);
  free( t->rgps);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( t->rgdTmp);
  free( t->rgBDz);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( --t->rgxmean);
  free( --t->rgxold); 
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( --t->rgxbestever); 
  free( --t->rgout); 
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( t->rgD);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  for (i = 0; i < N; ++i) {
    free( t->C[i]);
    free( t->B[i]);
  }
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  for (i = 0; i < t->sp.lambda; ++i) 
    free( --t->rgrgx[i]);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( t->rgrgx); 
  free( t->C);
  free( t->B);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( t->index);
  free( t->publicFitness);
  //printf("cmaes_exit:%d\n",exitpoint);exitpoint++;
  free( --t->rgFuncValue);
  free( --t->arFuncValueHist);
  random_exit (&t->rand);
  my_readpara_exit (&t->sp); 
} /* cmaes_exit() */

void my_readpara_exit(readpara_t *t)
{
  //int exitpoint; 
  //exitpoint=0;
  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;

  if (t->xstart != NULL) /* not really necessary */
    free( t->xstart);

  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;
  if (t->typicalX != NULL)
    free( t->typicalX);
  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;
  if (t->rgInitialStds != NULL)
    free( t->rgInitialStds);
  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;
  if (t->rgDiffMinChange != NULL)
    free( t->rgDiffMinChange);
  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;
  if (t->weights != NULL)
    free( t->weights);
  //printf("readpara_exit:%d\n",exitpoint);exitpoint++;
  free(t->rgsformat);

  free(t->rgpadr);

  free(t->rgskeyar);

  free(t->rgp2adr);

  free(t->weigkey);
}

