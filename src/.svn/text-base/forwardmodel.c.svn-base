#include "forwardmodel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define PI 3.1415926
#define DIPOLE_SOURCE_SEPARATION 0.01
#define DEBUG_FORWARD_MODEL 0


/** Features are parts of the dipole
  * The featureValues array should be structured as:
    Feature i*6+0: ith x value
    Feature i*6+1: ith y value
    Feature i*6+2: ith z value
    Feature i*6+3: ith theta value
    Feature i*6+4: ith phi value
    Feature i*6+5: ith mag value
   //This function will return NULL if numberFeatures%6 != 0
*/
double * forwardModel(int numberFeatures,const double * featureValues, double rho, double xmax, double ymax, double zmax)
{
  double * result;
  Dipole * dipoles;
  int i; 
  int j; 
  int observationCount;
  double observationZ[32];
  double observationX[]={2,2,2,2,2,2,2,4,4,4,4,4,4,4,6,6,6,6,6,6,6,8,8,8,8,8,8,8,10,10,10,10};
  double observationY[]={2,4,6,8,10,12,14,2,4,6,8,10,12,14,2,4,6,8,10,12,14,2,4,6,8,10,12,14,2,4,6,8};


  
  if (numberFeatures%6!=0)
   return NULL;

  observationCount=32;

  
  for (i=0;i<observationCount;i++)
  {
    observationZ[i]=zmax;
  }

  
  result = (double*)malloc(observationCount*sizeof(double));
  for (i=0;i<observationCount;i++)
  {
    result[i]=0.0;
  }


  dipoles = (Dipole*)malloc(numberFeatures/6*sizeof(Dipole));

  //Loop through the dipoles
  for (i=0;i<numberFeatures/6;i++)
  {
    Coord3d iPlusPosition; //Location of positive pole
    Coord3d iMinusPosition; //Location of negative pole
    double theta; //Theta in radians
    double phi; //Phi in radians
    double vCoeff; //Voltage coefficient
    double plusZDepth; //Depth of positive pole
    double minusZDepth; //Depth of negative pole

    dipoles[i].pos.x=featureValues[i*6];
    dipoles[i].pos.y=featureValues[i*6+1];
    dipoles[i].pos.z=featureValues[i*6+2];
    dipoles[i].theta=featureValues[i*6+3];
    dipoles[i].phi=featureValues[i*6+4];
    dipoles[i].mag=featureValues[i*6+5];

    theta = dipoles[i].theta*PI/180.0;
    phi = dipoles[i].phi*PI/180.0;

    //Now get the positive and negative dipole centers
    iPlusPosition = dipoles[i].pos;
    iMinusPosition = dipoles[i].pos;
    iPlusPosition.x+= DIPOLE_SOURCE_SEPARATION/2.0*sin(theta)*cos(phi);
    iMinusPosition.x-= DIPOLE_SOURCE_SEPARATION/2.0*sin(theta)*cos(phi);
    iPlusPosition.y+= DIPOLE_SOURCE_SEPARATION/2.0*sin(theta)*sin(phi);
    iMinusPosition.y-= DIPOLE_SOURCE_SEPARATION/2.0*sin(theta)*sin(phi);
    iPlusPosition.z+= DIPOLE_SOURCE_SEPARATION/2.0*cos(theta);
    iMinusPosition.z-= DIPOLE_SOURCE_SEPARATION/2.0*cos(theta);

    //Now set up the voltage coefficient
    vCoeff = dipoles[i].mag*rho/(4*PI);

    //Set up the depths
    plusZDepth=zmax-iPlusPosition.z;
    minusZDepth=zmax-iMinusPosition.z;
   
    if (DEBUG_FORWARD_MODEL)
    {
       printf("Dipole:%d;Plus:(%lf,%lf,%lf);Minus:(%lf,%lf,%lf);vCoeff:%lf;a_plus:%lf;a_minus:%lf\n",i,iPlusPosition.x,iPlusPosition.y,iPlusPosition.z,iMinusPosition.x,iMinusPosition.y,iMinusPosition.z,vCoeff,plusZDepth,minusZDepth);
    }

    //Loop through the observation points
    for (j=0;j<observationCount;j++)
    {
      double voltage;
      double voltage1,voltage2,voltage3;
      double radius_x;
      double radius_y;
      double x_c;
      double y_c;
      //Generate the Voltage for the observation point for  the positive pole
      voltage=0.0;
      radius_x=distance(iPlusPosition.x,observationX[j]);
      radius_y=distance(iPlusPosition.y,observationY[j]);
      x_c = xmax-iPlusPosition.x;
      y_c = ymax-iPlusPosition.y;
      voltage1=vz_reflection_sum(1.0, radius_x,radius_y,zmax,plusZDepth, 5);

      if (DEBUG_FORWARD_MODEL)
      {
        printf("Observation:%d;r_x:%lf;r_y:%lf;x_c:%lf;y_c:%lf\n",j,radius_x,radius_y,x_c,y_c);
      }

      voltage2=vxy_reflection_sum(1.0, radius_y,radius_x,ymax,plusZDepth,y_c,1);
      voltage3=vxy_reflection_sum(1.0, radius_x,radius_y,xmax,plusZDepth,x_c,1);
      voltage+=voltage1+voltage2+voltage3;
      voltage*=vCoeff;
      if (DEBUG_FORWARD_MODEL)
      {
        printf("Observation:%d;r_x:%lf;r_y:%lf;x_c:%lf;y_c:%lf;v1:%lf;v2:%lf;v3:%lf;vplus:%lf",j,radius_x,radius_y,x_c,y_c,voltage1,voltage2,voltage3,voltage);
      }

      //Add that voltage to the previous total for that observation point
      result[j]+=voltage;


      //Generate the Voltage for the observation point for  the negative pole
      voltage=0.0;
      radius_x=distance(iMinusPosition.x,observationX[j]);
      radius_y=distance(iMinusPosition.y,observationY[j]);
      x_c = xmax-iMinusPosition.x;
      y_c = ymax-iMinusPosition.y;
      voltage1=vz_reflection_sum(1.0,radius_x,radius_y,zmax,minusZDepth, 5);
      voltage2=vxy_reflection_sum(1.0, radius_y,radius_x,ymax,minusZDepth,y_c,1);
      voltage3=vxy_reflection_sum(1.0, radius_x,radius_y,xmax,minusZDepth,x_c,1);
      voltage+=voltage1+voltage2+voltage3;
      voltage*=-vCoeff;
      if (DEBUG_FORWARD_MODEL)
      {
        printf(";v1-:%lf;v2-:%lf;v3-:%lf;vminus:%lf\n",voltage1,voltage2,voltage3,voltage);
      }

      //Add that voltage to the previous total for that observation point
      result[j]+=voltage;


    }
  }

  return result;
}

double vz_reflection_sum(double kall, double r_x, double r_y, double z, double a, int numberTerms)
{
  int i;
  double result;
  result=0.0;
  for (i=1;i<=numberTerms;i++)
  {
    double numerator;
    numerator=0.5*pow(kall,i);
    result+=numerator/sqrt( ((2*i-2)*z+a)*((2*i-2)*z+a) + r_x*r_x + r_y*r_y); 
    result+=numerator/sqrt( ((2*i-1)*z+z-a)*((2*i-1)*z+z-a) + r_x*r_x + r_y*r_y );
    if (DEBUG_FORWARD_MODEL)
    {
      //printf("vz_reflection_sum, iteration:%d;result:%lf\n",i,result);
    }
  }
  return result;
}

double vxy_reflection_sum(double kall, double r_x, double r_y, double z, double a, double c, int numberTerms)
{
  int i;
  double result;
  result=0.0;

  for (i=1;i<=numberTerms;i++)
  {
     double numerator;
     numerator=0.5*pow(kall,i);
     result+=numerator/sqrt( (2*i*z - 2*c + r_x)*(2*i*z - 2*c + r_x) + r_y*r_y + a*a);
     result+=numerator/sqrt( ((2*i-2)*z + 2*c -r_x)*((2*i-2)*z + 2*c -r_x) + r_y*r_y + a*a);
     result+=numerator/sqrt( (2*i*z+r_x)*(2*i*z+r_x) +r_y*r_y + a*a);
     result+=numerator/sqrt( (2*i*z-r_x)*(2*i*z-r_x) +r_y*r_y + a*a);
  }

  result += 0.5*pow(kall,i)/sqrt(a*a+r_x*r_x+r_y*r_y);

  return result;
}
/**
 *  Distance function used for forward model.
*/
double distance(double x0, double x1)
{
  return sqrt((x0-x1)*(x0-x1));
}

const int MAX_LINE_SIZE=1024;
const int MAX_OBSERVATIONS=256;

#define DEBUG_CMAES 1

/**
  The objective function for comparing observations with model outputs.
*/
double objectiveFunction(int numberInputs, double * inputs, double * observations)
{
  double result;
  int i;
  result = 0.0;
  for (i=0;i<numberInputs;i++)
  {
	  result+=(inputs[i]-observations[i])*(inputs[i]-observations[i]);
  }
  return sqrt(result/numberInputs);//Matches Allan's objective function.
}

Inputs getObservations(const char * csvfile, int colWanted)
{
    Inputs result;
    FILE* file;
    char buffer[MAX_LINE_SIZE];
    double resultTemp[MAX_OBSERVATIONS];
    int row=0;
    int column=0;
    int j;
    int unused;
    char * lastFieldStart;

    file = fopen(csvfile,"r");
    if (file == NULL)
    {
        result.numberInputs=0;
        result.inputs=0;
    	fprintf(stderr,"Could not open input file. \n");
	return result;
    }
    if (DEBUG_CMAES)
    {
      printf("File opened\n");
    }

    //Get rid of header row.
    lastFieldStart  = fgets(buffer,MAX_LINE_SIZE,file);
    if (DEBUG_CMAES)
    {
      printf("Header removed:%s\n",lastFieldStart);
    }

    //Read one line at a time into buffer.
    while (fgets(buffer,MAX_LINE_SIZE,file)!=NULL)
    {
    	char * lastStart;
    	int i;

        if (DEBUG_CMAES)
        {
           printf("Line read:%s\n",buffer);
        }
    	i=0;
    	column = 0;
        lastStart=&buffer[0];
    	//Read through the buffer, finding each comma until we find the field we are interested in.
        while(buffer[i]!='\n'&&buffer[i]!='\r'&&buffer[i]!=EOF)
        {
          
          if (buffer[i]==' '||buffer[i]==',')
          {
        	  buffer[i]='\0';
        	  if (column==colWanted)
        	  {
        		unused=sscanf(lastStart,"%lf",&resultTemp[row]);
			break;
        	  }
		  
                  lastStart=&buffer[i+1];
        	  column++;
          }
          i++;
        }
        row++;
    }

	fclose(file);
	result.inputs = malloc(row*sizeof(double));
	for (j=0;j<row;j++)
	{
	    result.inputs[j]=resultTemp[j];
	}
	result.numberInputs=row;
	return result;
}



