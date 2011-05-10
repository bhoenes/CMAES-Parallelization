/*
 *
 *  Created on: Sep 21, 2010
 *      Author: dhakkarinen
 */

#ifndef FORWARDMODEL_H
#define FORWARDMODEL_H

#include "dipole.h"
#include "feature.h"
#include "inputs.h"

/*
 * Returns the value of evaluating the objective function for a given set of inputs (and observations)
 *
 */
double objectiveFunction(int numberInputs, double * inputs, double * observations);

/*
 * Returns the Input for a csv file for a specified column, ignoring the first line as the header.
 *
 */
Inputs getObservations(const char * csvfile, int colWanted);

/**
  Runs the forward model on the featureValues to generate the outputs that can then be used as inputs to the objective function.
*/
double* forwardModel(int numberFeatures,const double * featureValues, double rho, double xmax, double ymax, double zmax);


/**
 *  Distance function used for forward model.
*/
double distance(double x0, double x1);

/**
* Calculates the voltage from the z axis.
*/
double vz_reflection_sum(double kall, double r_x, double r_y, double z, double a, int numberTerms);

/**
* Calculates the voltage from the z axis.
*/
double vxy_reflection_sum(double kall, double r_x, double r_y, double z, double a, double c, int numberTerms);


#endif

