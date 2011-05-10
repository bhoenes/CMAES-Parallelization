/*
 * feature.h
 *
 *  Created on: Sep 21, 2010
 *      Author: dhakkarinen
 */

#ifndef FEATURE_H_
#define FEATURE_H_

struct feature_struct
{
	double min;
	double max;
	double mean;
	double sigma;
};

typedef struct feature_struct Feature;

#endif /* FEATURE_H_ */
