/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * cluster_data_outliers.c
 *
 * Code generation for function 'cluster_data_outliers'
 *
 */

/* Include files */


//#include "sort1.h"
//#include "zscore.h"
//#include "rand.h"
#include <stdbool.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <regex.h>
#include <string>
#include <cstring>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <numeric>
#include <stdlib.h>
#include "mat.h"

#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "common/distCalcMthds.h"
#include "common/zscore.h"



using namespace std;
/* Function Definitions */
float* cluster_data_outliers(int nrows, int ncols, float* X, char* distFun, float cutoff, int *newNElements)
{



	char method, dist;
	std::string distFunStr(distFun);

  /*  function [outlierMask,outlierIdxs,outlierZ,outlierDS] = cluster_data_outliers(X,distFun,cutoff) */
  /* CLUSTER_DATA_OUTLIERS - Finds outliers in [n,p] data X where p are features */
  /*  */
  /*  syntax: [outlierMask,outlierIdxs,outlierZ,outlierDS] = ... */
  /*            cluster_data_outliers(X,distFun,cutoff) */
  /*  */
  /*   X        [n,p] */
  /*   distFun  'mahal' or 'euc' */
  /*   cutoff   #stdevs  */
  /*  */


	if(cutoff <= 3){
		cutoff = 3;
	}

	if (distFunStr.compare("mahal") == 0 || distFunStr.compare("mahalanobis") == 0){

		dist = 'm';
	}
	else if (distFunStr.compare("euc") == 0 || distFunStr.compare("euclidean") == 0){

		dist = 'e';
	}

	// Find distance(X, X)

	//	int i/*, j, k*/;
	int transpose = 0;
	const int nelements = (transpose == 0) ? nrows : ncols;
	int** cmask;
	float** distmatrix = NULL;
	float* weight;
	float** data;

	weight = (float*)malloc(ncols * sizeof(float*));

	malloc2D(cmask, nrows, ncols, int);

	malloc2D(data, nrows, ncols, float);
	memcpy(data[0], X, nrows * ncols * sizeof(float));

	for (int row = 0; row < nrows; row++){
		for (int col = 0; col < ncols; col++){
			cmask[row][col] = data[row][col] == 0? 0 : 1;
		}
	}

	assert(weight != NULL);

	for (int i = 0; i < ncols; i++)
		weight[i] = 1.0;

	distmatrix =
	      distancematrix(nrows, ncols, data, cmask, weight, dist, transpose);
	if (!distmatrix) return NULL; /* Insufficient memory */


	float *Z = (float*) malloc(nelements * sizeof(float));


		for (int i = 0; i < nelements; i++){
			for (int j = 0; j < i; j++){
//				distmatrix_h [i * nrows + j] = distmatrix_h [j * nrows + i] =  distmatrix_m[TRI_COUNT(i)+j];
				Z[i] += distmatrix[i][j];
				Z[j] += distmatrix[i][j];
//				printf("%.4f ", distmatrix_m[TRI_COUNT(i)+j]);
			}
//			printf("\n");
			Z[i] /= nelements;
		}

	 zscore(nelements, Z);

	 int *outlierIdx = (int*) malloc(nelements * sizeof(int));


	 // compute  outlierIdx[k] = Z[k] >= cutoff?1:0;
	 for(int count=0; count<nelements; count++){
		 outlierIdx[count] = abs(Z[count]) > cutoff?1:0;

	 }


	 int numOutliers = 0;
	 for(int count=0; count<nelements; count++){
	 		 numOutliers += outlierIdx[count];
	 }

	 *newNElements = nelements-numOutliers;
	 float *newX;
	 newX = (float*) malloc((*newNElements)*ncols * sizeof(float));

	 if(numOutliers > 0){
		 for (int j = 0, j1=0 ; j < nelements; j++) {
			 if(outlierIdx[j] == 0){
				 memcpy(&newX[j1 * ncols], &X[j * ncols], ncols * sizeof(float));
				 j1++;
			 }
		 }
	 }
	 else{
		 memcpy(newX, X, (*newNElements)*ncols * sizeof(float));

	 }

	free(weight);
	free(cmask);
	for (int i = 1; i < nelements; i++) free(distmatrix[i]);
	free (distmatrix);
	free(data[0]);
	free(data);


	 //TODO: Debug code remove later!
//	 printf("\n-------------------------------------------------------------------------------------- \n\n");
//	 for (int m = 0; m < *newNElements; m++) {
//		 for (int n = 0; n < ncols; n++) {
//			 printf("%.3f ", newX[n + m * ncols]);
//		 }
//		 printf(" \n");
//	 }

	 return newX;

}

/* End of code generation (cluster_data_outliers.c) */
