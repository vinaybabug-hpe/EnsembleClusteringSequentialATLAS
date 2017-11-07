/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * cluster_util_ssw.c
 *
 * Code generation for function 'cluster_util_ssw'
 *
 */

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_set>
#include <set>
#include<cassert>

#include "common/indices_count.h"
#include "common/wrapperFuncs.h"

/* Function Definitions */
void cluster_util_ssw(float *X, int nSamples, int nFeatures, int  *I, double *ssw_data, int nClusters)
{


  /* CLUSTER_UTIL_SSW Compute the total SSE within clusters */
  /*  */
  /*    syntax:  [SSW,ssw] = cluster_util_ssw(X,I,metric) */
  /*  */
  /*    SSW is the total sum of squares within for all clusters (ergo, Total) */
  /*    ssw is an array with the sum squared error within, for each cluster */
  /*  */
  /*    X       : [n,d] data matrix */
  /*    I       : [n,1] column of cluster indices */
  /*   metric   : 'mean' (default) or 'median' */
  /*  set default for metric if not provided */
  /*  compute size variables */

  memset(ssw_data, 0, nClusters * sizeof(double));

  std::unordered_set<int> s(I, I + nSamples);
  std::set<int> nValues(s.begin(), s.end());

  int *valueList = (int*) malloc(nClusters * sizeof(int));
  indices_count(I, nSamples, valueList, nClusters);

  float *M_data = (float*) malloc(nFeatures * sizeof(float));
  assert(M_data);
  float *sw = (float*) malloc(nFeatures * nFeatures * sizeof(float));
  assert(sw);
//  printf("\n");
  for(int iValue = 0; iValue < nClusters; iValue++){
	  int label = *std::next(nValues.begin(), iValue);
//	  printf("%d ", valueList[iValue]);
	  /**
	   * number of points in the cluster
	   */
	  int n = valueList[iValue];
	 /**
	  * extract subset of data points in cluster k
	  */
	  float *x_data = (float*)malloc(n * nFeatures * sizeof(float));
	  int n_cnt = 0;
	  for (int j = 0; j < nSamples; j++) {
		  if(I[j] == label){
			  memcpy(&x_data[n_cnt * nFeatures], &X[j * nFeatures], nFeatures * sizeof(float));
			  n_cnt++;
		  }
	  }

	  /**
	   * cluster center
	   */
	  memset(M_data, 0, nFeatures * sizeof(float));
	  for(int col =0; col < nFeatures; col++){
		  float mean = 0;
		  int cntMean = 0;

		  for(int row =0; row < n; row++){

				  mean += X[row * nFeatures + col];


		  }
		  if(cntMean != 0)
			  M_data[col] = mean/cntMean;
		  else
			  M_data[col] = 0;
	  }

	  /**
	   * centered the data
	   */
	  for(int row =0; row < n; row++){
		  for(int col =0; col < nFeatures; col++){
			  x_data[row *nFeatures + col] = x_data[row *nFeatures + col] - M_data[col];
		  }
	  }

	  memset(sw, 0, nFeatures * nFeatures * sizeof(float));

	  float *transpose = (float *) malloc(n * nFeatures * sizeof(float));
	  assert(transpose);

	  // Compute transpose of centered data
	  matrixTranspose(transpose, x_data, n, nFeatures);

	  /**
	   * scatter matrix for the cluster
	   */
	  matrixMulCPU(sw, transpose, x_data, nFeatures, n, nFeatures);
	  /**
	   * sum sq error within cluster k
	   * essentially it is the sum of the diagonal elements of the matrix sw
	   */
	  float trace = 0;
	  for (int k = 0; k < nFeatures; k++) {
	    trace += sw[k + nFeatures * k];
	  }

	  ssw_data[iValue] = trace;

	  free(transpose);
	  free(x_data);
  }


//  printf("\n");
  free(sw);
  free(M_data);

  /*  total within-cluster sum squares  */
}

/* End of code generation (cluster_util_ssw.c) */
