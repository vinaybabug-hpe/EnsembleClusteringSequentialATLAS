/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * cluster_util_partition2cam.c
 *
 * Code generation for function 'cluster_util_partition2cam'
 *
 */


/* Include files */
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cstdio>



/* Function Definitions */
void cluster_util_partition2cam(int* I, int* A, int* N, int nSamples, int nPartitions)
{


  /* CLUSTER_UTIL_PARTITION2CAM Generates a co-association matrix from a partition  */
  /*  */
  /*  Syntax:  A = cluster_util_partition2cam(I)    */
  /*  */
  /*    INPUT ARGUMENTS */
  /*         I : [n,1] column of partition indices */
  /*  */
  /*    OUTPUT ARGUMENTS */
  /*    	    A :  Co-association matrix  [n,n] where n = num samples */
  /*                 Tracks which pairs of samples are in the same cluster. */
  /*                  */
  /*            N : Co-occurrence matrix [n,n]    */
  /*                Tracks which pairs are both present in the sample. If I */
  /*                does not contain zeros or NaNs, then all pairs are present */
  /*                and N is just the identity matrix.  If there are zeros or */
  /*                NaNs present (as a result of bootstrapping), then not all */
  /*                pairs of samples are present, N will contain zeros. */
  /*                  */
  /*  Author:         Lee I Newman */
  /*  Affiliation:    University of Michigan, Depts. Psychology, EECS */
  /*  email:          leenewm@umich.edu */
  /*  Website:        http://www.leenewman.org */
  /*  Created:        May 2009 */
  /*  Revised by/on:  person, date */
  /*  */
  /* ------------- BEGIN CODE -------------- */
  /* % Process input arguments and initialize variables */
  /*  ///////////////////////////////////////////////////////// */
  /*  size variables */
  /*  identify samples that have 0 or NaN in the partition */
  /* % Compute the co-association matrix (similarity format) */
  /*  initialize output variables with 1's on diagonal for non-zero/non-NaN */
  /*  elements and zeros everywhere else */

  memset(A, 0, nSamples * nSamples * sizeof(int));
  memset(N, 0, nSamples * nSamples * sizeof(int));



//  /*  co-occurrence matrix */
//  /*  loop over all pairs of points */
  for (int i = 0; i < nSamples; i++) {

    for (int j = i+1; j < nSamples; j++) {

    	if(I[i] == I[j] && I[i] != 0){
    		 A[j + nSamples * i] = 1;
    		 A[i + nSamples * j] = 1;
    	}

    	if(I[i]!= 0 && I[j] != 0){
    		N[j + nSamples * i] = 1;
    		N[i + nSamples * j] = 1;
    	}

    }

    A[i + nSamples * i] = 1.0;
   	N[i + nSamples * i] = 1.0;

    /*  end loop over j */
  }

  /*  end loop over i */
  /* ------------- END OF CODE -------------- */
}

/* End of code generation (cluster_util_partition2cam.c) */
