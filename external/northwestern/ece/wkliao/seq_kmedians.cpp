/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_kmeans.c  (sequential version)                        */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "common/distCalcMthds.h"
#include "common/wrapper.h"
#include "common/wrapperFuncs.h"

/*----< euclid_dist_2() >----------------------------------------------------*/


/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float **objects,      /* in: [numObjs][numCoords] */
                         int **mask,      /* in: [numObjs][numCoords] */
                         int oIndex,
                         float **clusters,     /* out: [numClusters][numCoords] */
                         int **cmask,     /* out: [numClusters][numCoords] */
                         float *weight,
                         char cdist
                         )
{
    int   index, i;
    float dist, min_dist;


    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = distancematrix(numCoords, objects, clusters, mask, cmask, weight, oIndex, 0, cdist, 0);

    for (i=1; i<numClusters; i++) {
        dist = distancematrix(numCoords, objects, clusters, mask, cmask, weight, oIndex, i, cdist, 0);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}


/* Some sample C code for the quickselect algorithm,
   taken from Numerical Recipes in C. */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

static float quickselect(float *arr, int n, int k) {
  unsigned long i,ir,j,l,mid;
  float a,temp;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}

static float median (int n, float x[])
/**
 * TODO: Running quick sort and taking median is too slow so changed
 * to use quick select as it has linear run time of O(n) or O(n^2) worst case
 *
 * Find the median of X(1), ... , X(N), using as much of the quicksort
 * algorithm as is needed to isolate it.
 * Uses cdpAdvancedQuickSort algorithm from nvidia 7.5 sample code
*/

{ int i, j;
  int nr = n / 2;
  int nl = nr - 1;
  int even = 0;


  if (n==2*nr) even = 1;
  if (n<3)
  { if (n<1) return 0.;
    if (n == 1) return x[0];
    return 0.5*(x[0]+x[1]);
  }

  if(even) return quickselect(x, n, nr);
  else return quickselect(x, n, nl);

//  run_qsort(n, x, 0);
//
//  if (even) return (0.5*(x[nl]+x[nr]));
//  return x[nr];

}



static void getclustermedians( int nrows, int ncolumns, int nclusters, int* clusterid, float** objects, float** dimClusters, float *cache) {

	int i, j, k, index;



	for (i = 0; i < nclusters; i++) {

		for (j = 0; j < ncolumns; j++) {
			int count = 0;
			for (k = 0; k < nrows; k++) {

				if (i == clusterid[k] && objects[k][j] != 0) {
					cache[count] = objects[k][j];
					count++;
				}

			}

			if (count > 0) {
				dimClusters[i][j] = median(count, cache);

			} else {
				dimClusters[i][j] = 0;
			}

		}
	}
}

/*----< seq_kmedians() >-------------------------------------------------------*/
/* return an array of cluster centers of size [numClusters][numCoords]       */
int seq_kmedians(char dist,
			   float **objects,      /* in: [numObjs][numCoords] */
			   int **mask,      /* in: [numObjs][numCoords] */
			   float *weight,
               int     numCoords,    /* no. features */
               int     numObjs,      /* no. objects */
               int     numClusters,  /* no. clusters */
               float   threshold,    /* % objects change membership */
               int    *membership,   /* out: [numObjs] */
               float **clusters     /* out: [numClusters][numCoords] */
               )

{
    int      i, j, index, loop=0;
//    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
//                                new cluster */
    float    delta;          /* % of objects change their clusters */
//    float  **newClusters;    /* [numClusters][numCoords] */

    int **cmask;     /* out: [numClusters][numCoords] */
    malloc2D(cmask, numClusters, numCoords, int);

    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    float *cache;

    cache = (float*) malloc(numObjs*sizeof(float));

    /* need to initialize newClusterSize and newClusters[0] to all 0 */
//    newClusterSize = (int*) calloc(numClusters, sizeof(int));
//    assert(newClusterSize != NULL);

//    newClusters    = (float**) malloc(numClusters * sizeof(float*));
//    assert(newClusters != NULL);
//    newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
//    assert(newClusters[0] != NULL);
//    for (i=1; i<numClusters; i++)
//        newClusters[i] = newClusters[i-1] + numCoords;

    for (i=0; i<numClusters; i++) {
    	for (j=0; j<numCoords; j++) {
    		clusters[i][j] = objects[i][j];
    		cmask[i][j] = clusters[i][j] == 0? 0 : 1;
    	}
    }

    do {
        delta = 0.0;
        for (i=0; i<numObjs; i++) {
            /* find the array index of nestest cluster center */
            index = find_nearest_cluster(numClusters, numCoords, objects, mask, i, clusters, cmask, weight, dist);

            /* if membership changes, increase delta by 1 */
            if (membership[i] != index) delta += 1.0;

            /* assign the membership to object i */
            membership[i] = index;

            /* update new cluster center : sum of objects located within */
//            newClusterSize[index]++;
//            for (j=0; j<numCoords; j++)
//                newClusters[index][j] += objects[i][j];
        }

//        /* average the sum and replace old cluster center with newClusters */
//        for (i=0; i<numClusters; i++) {
//            for (j=0; j<numCoords; j++) {
//                if (newClusterSize[i] > 0)
//                    clusters[i][j] = newClusters[i][j] / newClusterSize[i];
//                	cmask[i][j] = clusters[i][j] == 0? 0 : 1;
//                	newClusters[i][j] = 0.0;   /* set back to 0 */
//            }
//            newClusterSize[i] = 0;   /* set back to 0 */
//        }
            
        delta /= numObjs;

        getclustermedians(numObjs, numCoords, numClusters, membership, objects, clusters, cache);
        for (i=0; i<numClusters; i++) {
          	for (j=0; j<numCoords; j++) {
          		cmask[i][j] = clusters[i][j] == 0? 0 : 1;
          	}
          }

    } while (delta > threshold && loop++ < 500);

//    free(newClusters[0]);
//    free(newClusters);
//    free(newClusterSize);

    return 1;
}

