/*
 ============================================================================
 Name        : ensembleClust_pall_cuda.cu
 Author      : Vinay B Gavirangaswamy
 Version     :
 Copyright   :  This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.


    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.
 Description : CUDA related code goes into this file. You can customize it according to project needs. 
 ============================================================================
 */

#include <iostream>
#include <numeric>
#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <cassert>
#include <cstring>
#include <cassert>



#include "common/clustLib.h"
#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "common/distCalcMthds.h"
#include "timnugent/spectral/spectral.h"
#include "northwestern/ece/wkliao/kmeans.h"
#include "northwestern/ece/wkliao/kmedians.h"


bool isPow2Host(unsigned int x)
{
    return ((x&(x-1))==0);
}

unsigned int nextPow2Host(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}
////////////////////////////////////////////////////////////////////////////////
// Compute the number of threads and blocks to use for the given reduction kernel
// We set threads / block to the minimum of maxThreads and
////////////////////////////////////////////////////////////////////////////////
void getNumBlocksAndThreadsHost(/*int whichKernel,*/ int n, int maxBlocks,
	int maxThreads, int &blocks, int &threads, int *maxGridSize, int *maxThreadsPerBlock) {

	threads = (n < maxThreads) ? nextPow2Host(n) : maxThreads;
	blocks = (n + threads - 1) / threads;
}


/* ========================================================================= */

void seq_spectral_adapter(int nclusters, int nrows, int ncols, float* data1d, char method, char dist, int *clusterid)

{

//	int i/*, j, k*/;
	int transpose = 0;
	const int nelements = (transpose == 0) ? nrows : ncols;

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	 Spectral* P = new Spectral();
	 P->read_data(data1d, nrows, ncols);
	 P->set_gamma(172.05); // Set RBF gamma param
	 P->set_centers(nclusters); // Set K = 2 - number of centroids and number of eigenvectors
	 P->cluster(dist); // Cluster using K-means
	 P->write_data(clusterid);
	 delete P;

	return;
}
/* ========================================================================= */
void seq_spectral_adapter_full_distmat(int nrows, int ncols, int numclusters,float *distmatrix, int* idxs)
{
	if(numclusters == 1){
		memset(idxs, 0, nrows*sizeof(int));
		return;
	}
	fastsc_full_distmat(nrows, ncols, numclusters, distmatrix, idxs);
	return;
}
/* ========================================================================= */
void mex2seq_spectral_adapter_full_distmat(int nclusters, int nrows, int ncols, float* data1d,
	char* _method, char* _dist, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;

	// Assign lib specific parameters to link method
	// and distance function
	if(_method == NULL){
			method = 's';
		}
		else
		if (strcmp(_method, LNK_CODE_AVG) == 0) {
			method = 'a';
		}
		else if (strcmp(_method, LNK_CODE_CEN) == 0) {
			method = 'c';
		}
		else if (strcmp(_method, LNK_CODE_COM) == 0) {
			method = 'm';
		}
		else if (strcmp(_method, LNK_CODE_SIN) == 0) {
			method = 's';
		}
		else if (strcmp(_method, LNK_CODE_MED) == 0) {
			method = 'a'; // TODO: median linkage should be implemented later
		}
		else if (strcmp(_method, LNK_CODE_WAR) == 0) {
			method = 'w'; // TODO: ward linkage should be implemented later
		}
		else if (strcmp(_method, LNK_CODE_WEI) == 0) {
			method = 'a'; // TODO: weighted linkage should be implemented later
		}

		if (_dist == NULL) {
				dist = 'e';
		}
		else
		if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
			dist = 'e';
		}
		else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
			dist = 'l';
		}
		else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
			dist = 'b';
		}
		else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
			dist = 'c';
		}
		else if (strcmp(_dist, DIST_MTRC_ACOR) == 0) {
			dist = 'a';
		}
		else if (strcmp(_dist, DIST_MTRC_UCOR) == 0) {
			dist = 'u';
		}
		else if (strcmp(_dist, DIST_MTRC_AUCOR) == 0) {
			dist = 'x';
		}
		else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
			dist = 'o';
		}
		else if (strcmp(_dist, DIST_MTRC_KEN) == 0) {
			dist = 'k';
		}
		else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
			dist = 'm';
		}
		else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
			dist = 'j';
		}
		else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
			dist = 'h';
		}
		else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
			dist = 's';
		}
		else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
			dist = 'g';
		}

	seq_spectral_adapter_full_distmat(nrows, ncols, nclusters >= nrows ? nrows : nclusters, data1d, clusterid);

}
/* ========================================================================= */
void mex2seq_spectral_adapter(int nclusters, int nrows, int ncols, float* data1d,
	char* _method, char* _dist, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;

	// Assign lib specific parameters to link method
	// and distance function
	if(_method == NULL){
		method = 's';
	}
	else
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	}
	else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	}
	else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	}
	else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	}
	else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'w'; // TODO: ward linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (_dist == NULL) {
			dist = 'e';
	}
	else
	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'l';
	}
	else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	}
	else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	}
	else if (strcmp(_dist, DIST_MTRC_ACOR) == 0) {
		dist = 'a';
	}
	else if (strcmp(_dist, DIST_MTRC_UCOR) == 0) {
		dist = 'u';
	}
	else if (strcmp(_dist, DIST_MTRC_AUCOR) == 0) {
		dist = 'x';
	}
	else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	}
	else if (strcmp(_dist, DIST_MTRC_KEN) == 0) {
		dist = 'k';
	}
	else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	}
	else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	}
	else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	}
	else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	}
	else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	int** mask;
	float **data2d;
	malloc2D(mask, nrows, ncols, int);

	malloc2D(data2d, nrows, ncols, float);
	memcpy(data2d[0], data1d, nrows * ncols * sizeof(float));

	for (int row = 0; row < nrows; row++){
		for (int col = 0; col < ncols; col++){
			mask[row][col] = data2d[row][col] == 0? 0 : 1;
		}
	}

	float weight[ncols];

	for (int col = 0; col < ncols; col++){
		weight[col] = 1;
	}

	float **distmatrix =
	      distancematrix(nrows, ncols, data2d, mask, weight, dist, 0);
	    if (!distmatrix) return NULL; /* Insufficient memory */

	free(data2d[0]);
	free(data2d);
	free(mask[0]);
	free(mask);

	fastsc(nrows, ncols,nclusters, distmatrix, clusterid);

	 int i;
	 for (i = 1; i < nrows; i++) free(distmatrix[i]);
	 free (distmatrix);

//	seq_spectral_adapter(nclusters >= nrows ? nrows : nclusters, nrows, ncols, data1d, method, dist, clusterid);



}

/* ========================================================================= */

void mex2seq_kmeans_adapter(int nclusters, int nrows, int ncols, float** data2d,int **mask, float *weight,
	char* _method, char* _dist, float threshold, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;

	// Assign lib specific parameters to link method
	// and distance function
	if(_method == NULL){
		method = 's';
	}
	else
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	}
	else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	}
	else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	}
	else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	}
	else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'w'; // TODO: ward linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (_dist == NULL) {
			dist = 'e';
	}
	else
	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'l';
	}
	else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	}
	else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	}
	else if (strcmp(_dist, DIST_MTRC_ACOR) == 0) {
		dist = 'a';
	}
	else if (strcmp(_dist, DIST_MTRC_UCOR) == 0) {
		dist = 'u';
	}
	else if (strcmp(_dist, DIST_MTRC_AUCOR) == 0) {
		dist = 'x';
	}
	else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	}
	else if (strcmp(_dist, DIST_MTRC_KEN) == 0) {
		dist = 'k';
	}
	else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	}
	else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	}
	else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	}
	else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	}
	else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	float **clusters;
	malloc2D(clusters, nclusters, ncols, float);

//	printf("\n5. Calling CUDA Spectral Clustering Adapter...\n");
	seq_kmeans(dist,
				   data2d,      /* in: [numObjs][numCoords] */
				   mask,      /* in: [numObjs][numCoords] */
				   weight,
	               ncols,    /* no. features */
	               nrows,      /* no. objects */
	               nclusters,  /* no. clusters */
	               threshold,    /* % objects change membership */
	               clusterid,   /* out: [numObjs] */
	               clusters     /* out: [numClusters][numCoords] */
	               );


}

/* ========================================================================= */

void mex2seq_kmedians_adapter(int nclusters, int nrows, int ncols, float** data2d,int **mask, float *weight,
	char* _method, char* _dist, float threshold, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;

	// Assign lib specific parameters to link method
	// and distance function
	if(_method == NULL){
		method = 's';
	}
	else
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	}
	else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	}
	else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	}
	else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	}
	else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'w'; // TODO: ward linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (_dist == NULL) {
			dist = 'e';
	}
	else
	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'l';
	}
	else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	}
	else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	}
	else if (strcmp(_dist, DIST_MTRC_ACOR) == 0) {
		dist = 'a';
	}
	else if (strcmp(_dist, DIST_MTRC_UCOR) == 0) {
		dist = 'u';
	}
	else if (strcmp(_dist, DIST_MTRC_AUCOR) == 0) {
		dist = 'x';
	}
	else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	}
	else if (strcmp(_dist, DIST_MTRC_KEN) == 0) {
		dist = 'k';
	}
	else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	}
	else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	}
	else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	}
	else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	}
	else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	float **clusters;
	malloc2D(clusters, nclusters, ncols, float);

//	printf("\n5. Calling CUDA Spectral Clustering Adapter...\n");
	seq_kmedians(dist,
				   data2d,      /* in: [numObjs][numCoords] */
				   mask,      /* in: [numObjs][numCoords] */
				   weight,
	               ncols,    /* no. features */
	               nrows,      /* no. objects */
	               nclusters,  /* no. clusters */
	               threshold,    /* % objects change membership */
	               clusterid,   /* out: [numObjs] */
	               clusters     /* out: [numClusters][numCoords] */
	               );


}

/* ========================================================================= */
void mex2seq_gmm_adapter(int nclusters, int nrows, int ncols, float* data1d,
	char* _method, char* _dist, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;

	// Assign lib specific parameters to link method
	// and distance function
	if(_method == NULL){
		method = 's';
	}
	else
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	}
	else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	}
	else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	}
	else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	}
	else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'w'; // TODO: ward linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (_dist == NULL) {
			dist = 'e';
	}
	else
	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'l';
	}
	else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	}
	else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	}
	else if (strcmp(_dist, DIST_MTRC_ACOR) == 0) {
		dist = 'a';
	}
	else if (strcmp(_dist, DIST_MTRC_UCOR) == 0) {
		dist = 'u';
	}
	else if (strcmp(_dist, DIST_MTRC_AUCOR) == 0) {
		dist = 'x';
	}
	else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	}
	else if (strcmp(_dist, DIST_MTRC_KEN) == 0) {
		dist = 'k';
	}
	else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	}
	else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	}
	else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	}
	else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	}
	else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}



//	printf("\n5. Calling CUDA Spectral Clustering Adapter...\n");

	seq_gmm_main(nclusters, data1d, ncols, nrows, clusterid);


}

/* ========================================================================= */

void mex2seq_agglomerative_adapter(int nclusters, int nrows, int ncols, float** data1d,
	char* _method, char* _dist, int *clusterid)
	/* Perform hierarchical clustering on data */
{

	if(nclusters == 1){
		memset(clusterid,0, nrows*sizeof(int));
		return;
	}

	char method, dist;
	//int i, j;


	// Assign lib specific parameters to link method
	// and distance function
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	}
	else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	}
	else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	}
	else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	}
	else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'w'; // TODO: ward linkage should be implemented later
	}
	else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'e';
	}
	else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	}
	else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	}
	else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	}
	else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	}
	else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	}
	else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	}
	else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	}
	else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	agglom_adapter(nclusters >= nrows ? nrows : nclusters, nrows, ncols, data1d, method, dist, clusterid);


}

