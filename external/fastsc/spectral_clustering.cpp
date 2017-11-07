#include <iostream>
#include "rsymsol.h"
#include "ardsmat.h"
#include "ardssym.h"
#include "ardnsmat.h"
#include "ardsnsym.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include<cassert>

#ifdef __cplusplus
extern "C"
{
#endif
#include<cblas.h>
#ifdef __cplusplus
}
#endif



#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "northwestern/ece/wkliao/kmedians.h"

using namespace std;

int CPU_MULT(float *x, float *y, int n, int nnz, int* rowPtr, int* colIndex, float* val){
	float fone = 1.0;
	float fzero = 0.0;

	// y = Ax
	cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, fone, val, n, x, 1, fzero, y, 1);



	return 0;

}


void random_labels(std::vector<int>& labels, int n, int k) {

	for(int i = 0; i < n; i++) {
		labels[i] = rand() % k;
	}

}

void regular_labels(std::vector<int>& labels, int n, int k) {
	// Initialize by assigning nodes that are close in indexing order with the same label.
	int l = n/k;
	int count = 0;
	int cur = 0;
	for(int i = 0; i < n; i++) {
		labels[i] = cur;
		count++;
		if(count > l) {
			cur++;
			count = 0;
		}
	}
}

void fastsc(int nrows, int ncols, int numclusters,float **distmatrix, int* idxs)
{

	int n = nrows;
	int k = numclusters;
	int nnz = ((nrows * (nrows+1))/2)-nrows;


	std::vector<int> row(nnz), col(nnz);

	// Initialize the degree
	std::vector<float> degree(n, 0.0);

	// For unweighted graphs, edge weights are initilized to 1.0. Otherwise, revise the code to the specific graph representation.
	std::vector<float> val(nnz, 1.0);

	int count = 0;

	for(int crow=0; crow<nrows; crow++){
		for(int ccol=0; ccol<crow; ccol++){

				row[count] = crow;
				col[count] = ccol;
				val[count] =  distmatrix[crow][ccol];
				degree[row[count]] = degree[row[count]] + val[count];
//				cout<<"[" <<count<<","<<crow<<","<<ccol<<"] "<< row[count]<< " "<< col[count]<< " "<< val[count]<< " "<< degree[row[count]]<< " "<< distmatrix[crow][ccol]<<" "<<endl;
				count++;


//				row[count] = ccol;
//				col[count] = crow;
//				val[count] =  distmatrix[crow][ccol];
//				degree[row[count]] = degree[row[count]] + val[count];
//				cout<<"[" <<count<<","<<crow<<","<<ccol<<"] "<< row[count]<< " "<< col[count]<< " "<< val[count]<< " "<< degree[row[count]]<< " "<< distmatrix[crow][ccol]<<" "<<endl;
//				count++;




		}

	}


//		cout<<"Start computing normalized Graph Laplacian..."<<endl;
//		for(int i = 0; i < n; ++i) {
//			if (degree[i] < 1e-8) {
//				cout<<"Node " <<i<<" is an isolated node"<<endl;
//				cout<<"Please eliminate isolated nodes and try again!"<<endl;
//				exit(1);
//			}
//		}

		std::vector<float> degree_sqrt(n);

		// Normlize the edge weight of <i, j> by 1.0/sqrt(degree[i] * degree[j])
		for(int i = 0; i < n; ++i) {
			degree_sqrt[i] = sqrt(degree[i]);
		}

		for(int i = 0; i < nnz; ++i) {
			val[i] = val[i] / (degree_sqrt[col[i]] * degree_sqrt[row[i]]);
		}


		  int*    irow = row.data();       // pointer to an array that stores the row
		                      	  	  	   // indices of the nonzeros in A.
		  int*    pcol = col.data();       // pointer to an array of pointers to the
		  	  	  	  	  	  	  	  	   // beginning of each column of A in vector A.
		  float* A = val.data();          			   // pointer to an array that stores the
		  	  	  	  	  	   	   	   	   // nonzero elements of A.


		  ARrcSymStdEig<float> prob(n, k, "LM");
		  while (!prob.ArnoldiBasisFound()) {
			  prob.TakeStep();
			  if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {
				  CPU_MULT(prob.GetVector(), prob.PutVector(), n, nnz, irow, pcol, A);
			  }
		  }

		  // Finding eigenvalues and eigenvectors.
		  prob.FindEigenvectors();



//		cout<<"Completed computing the first smallest k eigenvectors!"<<endl;

		// Extract eigenvectors.
		// Rearrange the order such that values between i * k and (i+1)*k-1 are eigenmap for node indexed by i
//		cout<<"Start kmeans clustering algorithm on the k eigenvectors..."<<endl;
		// TODO: Re-write kmeans as this implementation is very slow! :0

		/*
		 * Variables to keep track of memory used on device.
		 * Depending on the current memory usage we need to restrict
		 * number of threads.
		 */
		ncols = k;


		float* weight;
		float **eigen_objects_host;
		int** mask;

		size_t memDataSz = 0, memWtSz = 0, memIdxSz=0;

		memDataSz = nrows * ncols * sizeof(float);

		malloc2D(mask, nrows, ncols, int);

		malloc2D(eigen_objects_host, nrows, ncols, float);
		assert(eigen_objects_host != NULL);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < k; ++j) {
				eigen_objects_host[i][j] = prob.Eigenvector(j, i);
			}
		}

		for (int row = 0; row < nrows; row++){
				for (int col = 0; col < ncols; col++){
					mask[row][col] = eigen_objects_host[row][col] == 0? 0 : 1;
				}
			}

		memWtSz = ncols * sizeof(float);

		weight = (float*)malloc(memWtSz);

		assert(weight != NULL);

		for(int i=0; i<ncols; i++)
			weight[i] = 1;


		float **clusters;
			malloc2D(clusters, k, ncols, float);
			//Apply K-means algorithm on the eigenvectors
			int iterations = 100;
			double threshold=0.000001;

			seq_kmedians('e',
						  eigen_objects_host,      /* in: [numObjs][numCoords] */
						  mask,      /* in: [numObjs][numCoords] */
						  weight,
			              ncols,    /* no. features */
			              nrows,      /* no. objects */
			              k,  /* no. clusters */
			              threshold,    /* % objects change membership */
			              idxs,   /* out: [numObjs] */
			              clusters     /* out: [numClusters][numCoords] */
			              );

		return;
}

void fastsc_full_distmat(int nrows, int ncols, int numclusters,float *distmatrix, int* idxs)
{

	int n = nrows;
	int k = numclusters;
	int nnz = ((nrows * (nrows+1))/2)-nrows;


	std::vector<int> row(nnz), col(nnz);

	// Initialize the degree
	std::vector<float> degree(n, 0.0);

	// For unweighted graphs, edge weights are initilized to 1.0. Otherwise, revise the code to the specific graph representation.
	std::vector<float> val(nnz, 1.0);

	int count = 0;

	for(int crow=0; crow<nrows; crow++){
		for(int ccol=0; ccol<crow; ccol++){

				row[count] = crow;
				col[count] = ccol;
				val[count] =  distmatrix[crow * nrows + ccol];
				degree[row[count]] = degree[row[count]] + val[count];
//				cout<<"[" <<count<<","<<crow<<","<<ccol<<"] "<< row[count]<< " "<< col[count]<< " "<< val[count]<< " "<< degree[row[count]]<< " "<< distmatrix[crow][ccol]<<" "<<endl;
				count++;


//				row[count] = ccol;
//				col[count] = crow;
//				val[count] =  distmatrix[crow][ccol];
//				degree[row[count]] = degree[row[count]] + val[count];
//				cout<<"[" <<count<<","<<crow<<","<<ccol<<"] "<< row[count]<< " "<< col[count]<< " "<< val[count]<< " "<< degree[row[count]]<< " "<< distmatrix[crow][ccol]<<" "<<endl;
//				count++;




		}

	}

	//		cout<<"Start computing normalized Graph Laplacian..."<<endl;
	//		for(int i = 0; i < n; ++i) {
	//			if (degree[i] < 1e-8) {
	//				cout<<"Node " <<i<<" is an isolated node"<<endl;
	//				cout<<"Please eliminate isolated nodes and try again!"<<endl;
	//				exit(1);
	//			}
	//		}

	std::vector<float> degree_sqrt(n);

	// Normlize the edge weight of <i, j> by 1.0/sqrt(degree[i] * degree[j])
	for(int i = 0; i < n; ++i) {
		degree_sqrt[i] = sqrt(degree[i]);
	}

	for(int i = 0; i < nnz; ++i) {
		val[i] = val[i] / (degree_sqrt[col[i]] * degree_sqrt[row[i]]);
	}


	int*    irow = row.data();       // pointer to an array that stores the row
	// indices of the nonzeros in A.
	int*    pcol = col.data();       // pointer to an array of pointers to the
	// beginning of each column of A in vector A.
	float* A = val.data();          			   // pointer to an array that stores the
	// nonzero elements of A.


//	ARdsSymMatrix<float> matrix(n, A, 'L');

	// Defining what we need: the k eigenvectors of A with largest magnitude.
//	ARluSymStdEig<float> dprob(k, matrix, "LM");

	// Finding eigenvalues and eigenvectors.
//	dprob.FindEigenvectors();
	// Printing eigenvalue solution.
	// Solution(prob);

	ARrcSymStdEig<float> prob(n, k, "LM");
	while (!prob.ArnoldiBasisFound()) {
		prob.TakeStep();
		if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {
			CPU_MULT(prob.GetVector(), prob.PutVector(), n, nnz, irow, pcol, A);
		}
	}

	// Finding eigenvalues and eigenvectors.
	prob.FindEigenvectors();


	//		cout<<"Completed computing the first smallest k eigenvectors!"<<endl;

	// Extract eigenvectors.
	// Rearrange the order such that values between i * k and (i+1)*k-1 are eigenmap for node indexed by i
	//		cout<<"Start kmeans clustering algorithm on the k eigenvectors..."<<endl;
	// TODO: Re-write kmeans as this implementation is very slow! :0

	/*
	 * Variables to keep track of memory used on device.
	 * Depending on the current memory usage we need to restrict
	 * number of threads.
	 */
	ncols = k;


	float* weight;
	float **eigen_objects_host;
	int** mask;

	size_t memDataSz = 0, memWtSz = 0, memIdxSz=0;

	memDataSz = nrows * ncols * sizeof(float);

	malloc2D(mask, nrows, ncols, int);

	malloc2D(eigen_objects_host, nrows, ncols, float);
	assert(eigen_objects_host != NULL);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < k; ++j) {
			eigen_objects_host[i][j] = prob.Eigenvector(j, i);
		}
	}

	for (int row = 0; row < nrows; row++){
		for (int col = 0; col < ncols; col++){
			mask[row][col] = eigen_objects_host[row][col] == 0? 0 : 1;
		}
	}

	memWtSz = ncols * sizeof(float);

	weight = (float*)malloc(memWtSz);

	assert(weight != NULL);

	for(int i=0; i<ncols; i++)
		weight[i] = 1;


	float **clusters;
	malloc2D(clusters, k, ncols, float);
	//Apply K-means algorithm on the eigenvectors
	int iterations = 100;
	double threshold=0.000001;

	seq_kmedians('e',
			eigen_objects_host,      /* in: [numObjs][numCoords] */
			mask,      /* in: [numObjs][numCoords] */
			weight,
			ncols,    /* no. features */
			nrows,      /* no. objects */
			k,  /* no. clusters */
			threshold,    /* % objects change membership */
			idxs,   /* out: [numObjs] */
			clusters     /* out: [numClusters][numCoords] */
	);

		return;
}

void hello_fastsc(){
	printf("\n~~~HELLO FROM FAST SC~~~\n");
}
