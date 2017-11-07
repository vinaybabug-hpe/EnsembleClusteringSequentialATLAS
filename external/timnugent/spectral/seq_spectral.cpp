// Spectral clustering, by Tim Nugent 2014
#include<cmath>
#include<cstdlib>

#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <map>
#include <math.h>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "timnugent/spectral/spectral.h"
#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "common/distCalcMthds.h"




using namespace Eigen;
using namespace std;
int Spectral::read_data(float *data1d, int nrows, int ncols){

	X.resize(nrows, ncols);
	// Read data
	for(int row=0; row < nrows; row++){
		for(int col = 0; col < ncols; col++){
			X(row,col) = data1d[row * ncols + col];
		}
	}

	//X = Map<MatrixXd>( data1d, nrows, ncols );

	return(1);

}

int Spectral::write_data(int *clusterid){


		for(unsigned int i = 0; i < X.rows(); i++){

			clusterid [i] = assignments[i];

		}	

		return(1);

}
/**
 * Returns distance or affinity between 2 objects
 */
double Spectral::kernel(const VectorXd& a, const VectorXd& b){

	switch(kernel_type){
	    case 2  :
	    	return(pow(a.dot(b)+constant,order));
	    default : 
	    	return(exp(-gamma*((a-b).squaredNorm())));
	}

}

/**
 * Returns distance or affinity between 2 objects
 */
//double Spectral::kernel(const VectorXd& a, const VectorXd& b, int row1, int row2, char dist){
//
//	int* mask;
//	float *weight;
//	float* cdata;
//	float* data;
//	int* cmask;
//	double distance;
//
//	int nrows = a.rows() > b.rows()? b.rows(): a.rows();
//	int ncolumns = a.cols() > b.cols()? b.cols(): a.cols();
//
////	cout << "nrows "<< nrows  << " ";
////	cout << "ncolumns "<< ncolumns  << " ";
//
//	mask = (int*) malloc(nrows * ncolumns * sizeof(int));
//	cmask = (int*) malloc(nrows * ncolumns * sizeof(int));
//
//	data = (float*) malloc(nrows * ncolumns * sizeof(float));
//	cdata = (float*) malloc(nrows * ncolumns * sizeof(float));
//
//	weight = (float*) malloc(nrows * sizeof(float));
//
//	for(int count = 0; count < nrows; count++)
//	{
//		data[count] = a[count];
//		cdata[count] = b[count];
//		mask[count] = data[count] == 0? 0: 1;
//		cmask[count] = cdata[count] == 0? 0: 1;
//		weight[count] = 1;
//	}
//
//
//	  distance = calcDistMetricCPU(dist, ncolumns, row1, row2, data, cdata, mask, cmask, weight, 0, 0, 0);
//
//	  return distance;
//}

void Spectral::generate_kernel_matrix(){

	// Fill kernel matrix
	K.resize(X.rows(),X.rows());
	for(unsigned int i = 0; i < X.rows(); i++){
		for(unsigned int j = i; j < X.rows(); j++){
			K(i,j) = K(j,i) = kernel(X.row(i),X.row(j));
//			if(i == 0) cout << K(i,j) << " ";

		}	
	}

	// Normalise kernel matrix	
	VectorXd d = K.rowwise().sum();
	for(unsigned int i = 0; i < d.rows(); i++){
		d(i) = 1.0/sqrt(d(i));
	}
	auto F = d.asDiagonal();
	MatrixXd l = (K * F);
	for(unsigned int i = 0; i < l.rows(); i++){
		for(unsigned int j = 0; j < l.cols(); j++){
			l(i,j) = l(i,j) * d(i);
		}
	}		
	K = l;

}

void Spectral::generate_kernel_matrix(char dist){

	int blockSize;   // The launch configurator returned block size
	int minGridSize; // The minimum grid size needed to achieve the
	// maximum occupancy for a full device launch
	int gridSize;    // The actual grid size needed, based on input size

	float** distmatrix = NULL;
	int **mask;
	float *weight;
	float **data1d;

	int nrows = X.rows(), ncolumns = X.cols();

	weight = (float*)malloc(ncolumns * sizeof(float*));

	malloc2D(mask, nrows, ncolumns, int);

	malloc2D(data1d, nrows, ncolumns, float);


	for (int row = 0; row < nrows; row++){
		mask[row] = (int*)malloc(ncolumns* sizeof(int));
		for (int col = 0; col < ncolumns; col++){
			data1d[row][col] = X(row,col);
			mask[row][col] = data1d[row][col] == 0? 0 : 1;
		}
	}

	assert(weight != NULL);

	for (int i = 0; i < ncolumns; i++)
		weight[i] = 1.0;

	distmatrix =
		      distancematrix(nrows, ncolumns, data1d, mask, weight, dist, 0);
		    if (!distmatrix) return NULL; /* Insufficient memory */


	double *resultC;                // NULL pointer

	// to convert data from MatrixXd to double* column major
	  //Map<MatrixXd>( resultC, resultEigen.rows(), resultEigen.cols() ) =   resultEigen;


	// Fill kernel matrix
	K.resize(X.rows(),X.rows());
	for(unsigned int i = 0; i < X.rows(); i++){
		for(unsigned int j = 0; j < i; j++){
			K(i,j) = K(j,i) =  distmatrix[i][j];//kernel(X.row(i),X.row(j), i, j,dist);
			//if(i == 0) cout << K(i,j) << " ";

		}
	}

	// Normalise kernel matrix
	VectorXd d = K.rowwise().sum();
	for(unsigned int i = 0; i < d.rows(); i++){
		d(i) = 1.0/sqrt(d(i));
	}
	auto F = d.asDiagonal();
	MatrixXd l = (K * F);
	for(unsigned int i = 0; i < l.rows(); i++){
		for(unsigned int j = 0; j < l.cols(); j++){
			l(i,j) = l(i,j) * d(i);
		}
	}
	K = l;

	free(weight);
	free(mask);
	for (int i = 1; i < nrows; i++) free(distmatrix[i]);
	free (distmatrix);
	free(data1d[0]);
	free(data1d);

}

void Spectral::eigendecomposition(){

	EigenSolver<MatrixXd> edecomp(K);
	eigenvalues = edecomp.eigenvalues().real();
	eigenvectors = edecomp.eigenvectors().real();
	cumulative.resize(eigenvalues.rows());
	vector<pair<double,VectorXd> > eigen_pairs; 
	double c = 0.0; 
	for(unsigned int i = 0; i < eigenvectors.cols(); i++){
		if(normalise){
			double norm = eigenvectors.col(i).norm();
			eigenvectors.col(i) /= norm;
		}
		eigen_pairs.push_back(make_pair(eigenvalues(i),eigenvectors.col(i)));
	}
	// http://stackoverflow.com/questions/5122804/sorting-with-lambda
	sort(eigen_pairs.begin(),eigen_pairs.end(), [](const pair<double,VectorXd> a, const pair<double,VectorXd> b) -> bool {return (a.first > b.first);} );

	if(centers > eigen_pairs.size()) centers = eigen_pairs.size();

	for(unsigned int i = 0; i < eigen_pairs.size(); i++){	
		eigenvalues(i) = eigen_pairs[i].first;
		c += eigenvalues(i);
		cumulative(i) = c;
		eigenvectors.col(i) = eigen_pairs[i].second;
	}

	/*
	cout << "Sorted eigenvalues:" << endl;
	for(unsigned int i = 0; i < eigenvalues.rows(); i++){
		if(eigenvalues(i) > 0){
			cout << "PC " << i+1 << ": Eigenvalue: " << eigenvalues(i);
			printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",eigenvalues(i)/eigenvalues.sum(),cumulative(i)/eigenvalues.sum());
			//cout << eigenvectors.col(i) << endl;
		}
	}
	cout << endl;
	*/
	MatrixXd tmp = eigenvectors;
	
	// Select top K eigenvectors where K = centers
	eigenvectors = tmp.block(0,0,tmp.rows(),centers);

}	

void Spectral::cluster(char dist){

	generate_kernel_matrix(dist);
	eigendecomposition();
	kmeans();

}

void Spectral::cluster(){

	generate_kernel_matrix();
	eigendecomposition();
	kmeans();

}

// Code adapted from https://github.com/pthimon/clustering
void Spectral::kmeans(){

	random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> rand_index(0,eigenvectors.rows()-1);
	MatrixXd centroids = MatrixXd::Zero(centers, eigenvectors.cols());
	MatrixXd old_centroids;

	// Randomly select unique centroids 
	vector<int> rands;
	while(rands.size() < centers){
		int r = rand_index(gen);
		bool tag = false;
		for (unsigned int j=0; j < rands.size(); ++j){
			if(r == rands[j]){
				tag = true;
				break;
			}
		}
		if(!tag){				
			centroids.row(rands.size()) = eigenvectors.row(r);
			rands.push_back(r);	
		}
	}

	MatrixXd id = MatrixXd::Identity(centers, centers);
	VectorXd minvals(eigenvectors.rows());
	// Matrix to map vectors to centroids
	MatrixXd post(eigenvectors.rows(), centers);

	int r, c;
	double old_e = 0;
	for (unsigned int n=0; n < max_iters; n++){
		old_centroids = centroids;

		// Calculate distances
		MatrixXd d2(eigenvectors.rows(), centers);
		for (unsigned int j = 0; j < centers; j++){
			for(int k=0; k < eigenvectors.rows(); k++) {
				d2(k,j) = (eigenvectors.row(k)-centroids.row(j)).squaredNorm();
			}
		}
		
		// Assign to nearest centroid
		for (unsigned int k = 0; k < eigenvectors.rows(); k++){
			// Get index of centroid
			minvals[k] = d2.row(k).minCoeff(&r, &c);
			// Set centroid
			post.row(k) = id.row(c);
		}

		// Adjust centeroids
		VectorXd num_points = post.colwise().sum();
		for(unsigned int j = 0; j < centers; j++){
			if(num_points(j) > 0) {
				MatrixXd s = MatrixXd::Zero(1,eigenvectors.cols());
				for(unsigned int k = 0; k < eigenvectors.rows(); k++){
					if(post(k,j) == 1){
						s += eigenvectors.row(k);
					}
				}
				centroids.row(j) = s/num_points[j];
			}
		}

		// Calculate error - total squared distance from centroids
		double e = minvals.sum();
		double ediff = fabs(old_e-e);
		double cdiff = (centroids-old_centroids).cwiseAbs().maxCoeff();
		// TODO: to check error
//		printf("Iterations %i : Error %2.4f : Error delta %2.4f : Centroid movement %2.4f\n",n+1,e,ediff,cdiff);
		if(n && cdiff < numeric_limits<double>::epsilon() && ediff < numeric_limits<double>::epsilon()){
			break;
		}
		old_e = e;
	}

	map<int,int> data_to_cluster;
	for (unsigned int j = 0; j < centers; j++){
		for (int k=0; k < eigenvectors.rows(); k++) {
			if (post(k,j) == 1) {
				data_to_cluster[k] = j+1;
			}
		}
	}
	for (int k=0; k < eigenvectors.rows(); k++){
		assignments.push_back(data_to_cluster[k]);
	}

}
