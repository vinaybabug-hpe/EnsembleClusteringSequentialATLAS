// Spectral clustering, by Tim Nugent 2014

#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

class Spectral{

public:
	Spectral() : centers(2), kernel_type(1), normalise(1), max_iters(1000), gamma(0.001), constant(1.0), order(2.0) {}
	explicit Spectral(MatrixXd& d) : centers(2), kernel_type(1), normalise(1), max_iters(1000), gamma(0.001), constant(1.0), order(2.0) {X = d;}
	int read_data(float *data1d, int nrows, int ncols);
	int write_data(int *clusterid);
	void set_centers(const unsigned int i){centers = i;};
	void set_kernel(const unsigned int i){kernel_type = i;};	
	void set_normalise(const unsigned int i){normalise = i;};
	void set_gamma(const double i){gamma = i;};
	void set_constant(const double i){constant = i;};
	void set_order(const double i){order = i;};
	void set_max_iters(const unsigned int i){max_iters = i;};
	void cluster();
	void cluster(char dist);
	const vector<int> &get_assignments() const {return assignments;};
private:
	void generate_kernel_matrix();
	void generate_kernel_matrix(char dist);
	double kernel(const VectorXd& a, const VectorXd& b);
	double kernel(const VectorXd& a, const VectorXd& b, int row1, int row2, char dist);
	void eigendecomposition();
	void kmeans();
	MatrixXd X, K, eigenvectors;
	VectorXd eigenvalues, cumulative;
	unsigned int centers, kernel_type, normalise, max_iters;
	double gamma, constant, order;
	vector<int> assignments;
};

#endif
