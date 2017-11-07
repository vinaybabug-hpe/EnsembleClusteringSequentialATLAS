/*
 * cuda_gmm.cu
 *
 *  Created on: Mar 23, 2016
 *      Author: vinaya
 *	   Version:
 *	 Copyright: This program is free software: you can redistribute it and/or modify
 *   			it under the terms of the GNU General Public License as published by
 *   			the Free Software Foundation, either version 3 of the License, or
 *   			(at your option) any later version.
 *
 *    			This program is distributed in the hope that it will be useful,
 *    			but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    			GNU General Public License for more details.
 *
 *
 *    			You should have received a copy of the GNU General Public License
 *  			along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * Description: 
 */

#include <cstdio>
#include <cfloat>

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h> // for clock(), clock_t, CLOCKS_PER_SEC

// includes, project
#include "rochester/gmm/gaussian.h"
#include "rochester/gmm/invert_matrix.h"


static void seed_clusters(float *data, clusters_t* clusters, int D, int M, int N) {
    float* variances = (float*) malloc(sizeof(float)*D);
    float* means = (float*) malloc(sizeof(float)*D);

    // Compute means
    for(int d=0; d < D; d++) {
        means[d] = 0.0;
        for(int n=0; n < N; n++) {
            means[d] += data[n*D+d];
        }
        means[d] /= (float) N;
    }

    // Compute variance of each dimension
    for(int d=0; d < D; d++) {
        variances[d] = 0.0;
        for(int n=0; n < N; n++) {
            variances[d] += data[n*D+d]*data[n*D+d];
        }
        variances[d] /= (float) N;
        variances[d] -= means[d]*means[d];
    }

    // Average variance
    float avgvar = 0.0;
    for(int d=0; d < D; d++) {
        avgvar += variances[d];
    }
    avgvar /= (float) D;

    // Initialization for random seeding and uniform seeding
    float fraction;
    int seed;
    if(M > 1) {
        fraction = (N-1.0f)/(M-1.0f);
    } else {
        fraction = 0.0;
    }
    srand(clock());

    for(int m=0; m < M; m++) {
        clusters->N[m] = (float) N / (float) M;
        clusters->pi[m] = 1.0f / (float) M;
        clusters->avgvar[m] = avgvar / COVARIANCE_DYNAMIC_RANGE;

        DEBUG("N: %.2f\tPi: %.2f\tAvgvar: %e\n",clusters->N[m],clusters->pi[m],clusters->avgvar[m]);

        // Choose cluster centers
        DEBUG("Means: ");
        #if UNIFORM_SEED
            for(int d=0; d < D; d++) {
                clusters->means[m*D+d] = data[((int)(m*fraction))*D+d];
                DEBUG("%.2f ",clusters->means[m*D+d]);
            }
        #else
            seed = rand() % N;
            DEBUG("Cluster %d seed = event #%d\n",m,seed);
            for(int d=0; d < D; d++) {
                clusters->means[m*D+d] = data[seed*D+d];
                DEBUG("%.2f ",clusters->means[m*D+d]);
            }
        #endif
        DEBUG("\n");

        // Set covariances to identity matrices
        for(int i=0; i < D; i++) {
            for(int j=0; j < D; j++) {
                if(i == j) {
                    clusters->R[m*D*D+i*D+j] = 1.0f;
                } else {
                    clusters->R[m*D*D+i*D+j] = 0.0f;
                }
            }
        }

        DEBUG("R:\n");
        for(int d=0; d < D; d++) {
            for(int e=0; e < D; e++)
                DEBUG("%.2f ",clusters->R[m*D*D+d*D+e]);
            DEBUG("\n");
        }
        DEBUG("\n");

    }

    free(variances);
    free(means);
}

static void constants(clusters_t* clusters, int M, int D) {
    float log_determinant;
    float* matrix = (float*) malloc(sizeof(float)*D*D);

    float sum = 0.0;
    for(int m=0; m < M; m++) {
        // Invert covariance matrix
        memcpy(matrix,&(clusters->R[m*D*D]),sizeof(float)*D*D);
        invert_cpu(matrix,D,&log_determinant);
        memcpy(&(clusters->Rinv[m*D*D]),matrix,sizeof(float)*D*D);

        // Compute constant
        clusters->constant[m] = -D*0.5f*logf(2.0f*PI) - 0.5f*log_determinant;
        DEBUG("Cluster %d constant: %e\n",m,clusters->constant[m]);

        // Sum for calculating pi values
        sum += clusters->N[m];
    }

    // Compute pi values
    for(int m=0; m < M; m++) {
        clusters->pi[m] = clusters->N[m] / sum;
    }

    free(matrix);
}

static void estep1(float* data, clusters_t* clusters, int D, int M, int N, float* likelihood) {
    clock_t start,finish;
    // Compute likelihood for every data point in each cluster
    float like;
    float* means;
    float* Rinv;
    start = clock();
    for(int m=0; m < M; m++) {
        means = (float*) &(clusters->means[m*D]);
        Rinv = (float*) &(clusters->Rinv[m*D*D]);

        for(int n=0; n < N; n++) {
            like = 0.0;
            #if DIAG_ONLY
                for(int i=0; i < D; i++) {
                    like += (data[i*N+n]-means[i])*(data[i*N+n]-means[i])*Rinv[i*D+i];
                }
            #else
                for(int i=0; i < D; i++) {
                    for(int j=0; j < D; j++) {
                        like += (data[i*N+n]-means[i])*(data[j*N+n]-means[j])*Rinv[i*D+j];
                    }
                }
            #endif
            clusters->memberships[m*N+n] = -0.5f * like + clusters->constant[m] + log(clusters->pi[m]);
        }
    }
    finish = clock();
    DEBUG("estep1: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

static void estep2(float* data, clusters_t* clusters, int D, int M, int N, float* likelihood) {
    clock_t start,finish;
    start = clock();
    float max_likelihood, denominator_sum;
    *likelihood = 0.0f;
    for(int n=0; n < N; n++) {
        // initial condition, maximum is the membership in first cluster
        max_likelihood = clusters->memberships[n];
        // find maximum likelihood for this data point
        for(int m=1; m < M; m++) {
            max_likelihood = fmaxf(max_likelihood,clusters->memberships[m*N+n]);
        }

//        printf("max_likelihood of event %d is : %.3f\n",n,max_likelihood);
        // Computes sum of all likelihoods for this event
        denominator_sum = 0.0f;
        for(int m=0; m < M; m++) {
            denominator_sum += exp(clusters->memberships[m*N+n] - max_likelihood);
        }

//        printf("denominator_sum of event %d is : %.3f likelihood: %.3f\n",n,denominator_sum, *likelihood);

        if(!isnan(denominator_sum)) {
        	denominator_sum = max_likelihood + log(denominator_sum);
        	*likelihood = *likelihood + denominator_sum;
        }
        else{
        	return;
        }

        // Divide by denominator to get each membership
        for(int m=0; m < M; m++) {
            clusters->memberships[m*N+n] = exp(clusters->memberships[m*N+n] - denominator_sum);
//            printf("Membership of event %d in cluster %d: %.3f\n",n,m,clusters->memberships[m*N+n]);
        }
    }
    finish = clock();
    DEBUG("estep2: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

static void mstep_n(float* data, clusters_t* clusters, int D, int M, int N) {
    DEBUG("mstep_n: D: %d, M: %d, N: %d\n",D,M,N);
    for(int m=0; m < M; m++) {
        clusters->N[m] = 0.0;
        // compute effective size of each cluster by adding up soft membership values
        for(int n=0; n < N; n++) {
            clusters->N[m] += clusters->memberships[m*N+n];
        }
    }
}

static void mstep_mean(float* data, clusters_t* clusters, int D, int M, int N) {
    DEBUG("mstep_mean: D: %d, M: %d, N: %d\n",D,M,N);
    for(int m=0; m < M; m++) {
        DEBUG("Cluster %d: ",m);
        for(int d=0; d < D; d++) {
            clusters->means[m*D+d] = 0.0;
            for(int n=0; n < N; n++) {
                clusters->means[m*D+d] += data[d*N+n]*clusters->memberships[m*N+n];
            }
            clusters->means[m*D+d] /= clusters->N[m];
            DEBUG("%f ",clusters->means[m*D+d]);
        }
        DEBUG("\n");
    }
}

static void mstep_covar(float* data, clusters_t* clusters, int D, int M, int N) {
    DEBUG("mstep_covar: D: %d, M: %d, N: %d\n",D,M,N);
    float sum;
    float* means;
    for(int m=0; m < M; m++) {
        means = &(clusters->means[m*D]);
        for(int i=0; i < D; i++) {
            for(int j=0; j <= i; j++) {
                #if DIAG_ONLY
                    if(i != j) {
                        clusters->R[m*D*D+i*D+j] = 0.0f;
                        clusters->R[m*D*D+j*D+i] = 0.0f;
                        continue;
                    }
                #endif
                sum = 0.0;
                for(int n=0; n < N; n++) {
                    sum += (data[i*N+n]-means[i])*(data[j*N+n]-means[j])*clusters->memberships[m*N+n];
                }
                if(clusters->N[m] >= 1.0f) {
                    clusters->R[m*D*D+i*D+j] = sum / clusters->N[m];
                    clusters->R[m*D*D+j*D+i] = sum / clusters->N[m];
                } else {
                    clusters->R[m*D*D+i*D+j] = 0.0f;
                    clusters->R[m*D*D+j*D+i] = 0.0f;
                }
            }
        }
    }
}

static void add_clusters(clusters_t clusters, int c1, int c2, clusters_t temp_cluster, int num_dimensions) {
    float wt1,wt2;

    wt1 = (clusters.N[c1]) / (clusters.N[c1] + clusters.N[c2]);
    wt2 = 1.0f - wt1;

    // Compute new weighted means
    for(int i=0; i<num_dimensions;i++) {
        temp_cluster.means[i] = wt1*clusters.means[c1*num_dimensions+i] + wt2*clusters.means[c2*num_dimensions+i];
    }

    // Compute new weighted covariance
    for(int i=0; i<num_dimensions; i++) {
        for(int j=i; j<num_dimensions; j++) {
            // Compute R contribution from cluster1
            temp_cluster.R[i*num_dimensions+j] = ((temp_cluster.means[i]-clusters.means[c1*num_dimensions+i])
                                                *(temp_cluster.means[j]-clusters.means[c1*num_dimensions+j])
                                                +clusters.R[c1*num_dimensions*num_dimensions+i*num_dimensions+j])*wt1;
            // Add R contribution from cluster2
            temp_cluster.R[i*num_dimensions+j] += ((temp_cluster.means[i]-clusters.means[c2*num_dimensions+i])
                                                    *(temp_cluster.means[j]-clusters.means[c2*num_dimensions+j])
                                                    +clusters.R[c2*num_dimensions*num_dimensions+i*num_dimensions+j])*wt2;
            // Because its symmetric...
            temp_cluster.R[j*num_dimensions+i] = temp_cluster.R[i*num_dimensions+j];
        }
    }

    // Compute pi
    temp_cluster.pi[0] = clusters.pi[c1] + clusters.pi[c2];

    // compute N
    temp_cluster.N[0] = clusters.N[c1] + clusters.N[c2];

    float log_determinant;
    // Copy R to Rinv matrix
    memcpy(temp_cluster.Rinv,temp_cluster.R,sizeof(float)*num_dimensions*num_dimensions);
    // Invert the matrix
    invert_cpu(temp_cluster.Rinv,num_dimensions,&log_determinant);
    // Compute the constant
    temp_cluster.constant[0] = (-num_dimensions)*0.5*logf(2*PI)-0.5*log_determinant;

    // avgvar same for all clusters
    temp_cluster.avgvar[0] = clusters.avgvar[0];
}

static float cluster_distance(clusters_t clusters, int c1, int c2, clusters_t temp_cluster, int num_dimensions) {
    // Add the clusters together, this updates pi,means,R,N and stores in temp_cluster
    add_clusters(clusters,c1,c2,temp_cluster,num_dimensions);

    return clusters.N[c1]*clusters.constant[c1] + clusters.N[c2]*clusters.constant[c2] - temp_cluster.N[0]*temp_cluster.constant[0];
}

static void copy_cluster(clusters_t dest, int c_dest, clusters_t src, int c_src, int num_dimensions) {
    dest.N[c_dest] = src.N[c_src];
    dest.pi[c_dest] = src.pi[c_src];
    dest.constant[c_dest] = src.constant[c_src];
    dest.avgvar[c_dest] = src.avgvar[c_src];
    memcpy(&(dest.means[c_dest*num_dimensions]),&(src.means[c_src*num_dimensions]),sizeof(float)*num_dimensions);
    memcpy(&(dest.R[c_dest*num_dimensions*num_dimensions]),&(src.R[c_src*num_dimensions*num_dimensions]),sizeof(float)*num_dimensions*num_dimensions);
    memcpy(&(dest.Rinv[c_dest*num_dimensions*num_dimensions]),&(src.Rinv[c_src*num_dimensions*num_dimensions]),sizeof(float)*num_dimensions*num_dimensions);
    // do we need to copy memberships?
}

static void writeCluster(FILE* f, clusters_t clusters, int c, int num_dimensions) {
    fprintf(f,"Probability: %f\n", clusters.pi[c]);
    fprintf(f,"N: %f\n",clusters.N[c]);
    fprintf(f,"Means: ");
    for(int i=0; i<num_dimensions; i++){
        fprintf(f,"%f ",clusters.means[c*num_dimensions+i]);
    }
    fprintf(f,"\n");
    fprintf(f,"Memberships: ");
       for(int i=0; i<num_dimensions; i++){
           fprintf(f,"%f ",clusters.memberships[c*num_dimensions+i]);
       }
    fprintf(f,"\n");

    fprintf(f,"\nR Matrix:\n");
    for(int i=0; i<num_dimensions; i++) {
        for(int j=0; j<num_dimensions; j++) {
            fprintf(f,"%f ", clusters.R[c*num_dimensions*num_dimensions+i*num_dimensions+j]);
        }
        fprintf(f,"\n");
    }
    fflush(f);
    /*
    fprintf(f,"\nR-inverse Matrix:\n");
    for(int i=0; i<num_dimensions; i++) {
        for(int j=0; j<num_dimensions; j++) {
            fprintf(f,"%.3f ", c->Rinv[i*num_dimensions+j]);
        }
        fprintf(f,"\n");
    }
    */
}

static void printCluster(clusters_t clusters, int c, int num_dimensions) {
    writeCluster(stdout,clusters,c,num_dimensions);
}

int seq_gmm_main(int desired_num_clusters, float* fcs_data_by_event, int num_dimensions, int num_events, int *clusterIdx){

	int original_num_clusters = desired_num_clusters, stop_number;

	// For profiling the seed kernel
	clock_t seed_start, seed_end, seed_total;

	// For profiling the regroup kernel
	clock_t regroup_start, regroup_end, regroup_total;
	int regroup_iterations = 0;

	// for profiling the reestimate_parameters kernel
	clock_t params_start, params_end, params_total;
	int params_iterations = 0;

	// for profiling the constants kernel
	clock_t constants_start, constants_end, constants_total;
	int constants_iterations = 0;

	// for profiling the GMM order reduction
	clock_t reduce_start, reduce_end, reduce_total;
	int reduce_iterations = 0;

	regroup_total = regroup_iterations = 0;
	params_total = params_iterations = 0;
	constants_total = constants_iterations = 0;
	reduce_total = reduce_iterations = 0;


	// Number of clusters to stop iterating at.
	if (desired_num_clusters == 0) {
		stop_number = 1;
	} else {
		stop_number = desired_num_clusters;
	}

	// Transpose the event data (allows coalesced access pattern in E-step kernel)
	// This has consecutive values being from the same dimension of the data
	// (num_dimensions by num_events matrix)
	float* fcs_data_by_dimension = (float*) malloc(sizeof(float) * num_events * num_dimensions);

	for (int e = 0; e < num_events; e++) {
		for (int d = 0; d < num_dimensions; d++) {
			fcs_data_by_dimension[d * num_events + e] = fcs_data_by_event[e * num_dimensions + d];
		}
	}



	PRINT("Number of events: %d\n", num_events);
	PRINT("Number of dimensions: %d\n\n", num_dimensions);

	PRINT("Starting with %d cluster(s), will stop at %d cluster(s).\n", original_num_clusters, stop_number);



	// Setup the cluster data structures on host
	clusters_t clusters;
	clusters.N = (float*) malloc(sizeof(float) * original_num_clusters);
	clusters.pi = (float*) malloc(sizeof(float) * original_num_clusters);
	clusters.constant = (float*) malloc(sizeof(float) * original_num_clusters);
	clusters.avgvar = (float*) malloc(sizeof(float) * original_num_clusters);
	clusters.means = (float*) malloc(sizeof(float) * num_dimensions * original_num_clusters);
	clusters.R = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	clusters.Rinv = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	clusters.memberships = (float*) malloc(sizeof(float) * num_events * original_num_clusters);
	if (!clusters.means || !clusters.R || !clusters.Rinv || !clusters.memberships) {
		printf("ERROR: Could not allocate memory for clusters.\n");
		return 1;
	}

	// Setup the cluster data structures on host
	clusters_t tempclusters;
	tempclusters.N = (float*) malloc(sizeof(float) * original_num_clusters);
	tempclusters.pi = (float*) malloc(sizeof(float) * original_num_clusters);
	tempclusters.constant = (float*) malloc(sizeof(float) * original_num_clusters);
	tempclusters.avgvar = (float*) malloc(sizeof(float) * original_num_clusters);
	tempclusters.means = (float*) malloc(sizeof(float) * num_dimensions * original_num_clusters);
	tempclusters.R = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	tempclusters.Rinv = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	tempclusters.memberships = (float*) malloc(sizeof(float) * num_events * original_num_clusters);
	if (!tempclusters.means || !tempclusters.R || !tempclusters.Rinv || !tempclusters.memberships) {
		printf("ERROR: Could not allocate memory for clusters.\n");
		return 1;
	}

	// Declare another set of clusters for saving the results of the best configuration
	clusters_t saved_clusters;
	saved_clusters.N = (float*) malloc(sizeof(float) * original_num_clusters);
	saved_clusters.pi = (float*) malloc(sizeof(float) * original_num_clusters);
	saved_clusters.constant = (float*) malloc(sizeof(float) * original_num_clusters);
	saved_clusters.avgvar = (float*) malloc(sizeof(float) * original_num_clusters);
	saved_clusters.means = (float*) malloc(sizeof(float) * num_dimensions * original_num_clusters);
	saved_clusters.R = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	saved_clusters.Rinv = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions * original_num_clusters);
	saved_clusters.memberships = (float*) malloc(sizeof(float) * num_events * original_num_clusters);
	if (!saved_clusters.means || !saved_clusters.R || !saved_clusters.Rinv || !saved_clusters.memberships) {
		printf("ERROR: Could not allocate memory for clusters.\n");
		return 1;
	}

	// Used as a temporary cluster for combining clusters in "distance" computations
	clusters_t scratch_cluster;
	scratch_cluster.N = (float*) malloc(sizeof(float));
	scratch_cluster.pi = (float*) malloc(sizeof(float));
	scratch_cluster.constant = (float*) malloc(sizeof(float));
	scratch_cluster.avgvar = (float*) malloc(sizeof(float));
	scratch_cluster.means = (float*) malloc(sizeof(float) * num_dimensions);
	scratch_cluster.R = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions);
	scratch_cluster.Rinv = (float*) malloc(sizeof(float) * num_dimensions * num_dimensions);
	scratch_cluster.memberships = (float*) malloc(sizeof(float) * num_events);

	DEBUG("Finished allocating memory on host for clusters.\n");

	float min_rissanen, rissanen;

	//////////////// Initialization done, starting kernels ////////////////
	DEBUG("Invoking seed_clusters kernel.\n");
	fflush(stdout);

	// seed_clusters sets initial pi values,
	// finds the means / covariances and copies it to all the clusters
	// TODO: Does it make any sense to use multiple blocks for this?
	seed_start = clock();
	seed_clusters(fcs_data_by_event, &clusters, num_dimensions, original_num_clusters, num_events);

	DEBUG("Invoking constants kernel.\n");
	// Computes the R matrix inverses, and the gaussian constant
	//constants_kernel<<<original_num_clusters, num_threads>>>(d_clusters,original_num_clusters,num_dimensions);
	constants(&clusters, original_num_clusters, num_dimensions);
	constants_iterations++;
	seed_end = clock();
	seed_total = seed_end - seed_start;

	// Calculate an epsilon value
	//int ndata_points = num_events*num_dimensions;
	float epsilon = (1 + num_dimensions + 0.5 * (num_dimensions + 1) * num_dimensions) * log((float) num_events * num_dimensions) * 0.01;
	float likelihood, old_likelihood;
	int iters;

	epsilon = 1e-6;
	PRINT("Gaussian.cu: epsilon = %f\n", epsilon);

	// Variables for GMM reduce order
	float distance, min_distance = 0.0;
	int min_c1, min_c2;
	int ideal_num_clusters;

	for (int num_clusters = original_num_clusters; num_clusters >= stop_number; num_clusters--) {
		/*************** EM ALGORITHM *****************************/

		// do initial regrouping
		// Regrouping means calculate a cluster membership probability
		// for each event and each cluster. Each event is independent,
		// so the events are distributed to different blocks
		// (and hence different multiprocessors)
		DEBUG("Invoking regroup (E-step) kernel with %d blocks.\n",NUM_BLOCKS);
		regroup_start = clock();
		estep1(fcs_data_by_dimension, &clusters, num_dimensions, num_clusters, num_events, &likelihood);
		PRINT("\nNEW ESTEP2\n");
		estep2(fcs_data_by_dimension, &clusters, num_dimensions, num_clusters, num_events, &likelihood);
		//estep2b(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events,&likelihood);
		regroup_end = clock();
		regroup_total += regroup_end - regroup_start;
		regroup_iterations++;
		DEBUG("Regroup Kernel Iteration Time: %f\n\n",((double)(regroup_end-regroup_start))/CLOCKS_PER_SEC);

		DEBUG("Likelihood: %e\n",likelihood);

		float change = epsilon * 2;

//		old_likelihood = likelihood;

		PRINT("Performing EM algorithm on %d clusters.\n", num_clusters);
		iters = 0;
		// This is the iterative loop for the EM algorithm.
		// It re-estimates parameters, re-computes constants, and then regroups the events
		// These steps keep repeating until the change in likelihood is less than some epsilon
        while(iters < MIN_ITERS || (fabs(change) > epsilon && iters < MAX_ITERS && change > 0)) {
            old_likelihood = likelihood;
            // Save the current cluster as we loose membership in stoping iteration
            // Save the cluster configuration somewhere
			memcpy(tempclusters.N, clusters.N, sizeof(float) * num_clusters);
			memcpy(tempclusters.pi, clusters.pi, sizeof(float) * num_clusters);
			memcpy(tempclusters.constant, clusters.constant, sizeof(float) * num_clusters);
			memcpy(tempclusters.avgvar, clusters.avgvar, sizeof(float) * num_clusters);
			memcpy(tempclusters.means, clusters.means, sizeof(float) * num_dimensions * num_clusters);
			memcpy(tempclusters.R, clusters.R, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
			memcpy(tempclusters.Rinv, clusters.Rinv, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
			memcpy(tempclusters.memberships, clusters.memberships, sizeof(float) * num_events * num_clusters);


            DEBUG("Invoking reestimate_parameters (M-step) kernel.\n");
            params_start = clock();
            // This kernel computes a new N, pi isn't updated until compute_constants though
            mstep_n(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events);
            mstep_mean(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events);
            mstep_covar(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events);
            params_end = clock();
            params_total += params_end - params_start;
            params_iterations++;
            DEBUG("Model M-Step Iteration Time: %f\n\n",((double)(params_end-params_start))/CLOCKS_PER_SEC);
            //return 0; // RETURN FOR FASTER PROFILING

            DEBUG("Invoking constants kernel.\n");
            // Inverts the R matrices, computes the constant, normalizes cluster probabilities
            constants_start = clock();
            constants(&clusters,num_clusters,num_dimensions);
            constants_end = clock();
            constants_total += constants_end - constants_start;
            constants_iterations++;
            DEBUG("Constants Kernel Iteration Time: %f\n\n",((double)(constants_end-constants_start))/CLOCKS_PER_SEC);

            DEBUG("Invoking regroup (E-step) kernel with %d blocks.\n",NUM_BLOCKS);
            regroup_start = clock();
            // Compute new cluster membership probabilities for all the events
            estep1(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events,&likelihood);
            PRINT("\nESTEP2 IN WHILE\n");
            estep2(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events,&likelihood);
            //estep2b(fcs_data_by_dimension,&clusters,num_dimensions,num_clusters,num_events,&likelihood);
            regroup_end = clock();
            regroup_total += regroup_end - regroup_start;
            regroup_iterations++;
            DEBUG("E-step Iteration Time: %f\n\n",((double)(regroup_end-regroup_start))/CLOCKS_PER_SEC);

            change = likelihood - old_likelihood;
            DEBUG("likelihood = %f\n",likelihood);
            DEBUG("Change in likelihood: %f\n",change);

            iters++;

        }
        // Copy the cluster back
		memcpy(clusters.N, tempclusters.N, sizeof(float) * num_clusters);
		memcpy(clusters.pi, tempclusters.pi, sizeof(float) * num_clusters);
		memcpy(clusters.constant, tempclusters.constant, sizeof(float) * num_clusters);
		memcpy(clusters.avgvar, tempclusters.avgvar, sizeof(float) * num_clusters);
		memcpy(clusters.means, tempclusters.means, sizeof(float) * num_dimensions * num_clusters);
		memcpy(clusters.R, tempclusters.R, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
		memcpy(clusters.Rinv, tempclusters.Rinv, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
		memcpy(clusters.memberships, tempclusters.memberships, sizeof(float) * num_events * num_clusters);

		// Calculate Rissanen Score
		rissanen = -likelihood + 0.5 * (num_clusters * (1 + num_dimensions + 0.5 * (num_dimensions + 1) * num_dimensions) - 1) * logf((float) num_events * num_dimensions);
		PRINT("\nRissanen Score: %e\n", rissanen);

		// Save the cluster data the first time through, so we have a base rissanen score and result
		// Save the cluster data if the solution is better and the user didn't specify a desired number
		// If the num_clusters equals the desired number, stop
		if (num_clusters == original_num_clusters || (rissanen < min_rissanen && desired_num_clusters == 0) || (num_clusters == desired_num_clusters)) {
			min_rissanen = rissanen;
			ideal_num_clusters = num_clusters;
			// Save the cluster configuration somewhere
			memcpy(saved_clusters.N, clusters.N, sizeof(float) * num_clusters);
			memcpy(saved_clusters.pi, clusters.pi, sizeof(float) * num_clusters);
			memcpy(saved_clusters.constant, clusters.constant, sizeof(float) * num_clusters);
			memcpy(saved_clusters.avgvar, clusters.avgvar, sizeof(float) * num_clusters);
			memcpy(saved_clusters.means, clusters.means, sizeof(float) * num_dimensions * num_clusters);
			memcpy(saved_clusters.R, clusters.R, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
			memcpy(saved_clusters.Rinv, clusters.Rinv, sizeof(float) * num_dimensions * num_dimensions * num_clusters);
			memcpy(saved_clusters.memberships, clusters.memberships, sizeof(float) * num_events * num_clusters);

//			printf("\n\n COPIED CLUSTER TO SAVED_CLUSTERS\n\n");
		}

		/**************** Reduce GMM Order ********************/
		reduce_start = clock();
		// Don't want to reduce order on the last iteration
		if (num_clusters > stop_number) {
			// First eliminate any "empty" clusters
			for (int i = num_clusters - 1; i >= 0; i--) {
				if (clusters.N[i] < 1.0) {
					DEBUG("Cluster #%d has less than 1 data point in it.\n",i);
					for (int j = i; j < num_clusters - 1; j++) {
						copy_cluster(clusters, j, clusters, j + 1, num_dimensions);
					}
					num_clusters--;
				}
			}

			min_c1 = 0;
			min_c2 = 1;
			DEBUG("Number of non-empty clusters: %d\n",num_clusters);
			// For all combinations of subclasses...
			// If the number of clusters got really big might need to do a non-exhaustive search
			// Even with 100*99/2 combinations this doesn't seem to take too long
			for (int c1 = 0; c1 < num_clusters; c1++) {
				for (int c2 = c1 + 1; c2 < num_clusters; c2++) {
					// compute distance function between the 2 clusters
					distance = cluster_distance(clusters, c1, c2, scratch_cluster, num_dimensions);

					// Keep track of minimum distance
					if ((c1 == 0 && c2 == 1) || distance < min_distance) {
						min_distance = distance;
						min_c1 = c1;
						min_c2 = c2;
					}
				}
			}

			PRINT("\nMinimum distance between (%d,%d). Combining clusters\n", min_c1, min_c2);
			// Add the two clusters with min distance together
			//add_clusters(&(clusters[min_c1]),&(clusters[min_c2]),scratch_cluster,num_dimensions);
			add_clusters(clusters, min_c1, min_c2, scratch_cluster, num_dimensions);
			// Copy new combined cluster into the main group of clusters, compact them
			//copy_cluster(&(clusters[min_c1]),scratch_cluster,num_dimensions);
			copy_cluster(clusters, min_c1, scratch_cluster, 0, num_dimensions);
			for (int i = min_c2; i < num_clusters - 1; i++) {
				//printf("Copying cluster %d to cluster %d\n",i+1,i);
				//copy_cluster(&(clusters[i]),&(clusters[i+1]),num_dimensions);
				copy_cluster(clusters, i, clusters, i + 1, num_dimensions);
			}

		} // GMM reduction block
		reduce_end = clock();
		reduce_total += reduce_end - reduce_start;
		reduce_iterations++;
	} // outer loop from M to 1 clusters


	// Open up the output file for cluster summary

	PRINT("\nFinal rissanen Score was: %f, with %d clusters.\n", min_rissanen, ideal_num_clusters);

	char* result_suffix = ".results";
	char* summary_suffix = ".summary";
	int filenamesize1 = strlen("gmm") + strlen(result_suffix) + 1;
	int filenamesize2 = strlen("gmm") + strlen(summary_suffix) + 1;
	char* result_filename = (char*) malloc(filenamesize1);
	char* summary_filename = (char*) malloc(filenamesize2);
	strcpy(result_filename, "gmm");
	strcpy(summary_filename, "gmm");
	strcat(result_filename, result_suffix);
	strcat(summary_filename, summary_suffix);

	PRINT("Summary filename: %s\n", summary_filename);
	PRINT("Results filename: %s\n", result_filename);



	// Open up the output file for cluster summary
//	FILE	* outf = fopen(summary_filename, "w");
//	if (!outf) {
//		printf("ERROR: Unable to open file '%s' for writing.\n", "gmm");
//	}
//
//	// Print the clusters with the lowest rissanen score to the console and output file
//	for (int c = 0; c < ideal_num_clusters; c++) {
//		//if(saved_clusters.N[c] == 0.0) {
//		//    continue;
//		//}
//		if (ENABLE_PRINT) {
//			// Output the final cluster stats to the console
//			PRINT("Cluster #%d\n", c);
//			printCluster(saved_clusters, c, num_dimensions);
//			PRINT("\n\n");
//		}
//
//		if (ENABLE_OUTPUT) {
//			// Output the final cluster stats to the output file
//			fprintf(outf, "Cluster #%d\n", c);
//			writeCluster(outf, saved_clusters, c, num_dimensions);
//			fprintf(outf, "\n\n");
//		}
//	}

	// Print profiling information
//	printf("Program Component\tTotal\tIters\tTime Per Iteration\n");
//	printf("        Seed Kernel:\t%7.4f\t%d\t%7.4f\n", seed_total / (double) CLOCKS_PER_SEC, 1, (double) seed_total / (double) CLOCKS_PER_SEC);
//	printf("      E-step Kernel:\t%7.4f\t%d\t%7.4f\n", regroup_total / (double) CLOCKS_PER_SEC, regroup_iterations, (double) regroup_total / (double) CLOCKS_PER_SEC / (double) regroup_iterations);
//	printf("      M-step Kernel:\t%7.4f\t%d\t%7.4f\n", params_total / (double) CLOCKS_PER_SEC, params_iterations, (double) params_total / (double) CLOCKS_PER_SEC / (double) params_iterations);
//	printf("   Constants Kernel:\t%7.4f\t%d\t%7.4f\n", constants_total / (double) CLOCKS_PER_SEC, constants_iterations,
//			(double) constants_total / (double) CLOCKS_PER_SEC / (double) constants_iterations);
//	printf("GMM Order Reduction:\t%7.4f\t%d\t%7.4f\n", reduce_total / (double) CLOCKS_PER_SEC, reduce_iterations, (double) reduce_total / (double) CLOCKS_PER_SEC / (double) reduce_iterations);
//
//	// Write profiling info to summary file
//	fprintf(outf, "Program Component\tTotal\tIters\tTime Per Iteration\n");
//	fprintf(outf, "        Seed Kernel:\t%7.4f\t%d\t%7.4f\n", seed_total / (double) CLOCKS_PER_SEC, 1, (double) seed_total / (double) CLOCKS_PER_SEC);
//	fprintf(outf, "      E-step Kernel:\t%7.4f\t%d\t%7.4f\n", regroup_total / (double) CLOCKS_PER_SEC, regroup_iterations,
//			(double) regroup_total / (double) CLOCKS_PER_SEC / (double) regroup_iterations);
//	fprintf(outf, "      M-step Kernel:\t%7.4f\t%d\t%7.4f\n", params_total / (double) CLOCKS_PER_SEC, params_iterations, (double) params_total / (double) CLOCKS_PER_SEC / (double) params_iterations);
//	fprintf(outf, "   Constants Kernel:\t%7.4f\t%d\t%7.4f\n", constants_total / (double) CLOCKS_PER_SEC, constants_iterations,
//			(double) constants_total / (double) CLOCKS_PER_SEC / (double) constants_iterations);
//	fprintf(outf, "GMM Order Reduction:\t%7.4f\t%d\t%7.4f\n", reduce_total / (double) CLOCKS_PER_SEC, reduce_iterations, (double) reduce_total / (double) CLOCKS_PER_SEC / (double) reduce_iterations);
//	fclose(outf);

	// Open another output file for the event level clustering results
//	FILE* fresults = fopen(result_filename, "w");
//
//	if (ENABLE_OUTPUT) {
//		for (int i = 0; i < num_events; i++) {
//			for (int d = 0; d < num_dimensions - 1; d++) {
//				fprintf(fresults, "%f,", fcs_data_by_event[i * num_dimensions + d]);
//			}
//			fprintf(fresults, "%f", fcs_data_by_event[i * num_dimensions + num_dimensions - 1]);
//			fprintf(fresults, "\t");
//			for (int c = 0; c < ideal_num_clusters - 1; c++) {
//				fprintf(fresults, "%f,", saved_clusters.memberships[c * num_events + i]);
//			}
//			fprintf(fresults, "%f", saved_clusters.memberships[(ideal_num_clusters - 1) * num_events + i]);
//			fprintf(fresults, "\n");
//		}
//	}
//	fclose(fresults);

	float current_likelihood = 0.0;
	for (int n = 0; n < num_events; n++) {
		clusterIdx[n] = 0;
		current_likelihood = 0.0;
		for (int m = 0; m < desired_num_clusters; m++) {

			current_likelihood = saved_clusters.memberships[m * num_events + n] > current_likelihood? saved_clusters.memberships[m * num_events + n]:current_likelihood;
			clusterIdx[n] = saved_clusters.memberships[m * num_events + n] >= current_likelihood? m:clusterIdx[n];

//			printf("Membership of event %d in cluster %d: %.3f\n", n, m, saved_clusters.memberships[m * num_events + n]);
		}
	}

	// cleanup host memory
//	free(fcs_data_by_event);
	free(fcs_data_by_dimension);
	free(clusters.N);
	free(clusters.pi);
	free(clusters.constant);
	free(clusters.avgvar);
	free(clusters.means);
	free(clusters.R);
	free(clusters.Rinv);
	free(clusters.memberships);

	free(saved_clusters.N);
	free(saved_clusters.pi);
	free(saved_clusters.constant);
	free(saved_clusters.avgvar);
	free(saved_clusters.means);
	free(saved_clusters.R);
	free(saved_clusters.Rinv);
	free(saved_clusters.memberships);

	free(scratch_cluster.N);
	free(scratch_cluster.pi);
	free(scratch_cluster.constant);
	free(scratch_cluster.avgvar);
	free(scratch_cluster.means);
	free(scratch_cluster.R);
	free(scratch_cluster.Rinv);
	free(scratch_cluster.memberships);

}
