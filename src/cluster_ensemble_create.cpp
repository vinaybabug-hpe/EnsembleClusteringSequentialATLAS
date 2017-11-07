/*
 ============================================================================
 Name        : SpectralClustering.cu
 Author      : Vinay B Gavirangaswamy
 Version     :
 Copyright   : Modified Source of Spectral Clustering CUDA Implementation by Jin Yu from UMD
 Description : CUDA compute reciprocals
 ============================================================================
 */

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
#include <unordered_set>


#include "mat.h"

#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "common/cluster_data_outliers.h"
#include "common/cluster_util_bootpartition2partition.h"
#include "common/cluster_ensemble2cam.h"
#include "common/cluster_util_indices2centers.h"
#include "common/indices_count.h"
#include "common/cluster_util_ssw.h"




using namespace std;
////////////////////////////////////////////////////////////////////////////////
// Main entry point.
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // NON MEX related variables
	char* en_name = "Ensemble Clustering Results using CUDA";
	char* en_desc;
	char *output_file;
	int ifield, nfields, field_num;
    int nStructElems = 0;
    int VERBOSE_Q, INCLUDE_REPsQ, INCLUDE_CENTERsQ;
    int nBootstrap = 0;
    int nReps = 0;
    char *name;
    double *kList;
    char **clustList;
    char **modelList;
    int nKList, nModelList;
    int m_x, n_x;
    int *bootIdxs;
    double *x_actual_data;
    float *x_transpose;
    int count_x_transpose;
	char *distmetric;
	char *centerfun;
	char *linkcode;
	char *outlierMetric;
	float outlierCutoff;

	int num_total_partitions = 0;
	int *data_partition_indices;


    // MEX related variables
    MATFile *mfPtr, *mfPtr2; /* MAT-file pointer */
    mxArray *aPST, *aX, *aBootIdxs;  /* mxArray pointer */
    mxArray *tmp, *pst_ensemble, *pst_extract;
	const mxArray *cell_element_ptr;
	mxArray *output_cell_array_ptr;

	// CUDA related variables


    if (argc < 4){
		cout << "Not enough input arguments!" << endl;
		cout << "The input format is: " << endl;
		cout << "1. PST .mat Filename" << endl;
		cout << "2. DATA .mat Filename" << endl;
		cout << "3. Ensemble Solution Description" << endl;
		cout << "4. Output Filename (.mat)" << endl;
		exit(1);
	}


    // TODO: Debug code remove later!
//    for(int count = 0; count < argc; count++){
//    	cout << "argv["<<count<<"]: " << argv[count] << endl;
//    }

    mfPtr = matOpen(argv[1], "r");
    if (mfPtr == NULL) {
        printf("Error opening file %s\n", argv[1]);
        return(EXIT_FAILURE);
    }

    aPST = matGetVariable(mfPtr, MEX_STRUCT_PST_NME);
    if (aPST == NULL) {
        printf("mxArray not found: %s\n", MEX_STRUCT_PST_NME);
        return(EXIT_FAILURE);
    }

    mfPtr2 = matOpen(argv[2], "r");
    if (mfPtr2 == NULL) {
    	printf("Error opening file %s\n", argv[2]);
    	return(EXIT_FAILURE);
    }

    aX = matGetVariable(mfPtr2, MEX_STRUCT_DATA_NME);
    if (aX == NULL) {
    	printf("mxArray not found: %s\n", MEX_STRUCT_PST_NME);
    	return(EXIT_FAILURE);
    }

    en_desc = argv[3];

    output_file = argv[4];

	nfields = mxGetNumberOfFields(aPST);
	nStructElems = mxGetNumberOfElements(aPST);
	/* check proper input and output */
	if (nStructElems != 1){
		printf("ENCLUST CUDA: Multiple PST structures \n Only one input is required.");
	}

	tmp = mxGetField(aPST, 0, MEX_STRUCT_PST_FIELD_VERBOSE);
	VERBOSE_Q = mxGetScalar(tmp);
	field_num = mxGetFieldNumber(aPST, MEX_STRUCT_PST_FIELD_NAME);
	tmp = mxGetFieldByNumber(aPST, 0, field_num);
	name = mxArrayToString(tmp);

	outlierCutoff = mxGetScalar(
				mxGetField(aPST, 0, MEX_STRUCT_OUTLIERCUTOFF_NME));

	outlierMetric = mxArrayToString(mxGetField(aPST, 0, MEX_STRUCT_OUTLIERMETRIC_NME));

//	cout<<"\nNO PROBLEM TILL HERE"<<endl;

	pst_ensemble = mxGetField(aPST, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE);

	/**
	 * Extract data from ENSEMBLE struct
	 * ENSEMBLE is inside PST
	 */
	INCLUDE_REPsQ = mxGetScalar(
			mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_INCLUDEREPSQ));
	INCLUDE_CENTERsQ = mxGetScalar(
			mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_INCLUDECENTERSQ));
	nBootstrap = mxGetScalar(mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_NBOOTSTRAPS));

	if (INCLUDE_REPsQ == 1) {
		nReps = mxGetScalar(mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_NREPS));
	} else {
		nReps = 1;
	}

	tmp = mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_KLIST);
	nKList = mxGetN(tmp);
	kList = mxGetPr(tmp);

	/*
	 * Get list of model and cluster to run.
	 */

	tmp = mxGetField(pst_ensemble, 0, MEX_STRUCT_PST_FIELD_ENSEMBLE_MODELLIST);
	nModelList = mxGetNumberOfElements(tmp);

	modelList = (char**) malloc(nModelList * sizeof(char*));

	for (int j = 0; j < nModelList; j++) {
		cell_element_ptr = mxGetCell(tmp, j);
		modelList[j] = mxArrayToString(cell_element_ptr);
		/*fprintf(ofp, "%s\n", modelList[j]);*/
	}


	m_x = mxGetM(aX);
	n_x = mxGetN(aX);
	x_actual_data = mxGetPr(aX);

    // size of memory required to store the matrix
    size_t mem_size = static_cast<size_t>(sizeof(float) * m_x* n_x);

	/*
	 * Transpose X as data stored in single dimensional array
	 * is row major not column
	 */

	x_transpose = (float*)malloc(mem_size);

	// Use x_transpose as temp buffer to convert double* to float*
//	std::copy(x_actual_data, x_actual_data + m_x* n_x, x_transpose);


	for (int i = 0; i < m_x; i++) {
		for (int j = 0; j < n_x; j++) {
			x_transpose[i*n_x+j] = x_actual_data[i + j * m_x];

		}
	}

	// Print PST to see
	printf("\n***************************************************************************\n");
	printf("\t\t\t\tPST\n");
	printf("***************************************************************************\n\n");

	printf("\n\tname: %s", name);
	printf("\n\tkList:\n");
	for (int j = 0; j < nKList; j++) {
		printf("\t%d\n", (int) kList[j]);
	}
	printf("\n\tModel List:\n");
	for (int j = 0; j < nModelList; j++) {
		printf("\t%s\n", modelList[j]);
	}
	printf("\n\tnReps : %d\n", nReps);
	printf("\n\toutlierCutoff: %f", outlierCutoff);
	printf("\n\toutlierMetric: %s", outlierMetric);

	num_total_partitions = nBootstrap * nModelList * nKList * nReps;

    printf("\n\tStarting Ensemble Clustering on CUDA for total of %d partitions...\n", num_total_partitions);
    double createen_time_spent = 0;


    /**
     * check to make sure max k is not too big for data set
     */

    double max_value = 0;

    for(int count=0;count < nKList; count++){
    	max_value = max_value < kList[count]? kList[count]: max_value;
    }


    if( max_value > m_x/3){
    	cout<<"\tThe largest value of k in the kList is too big for the data set. The largest k should be less than 1/3 of nSamples.\n"<<endl;
    	return (EXIT_FAILURE);
    }

    /**
     * Pre-Process data
     */
    printf("\n\tStarting pre-processing...\n");
    bool *outlierMask;
    int* outlierIdxs_data;
    int* outlierIdxs_size;
    float* outlierZ_data;
    int* outlierZ_size;
    outlierCutoff = 1;
    int nrows_no_outliers = 0;

    clock_t begin = clock();

    float *x_outlier_removed = cluster_data_outliers(m_x,
    						 n_x,
    						 x_transpose,
    						 outlierMetric,
    						 FLT_MAX/*outlierCutoff*/,
    						 &nrows_no_outliers);

    int nrows = nrows_no_outliers;
    int ncols = n_x;

    printf("\n\tDone pre-processing...\n");
    printf("\n\t%d objects were removed from original dataset...\n",m_x - nrows);

    clock_t end = clock();
    double preprocessing_time_spent = (double)(end - begin)/(double)CLOCKS_PER_SEC;

    printf("\n***************************************************************************\n");
    printf("Time spent in pre-processing for %d objects is %f\n", m_x, preprocessing_time_spent);
    printf("***************************************************************************\n\n");
    printf("\n\tGenerating boot indices...\n");

    mxArray *bootIdxsm = NULL;
    bootIdxsm = mxCreateNumericMatrix(nrows, nBootstrap, mxINT32_CLASS, mxREAL);
    bootIdxs = (int*)mxGetPr(bootIdxsm);
    int *bootIdxs2 = (int*)malloc(nrows * nBootstrap * sizeof(int));

    for(int row=0; row < nrows; row++){
    	for(int col=0; col < nBootstrap; col++){
    		bootIdxs2[row * nBootstrap + col] = (rand()%(nrows-0))+0;;
    		bootIdxs[row + col * nrows] = bootIdxs2[row * nBootstrap + col];
    	}
    }

	/*Do the clustering*/
	printf("\n***************************************************************************\n");
	printf("\t\t\tGenerating cluster ensemble\n");
	printf("***************************************************************************\n\n");

	data_partition_indices = (int*) malloc(num_total_partitions * nrows * sizeof(int));

	printf("\n\tRunning Ensemble with %d Bootstrap...\n", nBootstrap);

	//Setup the output .mat file structure.
	int total_models = 0;
	int index, nsubs=5;
	const char *field_names_e[] = {"name", "description","parameters","inputData","bootIndices", "partitions", "ensolution"};
	int num_fields_e = 7; // Basically length of [] field_names_e
	const char *field_names_partitions[] = {"model", "boot","k","rep","indices"};
	int num_fields_partitions = 5; // Basically length of [] field_names_partitions
	size_t subs[2];
	mxArray *partitions_array_ptr;
//	output_cell_array_ptr = mxCreateCellMatrix(nBootstrap*nModelList*nKList*nReps , nsubs);

	mwSize dims_partitions[2] = {1, nBootstrap*nModelList*nKList*nReps };
	mwSize dims_e[2] = {1, 1 };
	/* Create a 1-by-n array of structs. */
	partitions_array_ptr = mxCreateStructArray(2, dims_partitions, num_fields_partitions, field_names_partitions);

	output_cell_array_ptr = mxCreateStructArray(2, dims_e, num_fields_e, field_names_e);

	// Set partitions struct in output
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[0]), mxCreateString(en_name));
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[1]), mxCreateString(en_desc));
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[2]), aPST);
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[4]), bootIdxsm);

	mxArray *Xm = NULL;

	// copy data to E (output struct).
	Xm = mxCreateDoubleMatrix (nrows, ncols, mxREAL);
	double *X_p;
	X_p = mxGetPr(Xm);
	for (int i = 0; i < nrows; i++){
#pragma unroll
		for (int j = 0; j < ncols; j++){
			X_p[i + j * nrows] = x_outlier_removed[i * ncols + j];
		}
	}
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[3]), Xm);


	/*
	 * Below loops go through
	 * 1. nBootstrap data samples
	 * 2. nModelList e.g. 'kmeanseuc','kmeanscit','kmeanscor','kmeanscos',...
	 * 					  'spectraleuc','spectralcit','spectralcor','spectralcos',...
	 * 					  'aggeucwar','aggeucavg','aggeuccom',...
	 * 					  'aggcitavg','aggcitcom',...
	 * 					  'aggcoravg','aggcorcom',...
	 * 					  'aggcosavg','aggcoscom'
	 * 3. nKList: Number of clusters to be found in data for combinations of bootstrap data and model  e.g. 1, 2, 3..5,10
	 * 4. nReps: Repeat the model (combination of bootstrap data sample, model, #cluster) to avoid bias.
	 */
	for (int c_nBootstrap = 0; c_nBootstrap < nBootstrap; c_nBootstrap++) {

		/*
		 * Variables to keep track of memory used on device.
		 * Depending on the current memory usage we need to restrict
		 * number of threads.
		 */

		float* weight_host;
		float* data1d_host;
		//int *tri_idxs_host;
		size_t memDataSz = 0, memWtSz = 0/*, memIdxSz=0*/;
		memDataSz = nrows * ncols * sizeof(float);
		data1d_host = (float*) malloc(nrows*ncols*sizeof(float));
		assert(data1d_host!=NULL);

		/**
		 * Generate bootstrap data and write it
		 * to a temporary file to use by clustering
		 * algorithms.
		 * Create bootstrap data sets...
		 */
		int *b = (int*) malloc(nrows * sizeof(int));

		for (int j = 0; j < nrows; j++) {
			memcpy(	&data1d_host[j * ncols], &x_outlier_removed[bootIdxs2[j * nBootstrap + c_nBootstrap] * ncols], ncols * sizeof(float));
			b[j] = bootIdxs2[j * nBootstrap + c_nBootstrap];
		}


		//TODO: Debug code remove later!
//		printf("\n-------------------------------------------------------------------------------------- \n\n");
//		for (int m = 0; m < nrows; m++) {
//			for (int n = 0; n < ncols; n++) {
//
////				data1d_managed[n + m * ncols] = x_outlier_removed[bootIdxs2[m * nBootstrap + c_nBootstrap] * ncols+n];
////				printf("%.3f ", data1d_managed[n + m * ncols]);
//				printf("%.3f ", x_outlier_removed[(bootIdxs2[m * nBootstrap + c_nBootstrap] * ncols)+n]);
////				printf("%d ", bootIdxs2[m * nBootstrap + c_nBootstrap]>=nrows-1?1:0);
//			}
//			b[m] = bootIdxs2[m * nBootstrap + c_nBootstrap];
//			printf(" \n");
//		}

//		printf("\n NO PROBLEM TILL HERE \n");


		memWtSz = ncols * sizeof(float);

		weight_host = (float*) malloc(memWtSz);
		assert(weight_host != NULL);

//		memIdxSz = TRI_COUNT(nrows)*sizeof(int);
//
//		tri_idxs_host = (int*) malloc(memIdxSz);

		for (int i = 0; i < ncols; i++)
			weight_host[i] = 1.0;

//		int jCount = 0;
//		for (int i = 1; i < nrows; i++){
//			for (int j = 0; j < i; j++){
//				tri_idxs_host[jCount] = TRI_COUNT(i)+j;
//				jCount++;
//			}
//		}


		for (int c_nModels = 0; c_nModels < nModelList; c_nModels++) {
			for (int c_nKList = 0; c_nKList < nKList; c_nKList++) {
				for (int c_nReps = 1; c_nReps < nReps+1; c_nReps++) {

					mxArray *idxs_ptr = mxCreateNumericMatrix(1, nrows, mxINT32_CLASS, mxREAL);
					int *mxIdxs = (int*)mxGetData(idxs_ptr);

					int *idxs = (int*)malloc(nrows * sizeof(int));
					assert(idxs != NULL);
					assert(mxIdxs != NULL);

//					printf("%s %d %d %d \n", modelList[c_nModels], (int)kList[c_nKList], c_nBootstrap, c_nReps);

					distmetric = (char*) malloc(MODEL_STR_LEN * sizeof(char));
					centerfun = (char*) malloc(MODEL_STR_LEN * sizeof(char));
					linkcode = (char*) malloc(MODEL_STR_LEN * sizeof(char));

					clock_t begin = clock();

				     // take measurements for loop over kernel launches
//				     checkCudaErrors(cudaEventRecord(start, 0));

					// DO SPECTRAL
					if(getModelType(SPECTRAL_SHRT, modelList[c_nModels])){

						getDistMtrcnCntrFunBySpectral(SPECTRAL_SHRT, modelList[c_nModels], distmetric, centerfun);
//						cout<<distmetric << " : "<< centerfun<<endl;
						mex2seq_spectral_adapter((int)kList[c_nKList], nrows, ncols, data1d_host, NULL, distmetric, idxs);
					}
					// DO KMEANS
					else if(getModelType(KMEANS_SHRT, modelList[c_nModels]))
					{
						float threshold = 0.001;
						int loop_iterations;
						int** mask;
						float **data2d;
						malloc2D(mask, nrows, ncols, int);

						malloc2D(data2d, nrows, ncols, float);
						memcpy(data2d[0], data1d_host, nrows * ncols * sizeof(float));

						for (int row = 0; row < nrows; row++){
							for (int col = 0; col < ncols; col++){
								mask[row][col] = data2d[row][col] == 0? 0 : 1;
							}
						}

						getDistMtrcnCntrFunByKmeans(KMEANS_SHRT, modelList[c_nModels], distmetric, centerfun);
//						cout<<distmetric << " : "<< centerfun<<endl;
						mex2seq_kmeans_adapter((int)kList[c_nKList], nrows, ncols, data2d, mask, weight_host,NULL, distmetric, threshold, idxs);

						free(data2d[0]);
						free(data2d);
						free(mask[0]);
						free(mask);

					}
					// DO KMEDOID
					else if(getModelType(KMEDOIDS_SHRT, modelList[c_nModels]))
					{
						float threshold = 0.001;
						int loop_iterations;
						int** mask;
						float **data2d;
						malloc2D(mask, nrows, ncols, int);

						malloc2D(data2d, nrows, ncols, float);
						memcpy(data2d[0], data1d_host, nrows * ncols * sizeof(float));

						for (int row = 0; row < nrows; row++){
							for (int col = 0; col < ncols; col++){
								mask[row][col] = data2d[row][col] == 0? 0 : 1;
							}
						}

						getDistMtrcnCntrFunByKmedoid(KMEDOIDS_SHRT, modelList[c_nModels], distmetric, centerfun);
//						cout<<distmetric << " : "<< centerfun<<endl;
						mex2seq_kmedians_adapter((int)kList[c_nKList], nrows, ncols, data2d, mask, weight_host,NULL, distmetric, threshold, idxs);

						free(data2d[0]);
						free(data2d);
						free(mask[0]);
						free(mask);
					}
					// DO KMEDOID
					else if(getModelType(GMM_SHRT, modelList[c_nModels]))
					{
						mex2seq_gmm_adapter((int)kList[c_nKList], nrows, ncols, data1d_host, NULL, NULL, idxs);
					}
					// DO KMEDOID
					else if(getModelType(AGGLO_SHRT, modelList[c_nModels]))
					{
						float **data2d;
						malloc2D(data2d, nrows, ncols, float);
						memcpy(data2d[0], data1d_host, nrows * ncols * sizeof(float));


						getDistMtrcnCntrFunByAgglo(AGGLO_SHRT, modelList[c_nModels], distmetric, centerfun, linkcode);
						mex2seq_agglomerative_adapter((int)kList[c_nKList], nrows, ncols, data2d, linkcode, distmetric, idxs);

						free(data2d[0]);
						free(data2d);
					}

//			        checkCudaErrors(cudaEventRecord(stop, 0));
//			        checkCudaErrors(cudaEventSynchronize(stop));
//					float kernelTime;
//					checkCudaErrors(cudaEventElapsedTime(&kernelTime, start, stop));

					clock_t end = clock();
					double locat_time_spent = (double)(end - begin)/(double)CLOCKS_PER_SEC;
					createen_time_spent += locat_time_spent;
//					time_spent += kernelTime;
#ifdef MEASURE_TIME
					printf("\n\tTime spent in %+15s (model), %4d (bootstrap), %3d (k), %2d-Rep  is: %3.4f", modelList[c_nModels], c_nBootstrap, (int)kList[c_nKList], c_nReps, locat_time_spent);
#endif

//					for(int i = 0; i < nrows; i++){
//						printf("%d ", idxs[i]);
//					}
//					printf("\n");

					/**
					 * cluster indices for boot sample need to be converted to
					 * cluster indices for original data samples
					 * NB:  do this BEFORE computing the centers...so centers
					 * computed only based on original data, not boostrap data
					 */
					cluster_util_bootpartition2partition(nrows, b, idxs, mxIdxs);
					/**
					 * Copy data partition indices into a bigger storage matrix
					 * that will be used later to extract ensemble solution later.
					 */
					memcpy(&data_partition_indices[total_models * nrows], mxIdxs, nrows * sizeof(int));

					/* Place the Model string array into cell element (total_models,0). */
//					subs[0]=total_models; subs[1]=0;
//					index = mxCalcSingleSubscript(output_cell_array_ptr, nsubs, subs);
//					mxSetCell( output_cell_array_ptr, index, mxCreateString(modelList[c_nModels]) );

					mxSetFieldByNumber(partitions_array_ptr,total_models, mxGetFieldNumber(partitions_array_ptr,field_names_partitions[0]),mxCreateString(modelList[c_nModels]));


					/* Place the Bootstrap number into cell element (total_models,1). */
//					subs[0]=total_models; subs[1]=1;
//					index = mxCalcSingleSubscript(output_cell_array_ptr, nsubs, subs);
//					mxSetCell( output_cell_array_ptr, index, mxCreateDoubleScalar(c_nBootstrap) );

					mxSetFieldByNumber(partitions_array_ptr,total_models, mxGetFieldNumber(partitions_array_ptr,field_names_partitions[1]),mxCreateDoubleScalar(c_nBootstrap));
					/* Place the #cluster (k) into cell element (total_models,2). */
//					subs[0]=total_models; subs[1]=2;
//					index = mxCalcSingleSubscript(output_cell_array_ptr, nsubs, subs);
//					mxSetCell( output_cell_array_ptr, index, mxCreateDoubleScalar(kList[c_nKList]) );

					mxSetFieldByNumber(partitions_array_ptr,total_models, mxGetFieldNumber(partitions_array_ptr,field_names_partitions[2]), mxCreateDoubleScalar(kList[c_nKList]));

					/* Place the #repition into cell element (total_models,3). */
//					subs[0]=total_models; subs[1]=3;
//					index = mxCalcSingleSubscript(output_cell_array_ptr, nsubs, subs);
//					mxSetCell( output_cell_array_ptr, index, mxCreateDoubleScalar(c_nReps) );

					mxSetFieldByNumber(partitions_array_ptr,total_models, mxGetFieldNumber(partitions_array_ptr,field_names_partitions[3]), mxCreateDoubleScalar(c_nReps));
					/* Place the idx into cell element (total_models,4). */
//					subs[0]=total_models; subs[1]=4;
//					index = mxCalcSingleSubscript(output_cell_array_ptr, nsubs, subs);
//					mxSetCell( output_cell_array_ptr, index, idxs_ptr );

					mxSetFieldByNumber(partitions_array_ptr,total_models, mxGetFieldNumber(partitions_array_ptr,field_names_partitions[4]), idxs_ptr);
					//free(idxs);

					total_models++;

					free(distmetric);
					free(centerfun);
					free(linkcode);
					free(idxs);
				}
			}
		}


//		free(tri_idxs_host);

		free(weight_host);
		free(data1d_host);
		free(b);
	}

	// Set partitions struct in output
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[5]), partitions_array_ptr);

	// Print a message describing what the sample does.
	printf("\n\n***************************************************************************\n");
	printf("Time spent to Create Ensemble of Clusters for %d objects is %f\n", m_x, createen_time_spent);
	printf("***************************************************************************\n\n");

	/*----------------------------------START OF ENSEMBLE EXTRACTIN CODE----------------------*/

	printf("\n***************************************************************************\n");
	printf("\t\t\tExtracting Solution from ensemble\n");
	printf("***************************************************************************\n\n");

	pst_extract = mxGetField(aPST, 0, MEX_STRUCT_PST_FIELD_EXTRACT);

	tmp = mxGetField(pst_extract, 0, MEX_STRUCT_PST_FIELD_EXTRACT_KLIST);
	int extract_nKList = mxGetN(tmp);
	double *extract_kList = mxGetPr(tmp);

	const char *field_names_ensolution[] = {"dissimMatrix","coassocSoln"};
	int num_fields_ensolution = 2; // Basically length of [] field_names_partitions
	mxArray *ensolution_array_ptr;
	mwSize dims_ensolution[2] = {1, 1};

	const char *field_names_coassocSoln[] = {"k", "centers","clusterN","indices","clusterSSW"};
	int num_fields_coassocSoln = 5; // Basically length of [] field_names_partitions

	mxArray *coassocSoln_array_ptr;

	mwSize dims_coassocSoln[2] = {1, extract_nKList};

	/* Create a 1-by-n array of structs. */
	coassocSoln_array_ptr = mxCreateStructArray(2, dims_coassocSoln, num_fields_coassocSoln, field_names_coassocSoln);

	/* Create a 1-by-n array of structs. */
	ensolution_array_ptr = mxCreateStructArray(2, dims_ensolution, num_fields_ensolution, field_names_ensolution);

	mxArray *dissimMatrix_ptr = mxCreateNumericMatrix(nrows, nrows, mxDOUBLE_CLASS, mxREAL);
	double *mxdissimMatrix = mxGetData(dissimMatrix_ptr);

	float *A, *N;
	A = (float*) malloc(nrows * nrows * sizeof(float));
	N = (float*) malloc(nrows * nrows * sizeof(float));

	printf("\n\tConstructing co-association matrix from %d partitions...", num_total_partitions);
	cluster_ensemble2cam(data_partition_indices,  A, N, nrows, num_total_partitions);

	// TODO: SAVE A TO OUTPUT STRUCTURE

	for(int i =0; i < nrows; i++){
		for(int j=0; j< nrows; j++){
			mxdissimMatrix[i + j * nrows] = A[i * nrows + j];
//			printf("%f ", A[i * nrows + j]);
		}
//		printf("\n");
	}

	// save dissimMatrix
	mxSetFieldByNumber(ensolution_array_ptr,0, mxGetFieldNumber(ensolution_array_ptr,field_names_ensolution[0]), dissimMatrix_ptr);

	 /**
	  * use spectral clustering to get solutions for each k in kList
	  */
	   printf("\n\tClustering the co-association matrix...");
	   begin = clock();
	   for (int c_nKList = 0; c_nKList < extract_nKList; c_nKList++) {

		   mxArray *idxs_ptr = mxCreateNumericMatrix(1, nrows, mxINT32_CLASS, mxREAL);
		   int *mxIdxs = (int*)mxGetData(idxs_ptr);

		   // use spectral clustering to get solutions for each k in kList

		   mex2seq_spectral_adapter_full_distmat((int)extract_kList[c_nKList], nrows, nrows, A, NULL, "euc", mxIdxs);

//		   mex2seq_spectral_adapter((int)extract_kList[c_nKList], nrows, nrows, A, NULL, "euc", mxIdxs);

		   /**
		    * Computes cluster centers based on indices
		    */
		   int nClusters = unique_length(mxIdxs, nrows);
		   float *M = (float*) malloc(nClusters * ncols * sizeof(float));

		   cluster_util_indices2centers(x_outlier_removed, nrows, ncols, mxIdxs, "mean", nClusters, M);

		   mxArray *M_ptr = mxCreateNumericMatrix(nClusters, ncols, mxDOUBLE_CLASS, mxREAL);
		   double *M_data = mxGetData(M_ptr);

		   for(int i =0; i < nClusters; i++){
		   		for(int j=0; j< ncols; j++){
		   			M_data[i + j * nClusters] = M[i * ncols + j];

		   		}
		   	}
		   free(M);

		   /**
		    * Count how many objects were assigned to each labels
		    */
		   mxArray *valueList_ptr = mxCreateNumericMatrix(1, nClusters, mxINT32_CLASS, mxREAL);
		   int *valueList = (int*)mxGetData(valueList_ptr);
		   indices_count(mxIdxs, nrows, valueList, nClusters);

		   /**
		    *  Compute the sum squared error within, for each cluster
		    */
		   mxArray *ssw_ptr = mxCreateNumericMatrix(1, nClusters, mxDOUBLE_CLASS, mxREAL);
		   double *ssw_data = mxGetData(ssw_ptr);

		   cluster_util_ssw(x_outlier_removed, nrows, ncols, mxIdxs, ssw_data, nClusters);


		   /* Place the #cluster (k) into cell element (c_nKList,0). */
   		   mxSetFieldByNumber(coassocSoln_array_ptr,c_nKList, mxGetFieldNumber(coassocSoln_array_ptr,field_names_coassocSoln[0]), mxCreateDoubleScalar(extract_kList[c_nKList]));

		   /* Place the centers matrix into cell element (c_nKList,1). */
		   mxSetFieldByNumber(coassocSoln_array_ptr,c_nKList, mxGetFieldNumber(coassocSoln_array_ptr,field_names_coassocSoln[1]), M_ptr);

		   /* Place the number objects were assigned to each labels into cell element (c_nKList,2). */
		   mxSetFieldByNumber(coassocSoln_array_ptr,c_nKList, mxGetFieldNumber(coassocSoln_array_ptr,field_names_coassocSoln[2]), valueList_ptr);

		   /* Place the idx into cell element (c_nKList,3). */
		   mxSetFieldByNumber(coassocSoln_array_ptr,c_nKList, mxGetFieldNumber(coassocSoln_array_ptr,field_names_coassocSoln[3]), idxs_ptr);

		   /* Place the ssw into cell element (c_nKList,4). */
		   mxSetFieldByNumber(coassocSoln_array_ptr,c_nKList, mxGetFieldNumber(coassocSoln_array_ptr,field_names_coassocSoln[4]), ssw_ptr);

	   }

	   end = clock();
	   double extract_time_spent = (double)(end - begin)/(double)CLOCKS_PER_SEC;

		printf("\n\n***************************************************************************\n");
		printf("Time spent in Extracting Solution from ensemble for %d objects is %f\n", m_x, extract_time_spent);
		printf("***************************************************************************\n\n");

	   // Set coassocSoln struct in ensolution struct
	   	mxSetFieldByNumber(ensolution_array_ptr,0, mxGetFieldNumber(ensolution_array_ptr,field_names_ensolution[1]), coassocSoln_array_ptr);
	/*-----------------------------------END OF ENSEMBLE EXTRACTIN CODE-----------------------*/



	// Set ensolution struct in output
	mxSetFieldByNumber(output_cell_array_ptr,0, mxGetFieldNumber(output_cell_array_ptr,field_names_e[6]), ensolution_array_ptr);

	// TODO: might want to create a seperate method for output writing (but I am lazy for now)
	remove(output_file);
	MATFile *pmat = matOpen(output_file, "w");

	if (pmat == NULL) {
		printf("Error reopening file %s\n", output_file);
		return(EXIT_FAILURE);
	}
	int status = matPutVariable(pmat, "E", output_cell_array_ptr);
	if (status != 0) {
		printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
		return(EXIT_FAILURE);
	}

	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n",output_file);
		return(EXIT_FAILURE);
	}

    mxDestroyArray(aPST);
    mxDestroyArray(aX);


    if (matClose(mfPtr) != 0) {
        printf("Error closing file %s\n", argv[1]);
        return(EXIT_FAILURE);
    }
    if (matClose(mfPtr2) != 0) {
    	printf("Error closing file %s\n", argv[2]);
    	return(EXIT_FAILURE);
    }


    cout<<endl<<"\n\t...ensemble generation complete :)\n"<<endl;
    printf("\n\n***************************************************************************\n");
    printf("Time spent in Ensemble Clustering for %d objects is %f\n", m_x, preprocessing_time_spent + createen_time_spent + extract_time_spent);
    printf("***************************************************************************\n\n");
//    CUDA_CHECK_RETURN(cudaDeviceReset());
    free(x_transpose);
	free(modelList);
    free(data_partition_indices);
    free(bootIdxs2);
    return(EXIT_SUCCESS);
}

void testSpectral(){
    // Declare Variables
	string line;
	int testcase = 1;

	int nrows = 0, ncols = 0;
	int    *idxs;    /* [numObjs] */
	float *data1d;
	int numClusters = 0;

//	for (int row = 0; row < nrows; row++){
//		data1d[row] = (double*)malloc(nrows* sizeof(double));
//		for (int col = 0; col < ncols; col++){
//
//			data1d[row][col] = ((float)rand() / (RAND_MAX)) + 1;
//
////			cout << " "<<data1d[row][col];
//
//		}
////		cout<<endl;
//	}

	switch(testcase){
	case 1:{
		nrows = 10;
		ncols = 3;
		numClusters = 3;

		ifstream in("dataset/data10x3.txt");

		/* start the timer for the core computation -----------------------------*/
		/* membership: the cluster id for each data object */
		idxs = (int*)malloc(nrows * sizeof(int));
		assert(idxs != NULL);
		data1d = (float*)malloc(nrows*ncols * sizeof(float));
		assert(data1d != NULL);
		if (!in) {
			cout << "Cannot open file.\n";
			exit(EXIT_SUCCESS);
		}

		for(int i = 0; i < nrows; i++)
		{

			for(int j = 0; j < ncols; j++)
			{
				in >> data1d[i*ncols + j];

				       		         printf("%f ",data1d[i*ncols + j]);
			}
			       		        printf("\n");

		}
		in.close();

	}
	break;
	case 2:{
		nrows = 10;
		ncols = 3;
		numClusters = 3;

		ifstream in("dataset/data10x3-2.txt");

		/* start the timer for the core computation -----------------------------*/
		/* membership: the cluster id for each data object */
		idxs = (int*)malloc(nrows * sizeof(int));
		assert(idxs != NULL);
		data1d = (float*)malloc(nrows*ncols * sizeof(float));
		assert(data1d != NULL);
		if (!in) {
			cout << "Cannot open file.\n";
			exit(EXIT_SUCCESS);
		}

		for(int i = 0; i < nrows; i++)
		{

			for(int j = 0; j < ncols; j++)
			{
				in >> data1d[i*ncols + j];

				//       		         printf(" %f",data1d[i][j]);
			}
			//       		        printf("\n");

		}
		in.close();

	}
	break;
	case 3:{
			nrows = 10;
			ncols = 20;
			numClusters = 4;

			ifstream in("dataset/data10x20.txt");

			/* start the timer for the core computation -----------------------------*/
			/* membership: the cluster id for each data object */
			idxs = (int*)malloc(nrows * sizeof(int));
			assert(idxs != NULL);
			data1d = (float*)malloc(nrows*ncols * sizeof(float));
			assert(data1d != NULL);
			if (!in) {
				cout << "Cannot open file.\n";
				exit(EXIT_SUCCESS);
			}

			for(int i = 0; i < nrows; i++)
			{

				for(int j = 0; j < ncols; j++)
				{
					in >> data1d[i*ncols + j];

					       		         printf("%f ",data1d[i*ncols + j]);
				}
				       		        printf("\n");

			}
			in.close();

		}
		break;

	case 4:{
		nrows = 2048;
		ncols = 20;
		ifstream in("dataset/data2048x20.txt");

		/* start the timer for the core computation -----------------------------*/
		/* membership: the cluster id for each data object */
		idxs = (int*)malloc(nrows * sizeof(int));
		assert(idxs != NULL);
		data1d = (float*)malloc(nrows*ncols * sizeof(float));

		if (!in) {
			cout << "Cannot open file.\n";
			exit(EXIT_SUCCESS);
		}

		for(int i = 0; i < nrows; i++)
		{

			for(int j = 0; j < ncols; j++)
			{
				in >> data1d[i*ncols + j];

				//	         printf(" %f",data1d[i][j]);
			}
			//	        printf("\n");

		}
		in.close();
	}
	break;

	case 5:{
		nrows = 4096;
		ncols = 20;
		ifstream in("dataset/data4096x20.txt");

		/* start the timer for the core computation -----------------------------*/
		/* membership: the cluster id for each data object */
		idxs = (int*)malloc(nrows * sizeof(int));
		assert(idxs != NULL);
		data1d = (float*)malloc(nrows*ncols * sizeof(float));

		if (!in) {
			cout << "Cannot open file.\n";
			exit(EXIT_SUCCESS);
		}

		for(int i = 0; i < nrows; i++)
		{

			for(int j = 0; j < ncols; j++)
			{
				in >> data1d[i*ncols + j];

				//	         printf(" %f",data1d[i][j]);
			}
			//	        printf("\n");

		}
		in.close();

	}
	break;
	}


	clock_t begin = clock();


//	mex2cuw_spectral_adapter(numClusters, nrows, ncols, data1d, NULL, "euc", idxs);

	clock_t end = clock();
	double time_spent = (double)(end - begin)/(double)CLOCKS_PER_SEC;


//    // Print a message describing what the sample does.
    printf("\n***************************************************************************\n");
	printf("\nTime spent to cluster %d objects using Spectral Clustering is %f\n", nrows, time_spent);
//    printf("=%d blocks are launched!!! (%d from the GPU)\n", sum, sum-2);
    printf("***************************************************************************\n\n");

    for(int i = 0; i < nrows; i++){
    	printf("%d ", idxs[i]);
    }
    printf("\n");

    free(data1d);
    free(idxs);


//    exit(EXIT_SUCCESS);
}
