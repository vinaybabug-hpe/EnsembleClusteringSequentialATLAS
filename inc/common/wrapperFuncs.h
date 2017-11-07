/*
 ============================================================================
 Name        : wrapperFuncs.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Jan 30, 2016
 Version     : 1.0
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
 Description : 
 ============================================================================
 */


#ifndef WRAPPERFUNCS_H_
#define WRAPPERFUNCS_H_


#ifndef MIN
#define MIN(x,y) ((x <= y) ? x : y)
#endif

#ifndef MIN_IDX
#define MIN_IDX(x,y, idx_x, idx_y) ((x <= y) ? idx_x : idx_y)
#endif

#ifndef MAX
#define MAX(x,y) ((x >= y) ? x : y)
#endif

#if (__CUDA_ARCH__ < 200)
#define int_mult(x,y)	__mul24(x,y)
#else
#define int_mult(x,y)	x*y
#endif

#define malloc2D(name, xDim, yDim, type) do {               \
    name = (type **)malloc(xDim * sizeof(type *));          \
    assert(name != NULL);                                   \
    name[0] = (type *)malloc(xDim * yDim * sizeof(type));   \
    assert(name[0] != NULL);                                \
    for (size_t i = 1; i < xDim; i++)                       \
        name[i] = name[i-1] + yDim;                         \
} while (0)

int getModelType(char *clustMethod, char *methodList);
void getClusteringMethodsList(int length, char *methodList[], char **clustList);
int doKmeans(char **clustList);
int doKmedoids(char **clustList);
int doGMM(char **clustList);
int doSpectral(char **clustList);
int doAgglom(char **clustList);
int getMthdLstbyClust(char*clustMethod, char **methodList,
		int length, char **mthdLstbyClust);
int getClustMthdCnt(char*clustMethod, char **methodList,
		int length);
void getDistMtrcnCntrFunByKmeans(char*clustMthd, char *model,
		char *distmetric, char *centerfun);
void getDistMtrcnCntrFunByKmedoid(char*clustMthd, char *model,
		char *distmetric, char *centerfun);
void getDistMtrcnCntrFunByAgglo(char*clustMthd, char *model,
		char *distmetric, char *centerfun, char *linkcode);
void getDistMtrcnCntrFunBySpectral(char*clustMthd, char *model,
		char *distmetric, char *centerfun);

int unique_length(int *a, size_t len);

void mex2CudaWrapper();

void mex2seq_spectral_adapter(int nclusters,
							  int nrows,
							  int ncols,
							  float* data1d,
							  char* _method,
							  char* _dist,
							  int *clusterid);

void mex2seq_spectral_adapter_full_distmat(int nclusters, int nrows, int ncols, float* data1d,
	char* _method, char* _dist, int *clusterid);


void mex2seq_kmeans_adapter(int nclusters, int nrows, int ncols, float** data2d,int **mask, float *weight,
	char* _method, char* _dist, float threshold, int *clusterid);

void mex2seq_kmedians_adapter(int nclusters, int nrows, int ncols, float** data2d,int **mask, float *weight,
	char* _method, char* _dist, float threshold, int *clusterid);

int seq_gmm_main(int desired_num_clusters, float* fcs_data_by_event, int num_dimensions, int num_events, int *clusterIdx);


void mex2seq_gmm_adapter(int nclusters,
						int nrows,
						int ncols,
						float* data1d,
						char* _method,
						char* _dist,
						int *clusterid);

void agglom_adapter(int nclusters, int nrows, int ncols, float** data1d, char method, char dist, int *clusterid);

void mex2seq_agglomerative_adapter(int nclusters, int nrows, int ncols, float** data1d,
	char* _method, char* _dist, int *clusterid);

int run_qsort(unsigned int size, float *data, int debug);

bool isPow2Host(unsigned int x);

unsigned int nextPow2Host(unsigned int x);

void getNumBlocksAndThreadsHost(/*int whichKernel,*/ int n, int maxBlocks,
	int maxThreads, int &blocks, int &threads, int *maxGridSize, int *maxThreadsPerBlock);


void matrixTranspose(float *gold, float *idata, const  int size_x, const  int size_y);

void
matrixMulCPU(float *C, const float *A, const float *B, unsigned int hA, unsigned int wA, unsigned int wB);

#endif /* WRAPPERFUNCS_H_ */
