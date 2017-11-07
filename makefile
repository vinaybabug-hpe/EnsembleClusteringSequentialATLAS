NAME:=cluster_ensemble_create
LDIR=lib
ODIR = obj
SDIR = src
EXTERNAL=external
CPP_ODIR = obj/cpp
CU_ODIR = obj/cpp
#LAPACK_PATH =/public/apps/lapack/3.7.0/gcc.4.7.2/lib64
OPENBLAS_PATH=/public/apps/openblas/0.2.19/gcc.4.7.2/lib

MATLAB_INSTALL_PATH := /public/apps/mdcs/R2015b
ARPACKPP_DIR = $(HOME)/Documents/cuda_ensemble/methods_fitting_cuda/ClusterEnsembleCreateSeqATLAS/external/arpack++
ARPACKPP_LIB_DIR = /public/apps/arpack-ng/gcc.4.7.2/lib
SUITESPARSE_DIR_LOCAL = $(ARPACKPP_DIR)/external/SuiteSparse
MABLAB := $(MATLAB_INSTALL_PATH)
ATLAS_INC_DIR := /public/apps/atlas/3.11.11/avx/gcc.4.7.2/include
ATLAS_LIB_DIR := /public/apps/atlas/3.11.11/avx/gcc.4.7.2/lib

include  $(ARPACKPP_DIR)/Makefile.inc

INC := -Iinc -Iexternal -I$(MABLAB)/extern/include -I$(ARPACKPP_DIR)/include -I$(SUITESPARSE_DIR_LOCAL)/CHOLMOD/Include -I$(SUITESPARSE_DIR_LOCAL)/SuiteSparse_config -I$(ATLAS_INC_DIR)
LIBS := -L$(MATLAB_INSTALL_PATH)/bin/glnxa64 -L$(ARPACKPP_DIR) -L$(ATLAS_LIB_DIR) -L$(ARPACKPP_LIB_DIR) -L$(OPENBLAS_PATH)
#-L$(LAPACK_PATH)  

#put non cuda -l* stuff here
NONCU_LIBS:=-larpack -llapack -lopenblas -lgfortran -lcblas -lf77blas -latlas

CC = gcc
CXX = g++
MEX_EXE=mex

AR = ar

CFLAGS:=-fPIC -g -fno-inline -w -c -O3 -shared -fpermissive --std=c++11 -Wall -Wextra -Drestrict=__restrict


LFLAGS := -Wall
MEX_CFLAGS := "$(MATLAB_INSTALL_PATH)/extern/lib/glnxa64/mexFunction.map"
COMPFLAGS :=-fPIC -Wall

	
#solo: mkdir_o Lib_cuda.a
#	$(NVCC) -ccbin g++  $(INC) -L$(ODIR)/ -w -m64 -dc -gencode arch=compute_50,code=sm_50 -o "obj/SpectralClustering.o" -c "src/SpectralClustering.cu" $(LIBS) $(CUDA_LIBS) $(ALL_LIBS)
#	$(NVCC) -ccbin g++ -L$(ODIR)/ -m64 -gencode arch=compute_50,code=sm_50 -o $(NAME) "obj/SpectralClustering.o"  -lcudadevrt -lmycuda

all: clean mkdir_o 
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/spectral_clustering.o $(EXTERNAL)/fastsc/spectral_clustering.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/utility_functions.o $(SDIR)/utility_functions.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/zscore.o $(SDIR)/zscore.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_ssw.o $(SDIR)/cluster_util_ssw.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/indices_count.o $(SDIR)/indices_count.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_indices2centers.o $(SDIR)/cluster_util_indices2centers.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_partition2cam.o $(SDIR)/cluster_util_partition2cam.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_ensemble2cam.o $(SDIR)/cluster_ensemble2cam.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_data_outliers.o $(SDIR)/cluster_data_outliers.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_bootpartition2partition.o $(SDIR)/cluster_util_bootpartition2partition.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/mex2cuda.o $(SDIR)/mex2cuda.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/distCalcMthds.o $(SDIR)/distCalcMthds.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_spectral.o $(EXTERNAL)/timnugent/spectral/seq_spectral.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_kmeans.o $(EXTERNAL)/northwestern/ece/wkliao/seq_kmeans.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_kmedians.o $(EXTERNAL)/northwestern/ece/wkliao/seq_kmedians.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_gmm.o $(EXTERNAL)/rochester/seq_gmm.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/invert_matrix.o $(EXTERNAL)/rochester/invert_matrix.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/AggloClustCentroidSolo.o $(EXTERNAL)/cluster_3_0/AggloClustCentroidSolo.cpp $(LIBS) $(NONCU_LIBS)			
	$(CXX) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_ensemble_create.o $(SDIR)/cluster_ensemble_create.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) -Wl,-rpath,$(MATLAB_INSTALL_PATH)/bin/glnxa64 -o $(NAME) $(ODIR)/cluster_ensemble_create.o \
	$(ODIR)/seq_spectral.o $(ODIR)/seq_kmeans.o $(ODIR)/seq_kmedians.o $(ODIR)/seq_gmm.o $(ODIR)/invert_matrix.o  \
	$(ODIR)/AggloClustCentroidSolo.o $(ODIR)/cluster_ensemble2cam.o $(ODIR)/cluster_util_partition2cam.o \
	$(ODIR)/mex2cuda.o $(ODIR)/distCalcMthds.o $(ODIR)/utility_functions.o $(ODIR)/cluster_data_outliers.o $(ODIR)/zscore.o \
	$(ODIR)/cluster_util_bootpartition2partition.o $(ODIR)/cluster_util_indices2centers.o $(ODIR)/indices_count.o \
	$(ODIR)/cluster_util_ssw.o $(ODIR)/spectral_clustering.o \
	-lmat -lmx $(INC) $(NONCU_LIBS) $(LIBS) -L$(ODIR)/ 

debug: clean mkdir_o 
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/spectral_clustering.o $(EXTERNAL)/fastsc/spectral_clustering.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/utility_functions.o $(SDIR)/utility_functions.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/zscore.o $(SDIR)/zscore.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_ssw.o $(SDIR)/cluster_util_ssw.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/indices_count.o $(SDIR)/indices_count.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_indices2centers.o $(SDIR)/cluster_util_indices2centers.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_partition2cam.o $(SDIR)/cluster_util_partition2cam.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_ensemble2cam.o $(SDIR)/cluster_ensemble2cam.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_data_outliers.o $(SDIR)/cluster_data_outliers.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_util_bootpartition2partition.o $(SDIR)/cluster_util_bootpartition2partition.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/mex2cuda.o $(SDIR)/mex2cuda.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/distCalcMthds.o $(SDIR)/distCalcMthds.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_spectral.o $(EXTERNAL)/timnugent/spectral/seq_spectral.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_kmeans.o $(EXTERNAL)/northwestern/ece/wkliao/seq_kmeans.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_kmedians.o $(EXTERNAL)/northwestern/ece/wkliao/seq_kmedians.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/seq_gmm.o $(EXTERNAL)/rochester/seq_gmm.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/invert_matrix.o $(EXTERNAL)/rochester/invert_matrix.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/AggloClustCentroidSolo.o $(EXTERNAL)/cluster_3_0/AggloClustCentroidSolo.cpp $(LIBS) $(NONCU_LIBS)			
	$(CXX) $(DFLAGS) $(CFLAGS) $(INC) -L$(ODIR)/ -c -o $(ODIR)/cluster_ensemble_create.o $(SDIR)/cluster_ensemble_create.cpp $(LIBS) $(NONCU_LIBS)
	$(CXX) $(DFLAGS) -Wl,-rpath,$(MATLAB_INSTALL_PATH)/bin/glnxa64 -o $(NAME) $(ODIR)/cluster_ensemble_create.o \
	$(ODIR)/seq_spectral.o $(ODIR)/seq_kmeans.o $(ODIR)/seq_kmedians.o $(ODIR)/seq_gmm.o $(ODIR)/invert_matrix.o  \
	$(ODIR)/AggloClustCentroidSolo.o $(ODIR)/cluster_ensemble2cam.o $(ODIR)/cluster_util_partition2cam.o \
	$(ODIR)/mex2cuda.o $(ODIR)/distCalcMthds.o $(ODIR)/utility_functions.o $(ODIR)/cluster_data_outliers.o $(ODIR)/zscore.o \
	$(ODIR)/cluster_util_bootpartition2partition.o $(ODIR)/cluster_util_indices2centers.o $(ODIR)/indices_count.o \
	$(ODIR)/cluster_util_ssw.o $(ODIR)/spectral_clustering.o \
	-lmat -lmx $(INC) $(NONCU_LIBS) $(LIBS) -L$(ODIR)/ 


$(ODIR)/%.o: $(SDIR)/%.c mkdir_o
	$(CC) -c $(CFLAGS) $(INC) -o $@ $<  

$(ODIR)/%.ou: $(SDIR)/%.cu mkdir_o
	$(NVCC) -c $(CFLAGS) $(INC) -o $@ $<


mkdir_o:	
	@mkdir -p $(ODIR)
	@mkdir -p $(ODIR)/common	
	
	
	
clean:
	rm -f $(NAME) $(ODIR)/*.o $(ODIR)/*/*.o $(ODIR)/*/*.*.o $(ODIR)/*.*.o $(OUT) $(ODIR)/*/*.a $(ODIR)/*.a *.mexa64
	
	