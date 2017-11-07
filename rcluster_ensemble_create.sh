#!/bin/bash 

#make sure to change blas and lapack to correct alternative
# sudo update-alternatives --config liblapack.so.3
# sudo update-alternatives --config libblas.so.3

#call executable with .mat files as arguments
#./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_50.mat "Sample Description" dataset/testcase1/output_idxs_50.mat
./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_100x50.mat "Sample Description" dataset/testcase1/output_idxs.mat
#./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_200x50.mat "Sample Description" dataset/testcase1/output_idxs_200_50.mat
#./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_300x50.mat "Sample Description" dataset/testcase1/output_idxs_300_50.mat
#./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_1000x50.mat "Sample Description" dataset/testcase1/output_idxs_1000_50.mat
#./cluster_ensemble_create dataset/testcase1/PST_all.mat dataset/testcase1/X_4096x100.mat "Sample Description" dataset/testcase1/output_idxs_4096_100.mat


# Risky Decision Making Data From IGT
#./cluster_ensemble_create dataset/testcase1/PST_all_nBoot_1.mat dataset/testcase1/ens_combined_4x5_512.mat "Risky Decision Making Data From IGT" dataset/testcase1/output_ens_combined_4x5_512.mat