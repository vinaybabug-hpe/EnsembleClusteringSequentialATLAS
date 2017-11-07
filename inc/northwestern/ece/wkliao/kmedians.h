/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         kmedians.h   (an OpenMP version)                          */
/*   Description:  header file for a simple k-means clustering program       */
/*                 Modified by Vinay B Gavirangaswamy for KMedians 			 */
/*                 Western Michigan University					 			 */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*            email: vinay.b.gavirangaswamy@wmich.edu                        */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _H_KMEDIANS
#define _H_KMEANS



//int omp_kmeans(int, float**, int, int, int, float, int*, float**);
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
               );

//float** file_read(int, char*, int*, int*);
//int     file_write(char*, int, int, int, float**, int*, int);
//
//int read_n_objects(int, char*, int, int, float**);
//
//int check_repeated_clusters(int, int, float**);
//
//double  wtime(void);

extern int _debug;

#endif
