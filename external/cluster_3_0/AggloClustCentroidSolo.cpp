/*
 ============================================================================
 Name        : AggloClustCentroidSolo.cpp
 Author      : Vinay B Gavirangaswamy
 Created on	 : Aug 21, 2015
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

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cassert>
#include <cfloat>
#include <cmath>
#include "mat.h"

#include "common/clustLib.h"
#include "common/wrapper.h"
#include "common/wrapperFuncs.h"
#include "common/fibo.h"
#include "common/distCalcMthds.h"



/* ********************************************************************** */

static const float* sortdata = NULL; /* used in the quicksort algorithm */

/* ---------------------------------------------------------------------- */

static
int compare(const void* a, const void* b)
/* Helper function for sort. Previously, this was a nested function under
 * sort, which is not allowed under ANSI C.
 */
{ const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const float term1 = sortdata[i1];
  const float term2 = sortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

/* ---------------------------------------------------------------------- */

static void sort(int n, const float data[], int index[])
/* Sets up an index table given the data, such that data[index[]] is in
 * increasing order. Sorting is done on the indices; the array data
 * is unchanged.
 */
{ int i;
  sortdata = data;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), compare);
}

/* ********************************************************************* */
static float* getrank (int n, float data[])
/* Calculates the ranks of the elements in the array data. Two elements with
 * the same value get the same rank, equal to the average of the ranks had the
 * elements different values. The ranks are returned as a newly allocated
 * array that should be freed by the calling routine. If getrank fails due to
 * a memory allocation error, it returns NULL.
 */
{ int i;
  float* rank;
  int* index;
  rank = (float*) malloc(n*sizeof(float));
  if (!rank) return NULL;
  index = (int*)malloc(n*sizeof(int));
  if (!index)
  { free(rank);
    return NULL;
  }
  /* Call sort to get an index table */
  sort (n, data, index);
  /* Build a rank table */
  for (i = 0; i < n; i++) rank[index[i]] = i;
  /* Fix for equal ranks */
  i = 0;
  while (i < n)
  { int m;
    float value = data[index[i]];
    int j = i + 1;
    while (j < n && data[index[j]] == value) j++;
    m = j - i; /* number of equal ranks found */
    value = rank[index[i]] + (m-1)/2.;
    for (j = i; j < i + m; j++) rank[index[j]] = value;
    i += m;
  }
  free (index);
  return rank;
}


/* ******************************************************************** */
static int
makedatamask(int nrows, int ncols, float*** pdata, int*** pmask)
{ int i;
  float** data;
  int** mask;
  data = (float**)malloc(nrows*sizeof(float*));
  if(!data) return 0;
  mask = (int**)malloc(nrows*sizeof(int*));
  if(!mask)
  { free(data);
    return 0;
  }
  for (i = 0; i < nrows; i++)
  { data[i] = (float*)malloc(ncols*sizeof(float));
    if(!data[i]) break;
    mask[i] = (int*)malloc(ncols*sizeof(int));
    if(!mask[i])
    { free(data[i]);
      break;
    }
  }
  if (i==nrows) /* break not encountered */
  { *pdata = data;
    *pmask = mask;
    return 1;
  }
  *pdata = NULL;
  *pmask = NULL;
  nrows = i;
  for (i = 0; i < nrows; i++)
  { free(data[i]);
    free(mask[i]);
  }
  free(data);
  free(mask);
  return 0;
}
/* ******************************************************************** */

static
int nodecompare(const void* a, const void* b)
/* Helper function for qsort. */
{ const Node* node1 = (const Node*)a;
  const Node* node2 = (const Node*)b;
  const float term1 = node1->distance;
  const float term2 = node2->distance;
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

/* ******************************************************************** */

static
float find_closest_pair(int n, float** distmatrix, int* ip, int* jp)
/*
This function searches the distance matrix to find the pair with the shortest
distance between them. The indices of the pair are returned in ip and jp; the
distance itself is returned by the function.

n          (input) int
The number of elements in the distance matrix.

distmatrix (input) float**
A ragged array containing the distance matrix. The number of columns in each
row is one less than the row index.

ip         (output) int*
A pointer to the integer that is to receive the first index of the pair with
the shortest distance.

jp         (output) int*
A pointer to the integer that is to receive the second index of the pair with
the shortest distance.
*/
{ int i, j;
  float temp;
  float distance = distmatrix[1][0];




  *ip = 1;
  *jp = 0;
  for (i = 1; i < n; i++)
  { for (j = 0; j < i; j++)
    { temp = distmatrix[i][j];

      if (temp<distance)
      { distance = temp;
        *ip = i;
        *jp = j;

      }
    }

  }
  return distance;
}
/* ******************************************************************** */

void cuttree (int nelements, Node* tree, int nclusters, int *clusterid)

/*
Purpose
=======

The cuttree routine takes the output of a hierarchical clustering routine, and
divides the elements in the tree structure into clusters based on the
hierarchical clustering result. The number of clusters is specified by the user.

Arguments
=========

nelements      (input) int
The number of elements that were clustered.

tree           (input) Node[nelements-1]
The clustering solution. Each node in the array describes one linking event,
with tree[i].left and tree[i].right representig the elements that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1).

nclusters      (input) int
The number of clusters to be formed.

clusterid      (output) int[nelements]
The number of the cluster to which each element was assigned. Space for this
array should be allocated before calling the cuttree routine. If a memory
error occured, all elements in clusterid are set to -1.

========================================================================
*/
{ int i, j, k;
  int icluster = 0;
  const int n = nelements-nclusters; /* number of nodes to join */
  int* nodeid;
  for (i = nelements-2; i >= n; i--)
  { k = tree[i].left;
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
    k = tree[i].right;
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
  }
  nodeid = (int*) malloc(n*sizeof(int));
  if(!nodeid)
  { for (i = 0; i < nelements; i++) clusterid[i] = -1;
    return;
  }
  for (i = 0; i < n; i++) nodeid[i] = -1;
  for (i = n-1; i >= 0; i--)
  { if(nodeid[i]<0)
    { j = icluster;
      nodeid[i] = j;
      icluster++;
    }
    else j = nodeid[i];
    k = tree[i].left;
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
    k = tree[i].right;
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
  }
//  free(nodeid);
  return;
}

/* ******************************************************************** */

static
Node* pwlcluster (int nrows, int ncolumns, float** data, int** mask,
  float weight[], float** distmatrix, char dist, int transpose)

/*

Purpose
=======

The pclcluster routine performs clustering using pairwise ward-linking
on a given set of gene expression data, using the distance metric given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) float[nrows][ncolumns]
The array containing the gene expression data.

mask       (input) int[nrows][ncolumns]
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight     (input) float[ncolumns] if transpose==0;
                   float[nrows]    if transpose==1
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, and nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) float**
The distance matrix. This matrix is precalculated by the calling routine
treecluster. The pclcluster routine modifies the contents of distmatrix, but
does not deallocate it.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pclcluster returns NULL.
========================================================================
*/
{ int i, j;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  int inode;
  const int ndata = transpose ? nrows : ncolumns;
  const int nnodes = nelements - 1;

  Node* result;
  int* number;
  float** newdata;
  int** newmask;
  int* distid = (int*)malloc(nelements*sizeof(int));
  if(!distid) return NULL;
  result = (Node*)malloc(nnodes*sizeof(Node));
  if(!result)
  { free(distid);
    return NULL;
  }

  if(!makedatamask(nelements, ndata, &newdata, &newmask))
  { free(result);
    free(distid);
    return NULL;
  }

  number = (int*) malloc(nelements*sizeof(int));
    if(!number)
    {
  	printf( "10. malloc() failed");
  	free(result);
  	    free(distid);
      return NULL;
    }

  for (i = 0; i < nelements; i++){ distid[i] = i; number[j] = 1;}
  /* To remember which row/column in the distance matrix contains what */

  /* Storage for node data */
  if (transpose)
  { for (i = 0; i < nelements; i++)
    { for (j = 0; j < ndata; j++)
      { newdata[i][j] = data[j][i];
        newmask[i][j] = mask[j][i];
      }
    }
    data = newdata;
    mask = newmask;
  }
  else
  { for (i = 0; i < nelements; i++)
    { memcpy(newdata[i], data[i], ndata*sizeof(float));
      memcpy(newmask[i], mask[i], ndata*sizeof(int));
    }
    data = newdata;
    mask = newmask;
  }

  for (inode = 0; inode < nnodes; inode++)
  { /* Find the pair with the shortest distance */
    int is = 1;
    int js = 0;
    result[inode].distance = find_closest_pair(nelements-inode, distmatrix, &is, &js);
    result[inode].left = distid[js];
    result[inode].right = distid[is];

    /* Fix the distances */
     int sum = number[is] + number[js];

    /* Make node js the new node */
    for (i = 0; i < ndata; i++)
    { data[js][i] = data[js][i]*mask[js][i] + data[is][i]*mask[is][i];
      mask[js][i] += mask[is][i];
      if (mask[js][i]) data[js][i] /= mask[js][i];
    }
    free(data[is]);
    free(mask[is]);
    data[is] = data[nnodes-inode];
    mask[is] = mask[nnodes-inode];

    /* Fix the distances */
    distid[is] = distid[nnodes-inode];
    for (i = 0; i < is; i++)
      distmatrix[is][i] = distmatrix[nnodes-inode][i];
    for (i = is + 1; i < nnodes-inode; i++)
      distmatrix[i][is] = distmatrix[nnodes-inode][i];

    /* Update number of elements in the clusters */
    number[js] = sum;
    number[is] = number[nnodes-inode];

    distid[js] = -inode-1;
    for (i = 0; i < js; i++){
      distmatrix[js][i] = distancematrix(ndata,data,data,mask,mask,weight,js,i,dist, 0);
      distmatrix[js][i] = ((float)(number[is]*number[js])/(float)(number[is]+number[js]))* fabs(distmatrix[js][i]*distmatrix[js][i]);
    }
    for (i = js + 1; i < nnodes-inode; i++){
      distmatrix[i][js] = distancematrix(ndata,data,data,mask,mask,weight,js,i,dist, 0);
      distmatrix[i][js] = ((float)(number[is]*number[js])/(float)(number[is]+number[js]))* fabs(distmatrix[i][js]*distmatrix[i][js]);
    }
  }

  /* Free temporarily allocated space */
  free(data[0]);
  free(mask[0]);
  free(data);
  free(mask);
  free(distid);
  free(number);

  return result;
}

/* ******************************************************************** */
static
Node* pclcluster (int nrows, int ncolumns, float** data, int** mask,
  float weight[], float** distmatrix, char dist, int transpose)

/*

Purpose
=======

The pclcluster routine performs clustering using pairwise centroid-linking
on a given set of gene expression data, using the distance metric given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) float[nrows][ncolumns]
The array containing the gene expression data.

mask       (input) int[nrows][ncolumns]
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight     (input) float[ncolumns] if transpose==0;
                   float[nrows]    if transpose==1
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, and nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) float**
The distance matrix. This matrix is precalculated by the calling routine
treecluster. The pclcluster routine modifies the contents of distmatrix, but
does not deallocate it.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pclcluster returns NULL.
========================================================================
*/
{ int i, j;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  int inode;
  const int ndata = transpose ? nrows : ncolumns;
  const int nnodes = nelements - 1;



  Node* result;
  float** newdata;
  int** newmask;
  int* distid = (int*)malloc(nelements*sizeof(int));
  if(!distid) return NULL;
  result = (Node*)malloc(nnodes*sizeof(Node));
  if(!result)
  { free(distid);
    return NULL;
  }
  if(!makedatamask(nelements, ndata, &newdata, &newmask))
  { free(result);
    free(distid);
    return NULL;
  }

  for (i = 0; i < nelements; i++) distid[i] = i;
  /* To remember which row/column in the distance matrix contains what */

  /* Storage for node data */
  if (transpose)
  { for (i = 0; i < nelements; i++)
    { for (j = 0; j < ndata; j++)
      { newdata[i][j] = data[j][i];
        newmask[i][j] = mask[j][i];
      }
    }
    data = newdata;
    mask = newmask;
  }
  else
  { for (i = 0; i < nelements; i++)
    { memcpy(newdata[i], data[i], ndata*sizeof(float));
      memcpy(newmask[i], mask[i], ndata*sizeof(int));
    }
    data = newdata;
    mask = newmask;
  }
//  for (int inode = 0, iteration=0; inode < nnodes, iteration < 11; inode++, iteration++)
  for (inode = 0; inode < nnodes; inode++)
  { /* Find the pair with the shortest distance */
    int is = 1;
    int js = 0;
    result[inode].distance = find_closest_pair(nelements-inode, distmatrix, &is, &js);
    result[inode].left = distid[js];
    result[inode].right = distid[is];


    /* Make node js the new node */
    for (i = 0; i < ndata; i++)
    { data[js][i] = data[js][i]*mask[js][i] + data[is][i]*mask[is][i];
      mask[js][i] += mask[is][i];
      if (mask[js][i]) data[js][i] /= mask[js][i];
    }


    free(data[is]);
    free(mask[is]);
    data[is] = data[nnodes-inode];
    mask[is] = mask[nnodes-inode];


    /* Fix the distances */
    distid[is] = distid[nnodes-inode];
    for (i = 0; i < is; i++)
      distmatrix[is][i] = distmatrix[nnodes-inode][i];


    for (i = is + 1; i < nnodes-inode; i++)
      distmatrix[i][is] = distmatrix[nnodes-inode][i];


    distid[js] = -inode-1;
    for (i = 0; i < js; i++)
      distmatrix[js][i] = distancematrix(ndata,data,data,mask,mask,weight,js,i,dist,0);


    for (i = js + 1; i < nnodes-inode; i++)
      distmatrix[i][js] = distancematrix(ndata,data,data,mask,mask,weight,js,i,dist,0);


  }

  /* Free temporarily allocated space */
  free(data[0]);
  free(mask[0]);
  free(data);
  free(mask);
  free(distid);

  return result;
}


/* ******************************************************************** */

static Node* pmlcluster  (int nrows, int ncolumns, float** data, int** mask,
		  float weight[], float** distmatrix, char dist, int transpose)
/*

Purpose
=======

The pmlcluster routine performs clustering using pairwise maximum- (complete-)
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) float**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pmlcluster returns NULL.
========================================================================
*/
{ int j;
  int n;
  int* clusterid;
  Node* result;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  clusterid = (int*)malloc(nelements*sizeof(int));
  if(!clusterid){printf( "7. malloc() failed"); return NULL;}
  result = (Node*) malloc((nelements-1)*sizeof(Node));
  if (!result)
  {
	printf( "8. malloc() failed");
	free(clusterid);
    return NULL;
  }



  /* Setup a list specifying to which cluster a gene belongs */
  for (j = 0; j < nelements; j++) clusterid[j] = j;

  int iteration;
//  for (n = nelements, iteration=0; n > 1 && iteration ==0; n--, iteration++)
  for (n = nelements; n > 1; n--)
  {

	int is = 1;
    int js = 0;
    result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

//	printf("\n[%d] %f -------->[%dx%d]",n, result[nelements-n].distance, is, js);

    /* Fix the distances */
    for (j = 0; j < js; j++)
      distmatrix[js][j] = std::max(distmatrix[is][j],distmatrix[js][j]);
    for (j = js+1; j < is; j++)
      distmatrix[j][js] = std::max(distmatrix[is][j],distmatrix[j][js]);
    for (j = is+1; j < n; j++)
      distmatrix[j][js] = std::max(distmatrix[j][is],distmatrix[j][js]);

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update clusterids */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
  }
  free(clusterid);

  return result;
}

/* ******************************************************************* */

static Node* palcluster  (int nrows, int ncolumns, float** data, int** mask,
		  float weight[], float** distmatrix, char dist, int transpose)
/*
Purpose
=======

The palcluster routine performs clustering using pairwise average
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) float**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, palcluster returns NULL.
========================================================================
*/
{ int i, j;
  int n;
  int* clusterid;
  int* number;
  Node* result;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  clusterid = (int*) malloc(nelements*sizeof(int));
  if(!clusterid) {printf( "9. malloc() failed"); return NULL;}
  number = (int*) malloc(nelements*sizeof(int));
  if(!number)
  {
	printf( "10. malloc() failed");
	free(clusterid);
    return NULL;
  }
  result = (Node*) malloc((nelements-1)*sizeof(Node));
  if (!result)
  {
	printf( "11. malloc() failed");
	free(clusterid);
    free(number);
    return NULL;
  }

  /* Setup a list specifying to which cluster a gene belongs, and keep track
   * of the number of elements in each cluster (needed to calculate the
   * average). */
  for (j = 0; j < nelements; j++)
  { number[j] = 1;
    clusterid[j] = j;
  }



  for (n = nelements; n > 1; n--)
  {


	int sum;
    int is = 1;
    int js = 0;
    result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

    /* Save result */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];

    /* Fix the distances */
    sum = number[is] + number[js];
    for (j = 0; j < js; j++)
    { distmatrix[js][j] = distmatrix[is][j]*number[is]
                        + distmatrix[js][j]*number[js];
      distmatrix[js][j] /= sum;
    }
    for (j = js+1; j < is; j++)
    { distmatrix[j][js] = distmatrix[is][j]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }
    for (j = is+1; j < n; j++)
    { distmatrix[j][js] = distmatrix[j][is]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update number of elements in the clusters */
    number[js] = sum;
    number[is] = number[n-1];

    /* Update clusterids */
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];


  }
  free(clusterid);
  free(number);

  return result;
}


/* ******************************************************************** */

static
Node* pslcluster (int nrows, int ncolumns, float** data, int** mask,
  float weight[], float** distmatrix, char dist, int transpose)

/*

Purpose
=======

The pslcluster routine performs single-linkage hierarchical clustering, using
either the distance matrix directly, if available, or by calculating the
distances from the data array. This implementation is based on the SLINK
algorithm, described in:
Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link
cluster method. The Computer Journal, 16(1): 30-34.
The output of this algorithm is identical to conventional single-linkage
hierarchical clustering, but is much more memory-efficient and faster. Hence,
it can be applied to large data sets, for which the conventional single-
linkage algorithm fails due to lack of memory.


Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) float[nrows][ncolumns]
The array containing the gene expression data.

mask       (input) int[nrows][ncolumns]
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, and nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) float**
The distance matrix. If the distance matrix is passed by the calling routine
treecluster, it is used by pslcluster to speed up the clustering calculation.
The pslcluster routine does not modify the contents of distmatrix, and does
not deallocate it. If distmatrix is NULL, the pairwise distances are calculated
by the pslcluster routine from the gene expression data (the data and mask
arrays) and stored in temporary arrays. If distmatrix is passed, the original
gene expression data (specified by the data and mask arguments) are not needed
and are therefore ignored.


Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pslcluster returns NULL.

========================================================================
*/
{ int i, j, k;
  const int nelements = transpose ? ncolumns : nrows;
  const int nnodes = nelements - 1;
  int* vector;
  float* temp;
  int* index;
  Node* result;
  temp = (float*) malloc(nnodes*sizeof(float));
  if(!temp) { printf( "3. malloc() failed"); return NULL;}
  index = (int*) malloc(nelements*sizeof(int));
  if(!index)
  {
	  printf( "4. malloc() failed");
	  free(temp);
    return NULL;
  }
  vector = (int*) malloc(nnodes*sizeof(int));
  if(!vector)
  {
	printf( "5. malloc() failed");
	free(index);
    free(temp);
    return NULL;
  }
  result = (Node*) malloc(nelements*sizeof(Node));
  if(!result)
  {
	printf( "6. malloc() failed");
	free(vector);
    free(index);
    free(temp);
    return NULL;
  }

  for (i = 0; i < nnodes; i++) vector[i] = i;

  if(distmatrix)
  { for (i = 0; i < nrows; i++)
    { result[i].distance = DBL_MAX;
      for (j = 0; j < i; j++) temp[j] = distmatrix[i][j];
      for (j = 0; j < i; j++)
      { k = vector[j];
        if (result[j].distance >= temp[j])
        { if (result[j].distance < temp[k]) temp[k] = result[j].distance;
          result[j].distance = temp[j];
          vector[j] = i;
        }
        else if (temp[j] < temp[k]) temp[k] = temp[j];
      }
      for (j = 0; j < i; j++)
      {
        if (result[j].distance >= result[vector[j]].distance) vector[j] = i;
      }
    }
  }

  free(temp);

  for (i = 0; i < nnodes; i++) result[i].left = i;
  qsort(result, nnodes, sizeof(Node), nodecompare);

  for (i = 0; i < nelements; i++) index[i] = i;
  for (i = 0; i < nnodes; i++)
  { j = result[i].left;
    k = vector[j];
    result[i].left = index[j];
    result[i].right = index[k];
    index[k] = -i-1;
  }
  free(vector);
  free(index);


  //result = (Node*) realloc(result, nnodes*sizeof(Node));

  return result;
}


/* ******************************************************************* */

Node* treecluster (int nrows, int ncolumns, float** data, int** mask,
  float weight[], int transpose, char dist, char method, float** distmatrix)
{ Node* result = NULL;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  const int ldistmatrix = (distmatrix==NULL) ? 1 : 0;

  if (nelements < 2) return NULL;

  /* Calculate the distance matrix if the user didn't give it */
  if(ldistmatrix)
  { distmatrix =
      distancematrix(nrows, ncolumns, data, mask, weight, dist, transpose);
    if (!distmatrix) return NULL; /* Insufficient memory */


  }

  switch(method)
  {
  	  case 's':
 	        result = pslcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist, transpose);
      break;
  	  case 'm':
  	        result = pmlcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist, transpose);
  	  break;
  	  case 'c':
	        result = pclcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist, transpose);
      break;
  	  case 'a':
  		  	result = palcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist, transpose);
      break;
      case 'w':
    	  	result = pwlcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist, transpose);
      break;
  }

//    for (int i = 0; i < nelements; i++){
//  	  for (int j = 0; j < i; j++){
//  		if(distmatrix[i][j]==0)printf("%.4f ", distmatrix[i][j]);
//  		else printf("     ");
//  	  }
//  	  	printf("\n");
//    }


  /* Deallocate space for distance matrix, if it was allocated by treecluster */
  if(ldistmatrix)
  { int i;
    for (i = 1; i < nelements; i++) free(distmatrix[i]);
    free (distmatrix);
  }


  return result;
}


/* ========================================================================= */

void agglom_adapter(int nclusters, int nrows, int ncols, float** data1d, char method, char dist, int *clusterid)
/* Perform hierarchical clustering on data */
{
	int i, j;
	int transpose = 0;
	const int nelements = (transpose == 0) ? nrows : ncols;
	Node* tree = NULL;
	float weight[ncols];
	int** cmask;


	for (int col = 0; col < ncols; col++){
		weight[col] = 1;
	}

	cmask = (int**)malloc(nrows * sizeof(int*));

	for (int row = 0; row < nrows; row++){
		cmask[row] = (int*)malloc(ncols* sizeof(int));
		for (int col = 0; col < ncols; col++){
			cmask[row][col] = data1d[row][col] == 0? 0 : 1;
		}
	}

//	tree = (Node*) malloc(nelements*sizeof(Node));
//	  if(!tree)
//	  {
//		printf( "2. malloc() failed");
//	    return;
//	  }



	tree = treecluster(nrows, ncols, data1d, cmask, weight, 0, dist, method, NULL);


	if (!tree) { /* Indication that the treecluster routine failed */
		printf("\ntreecluster routine failed due to insufficient memory\n");
		//free(weight);
		return;
	}


//	printf("\n=============== Cutting a hierarchical clustering tree ==========\n");


	cuttree(nrows, tree, nclusters, clusterid);

//	for(int i=0; i<nrows; i++)
//			printf("%2d ", clusterid[i]);
//		printf("\n");



	return;
}

