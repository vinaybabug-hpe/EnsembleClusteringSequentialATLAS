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

void sort(int n, const float data[], int index[])
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

/* ********************************************************************* */

static
float euclid (int n, float** data1, float** data2, int** mask1, int** mask2,
  const float weight[], int index1, int index2, int transpose)

/*
Purpose
=======

The euclid routine calculates the weighted Euclidean distance between two
rows or columns in a matrix.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.

============================================================================
*/
{ float result = 0.;
  float tweight = 0;
  int i;
  if (transpose==0) /* Calculate the distance between two rows */
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term = data1[index1][i] - data2[index2][i];
        result += weight[i]*term*term;
        tweight += weight[i];
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term = data1[i][index1] - data2[i][index2];
        result += weight[i]*term*term;
        tweight += weight[i];
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result /= tweight;
  return result;
}

/* ********************************************************************* */

static
float cityblock (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)

/*
Purpose
=======

The cityblock routine calculates the weighted "City Block" distance between
two rows or columns in a matrix. City Block distance is defined as the
absolute value of X1-X2 plus the absolute value of Y1-Y2 plus..., which is
equivalent to taking an "up and over" path.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.

============================================================================ */
{ float result = 0.;
  float tweight = 0;
  int i;
  if (transpose==0) /* Calculate the distance between two rows */
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term = data1[index1][i] - data2[index2][i];
        result = result + weight[i]*fabs(term);
        tweight += weight[i];
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term = data1[i][index1] - data2[i][index2];
        result = result + weight[i]*fabs(term);
        tweight += weight[i];
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result /= tweight;
  return result;
}

/* ********************************************************************* */

static
float correlation (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The correlation routine calculates the weighted Pearson distance between two
rows or columns in a matrix. We define the Pearson distance as one minus the
Pearson correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ float result = 0.;
  float sum1 = 0.;
  float sum2 = 0.;
  float denom1 = 0.;
  float denom2 = 0.;
  float tweight = 0.;
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term1 = data1[index1][i];
        float term2 = data2[index2][i];
        float w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term1 = data1[i][index1];
        float term2 = data2[i][index2];
        float w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result -= sum1 * sum2 / tweight;
  denom1 -= sum1 * sum1 / tweight;
  denom2 -= sum2 * sum2 / tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
float acorrelation (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The acorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ float result = 0.;
  float sum1 = 0.;
  float sum2 = 0.;
  float denom1 = 0.;
  float denom2 = 0.;
  float tweight = 0.;
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term1 = data1[index1][i];
        float term2 = data2[index2][i];
        float w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term1 = data1[i][index1];
        float term2 = data2[i][index2];
        float w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result -= sum1 * sum2 / tweight;
  denom1 -= sum1 * sum1 / tweight;
  denom2 -= sum2 * sum2 / tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
float ucorrelation (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The ucorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the uncentered version of the Pearson correlation. In the
uncentered Pearson correlation, a zero mean is used for both vectors even if
the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ float result = 0.;
  float denom1 = 0.;
  float denom2 = 0.;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term1 = data1[index1][i];
        float term2 = data2[index2][i];
        float w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term1 = data1[i][index1];
        float term2 = data2[i][index2];
        float w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  if (!flag) return 0.;
  if (denom1==0.) return 1.;
  if (denom2==0.) return 1.;
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
float uacorrelation (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The uacorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the uncentered version of the
Pearson correlation. In the uncentered Pearson correlation, a zero mean is used
for both vectors even if the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ float result = 0.;
  float denom1 = 0.;
  float denom2 = 0.;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { float term1 = data1[index1][i];
        float term2 = data2[index2][i];
        float w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { float term1 = data1[i][index1];
        float term2 = data2[i][index2];
        float w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  if (!flag) return 0.;
  if (denom1==0.) return 1.;
  if (denom2==0.) return 1.;
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* *********************************************************************  */

static
float spearman (int n, float** data1, float** data2, int** mask1,
  int** mask2, const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The spearman routine calculates the Spearman distance between two rows or
columns. The Spearman distance is defined as one minus the Spearman rank
correlation.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
These weights are ignored, but included for consistency with other distance
measures.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ int i;
  int m = 0;
  float* rank1;
  float* rank2;
  float result = 0.;
  float denom1 = 0.;
  float denom2 = 0.;
  float avgrank;
  float* tdata1;
  float* tdata2;
  tdata1 = (float*)malloc(n*sizeof(float));
  if(!tdata1) return 0.0; /* Memory allocation error */
  tdata2 = (float*)malloc(n*sizeof(float));
  if(!tdata2) /* Memory allocation error */
  { free(tdata1);
    return 0.0;
  }
  if (transpose==0)
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { tdata1[m] = data1[index1][i];
        tdata2[m] = data2[index2][i];
        m++;
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { tdata1[m] = data1[i][index1];
        tdata2[m] = data2[i][index2];
        m++;
      }
    }
  }
  if (m==0)
  { free(tdata1);
    free(tdata2);
    return 0;
  }
  rank1 = getrank(m, tdata1);
  free(tdata1);
  if(!rank1)
  { free(tdata2);
    return 0.0; /* Memory allocation error */
  }
  rank2 = getrank(m, tdata2);
  free(tdata2);
  if(!rank2) /* Memory allocation error */
  { free(rank1);
    return 0.0;
  }
  avgrank = 0.5*(m-1); /* Average rank */
  for (i = 0; i < m; i++)
  { const float value1 = rank1[i];
    const float value2 = rank2[i];
    result += value1 * value2;
    denom1 += value1 * value1;
    denom2 += value2 * value2;
  }
  /* Note: denom1 and denom2 cannot be calculated directly from the number
   * of elements. If two elements have the same rank, the squared sum of
   * their ranks will change.
   */
  free(rank1);
  free(rank2);
  result /= m;
  denom1 /= m;
  denom2 /= m;
  result -= avgrank * avgrank;
  denom1 -= avgrank * avgrank;
  denom2 -= avgrank * avgrank;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* *********************************************************************  */

static
float kendall (int n, float** data1, float** data2, int** mask1, int** mask2,
  const float weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The kendall routine calculates the Kendall distance between two
rows or columns. The Kendall distance is defined as one minus Kendall's tau.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) float array
The data array containing the first vector.

data2  (input) float array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) float[n]
These weights are ignored, but included for consistency with other distance
measures.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ int con = 0;
  int dis = 0;
  int exx = 0;
  int exy = 0;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  float denomx;
  float denomy;
  float tau;
  int i, j;
  if (transpose==0)
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { for (j = 0; j < i; j++)
        { if (mask1[index1][j] && mask2[index2][j])
          { float x1 = data1[index1][i];
            float x2 = data1[index1][j];
            float y1 = data2[index2][i];
            float y2 = data2[index2][j];
            if (x1 < x2 && y1 < y2) con++;
            if (x1 > x2 && y1 > y2) con++;
            if (x1 < x2 && y1 > y2) dis++;
            if (x1 > x2 && y1 < y2) dis++;
            if (x1 == x2 && y1 != y2) exx++;
            if (x1 != x2 && y1 == y2) exy++;
            flag = 1;
          }
        }
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { for (j = 0; j < i; j++)
        { if (mask1[j][index1] && mask2[j][index2])
          { float x1 = data1[i][index1];
            float x2 = data1[j][index1];
            float y1 = data2[i][index2];
            float y2 = data2[j][index2];
            if (x1 < x2 && y1 < y2) con++;
            if (x1 > x2 && y1 > y2) con++;
            if (x1 < x2 && y1 > y2) dis++;
            if (x1 > x2 && y1 < y2) dis++;
            if (x1 == x2 && y1 != y2) exx++;
            if (x1 != x2 && y1 == y2) exy++;
            flag = 1;
          }
        }
      }
    }
  }
  if (!flag) return 0.;
  denomx = con + dis + exx;
  denomy = con + dis + exy;
  if (denomx==0) return 1;
  if (denomy==0) return 1;
  tau = (con-dis)/sqrt(denomx*denomy);
  return 1.-tau;
}

/* *********************************************************************  */

static float(*setmetric(char dist))
  (int, float**, float**, int**, int**, const float[], int, int, int)
{ switch(dist)
  { case 'e': return &euclid;
    case 'b': return &cityblock;
    case 'c': return &correlation;
    case 'a': return &acorrelation;
    case 'u': return &ucorrelation;
    case 'x': return &uacorrelation;
    case 's': return &spearman;
    case 'k': return &kendall;
    default: return &euclid;
  }
  return NULL; /* Never get here */
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

float distancematrix (int n, float** data1, float** data2, int** mask1, int** mask2,
		  const float weight[], int index1, int index2,char dist, int transpose)

{
	/* Set the metric function as indicated by dist */
	float (*metric)(int, float**, float**, int**, int**, const float[], int, int, int) = setmetric(dist);
  return metric(n,data1,data2,mask1,mask2,weight,index1,index2,transpose);
}
/* ******************************************************************** */

float* distancematrix (int nrows, int ncolumns, float** data,
  int** mask, float weights[], char dist, int row,int transpose)

{ /* First determine the size of the distance matrix */
  const int n = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;
  float* matrix;

  /* Set the metric function as indicated by dist */
  float (*metric)
    (int, float**, float**, int**, int**, const float[], int, int, int) =
       setmetric(dist);

  if (n < 2) return NULL;

  /* Set up the ragged array */
  matrix = (float*) malloc(n*sizeof(float));
  if(matrix==NULL) return NULL; /* Not enough memory available */


  /* Calculate the distances and save them in the ragged array */
  for (i = 0; i < n; i++){
	  if(i == row){
		  matrix[i] = FLT_MAX;
	  }else{
		  matrix[i]=metric(ndata,data,data,mask,mask,weights,row,i,transpose);
	  }
  }

  return matrix;
}

/* ******************************************************************** */


float** distancematrix (int nrows, int ncolumns, float** data,
  int** mask, float weights[], char dist, int transpose)

{ /* First determine the size of the distance matrix */
  const int n = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;
  float** matrix;

  /* Set the metric function as indicated by dist */
  float (*metric)
    (int, float**, float**, int**, int**, const float[], int, int, int) =
       setmetric(dist);

  if (n < 2) return NULL;

  /* Set up the ragged array */
  matrix = (float**) malloc(n*sizeof(float*));
  if(matrix==NULL) return NULL; /* Not enough memory available */
  matrix[0] = NULL;
  /* The zeroth row has zero columns. We allocate it anyway for convenience.*/
  for (i = 1; i < n; i++)
  { matrix[i] = (float*)malloc(i*sizeof(float));
    if (matrix[i]==NULL) break; /* Not enough memory available */
  }
  if (i < n) /* break condition encountered */
  { j = i;
    for (i = 1; i < j; i++) free(matrix[i]);
    return NULL;
  }

  /* Calculate the distances and save them in the ragged array */
  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++)
      matrix[i][j]=metric(ndata,data,data,mask,mask,weights,i,j,transpose);

  return matrix;
}
/* ******************************************************************** */

float* distancematrix (int nrows, int ncolumns, float** centroids, float** data,
  int** mask, float weights[], int *active,char *is_centroid, char dist, int row,int transpose)

{ /* First determine the size of the distance matrix */
  const int n = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;
  float* matrix;

  /* Set the metric function as indicated by dist */
  float (*metric)
    (int, float**, float**, int**, int**, const float[], int, int, int) =
       setmetric(dist);

  if (n < 2) return NULL;

  /* Set up the ragged array */
  matrix = (float*) malloc(n*sizeof(float));
  if(matrix==NULL) return NULL; /* Not enough memory available */


  /* Calculate the distances and save them in the ragged array */
  for (i = 0; i < n; i++){
	  if(i == row || active[i] == row){
		  matrix[i] = FLT_MAX;
	  }else{
		  if(is_centroid[i] == 1 || active[i] != i){


			  matrix[i] = metric(ndata,centroids,centroids,mask,mask,weights,row,active[i],transpose);

//			  printf("\n[%d to %d] calc with centroid %f", row, i, matrix[i]);

		  }
		  else{

			  matrix[i]=metric(ndata,centroids,data,mask,mask,weights,row, i,transpose);

//			  printf("\n[%d to %d] calc with data %f active[i] = %d", row, i, matrix[i], active[i]);

		  }

	  }
  }

  return matrix;
}



