/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * cluster_util_bootpartition2partition.c
 *
 * Code generation for function 'cluster_util_bootpartition2partition'
 *
 */

/* Include files */


/* Function Definitions */
void cluster_util_bootpartition2partition(int n, int *b, int* bi, int *i)
{


  /* CLUSTER_UTIL_BOOTPARTITION2PARTITION - Converts boot partition to sample partition */
  /*  */
  /*    Syntax:  M = cluster_util_bootpartition2partition(X,I,metric) */
  /*  */
  /*    B      :[n,1] column of bootstrap data set indices */
  /*                  where elements are indices of original samples. E.g. if */
  /*                  B(3)= 220  then sample 220 from the original data set is */
  /*                  the third sample in the boostrapped data set */
  /*  */
  /*    BI     : [1,n] or[n,1] array of cluster labels/indices assigned to each */
  /*             boostrap data element by a clustering algorithm. */
  /*  */
  /*    I      : [1,n] or [n,1] array of cluster labels for each element of */
  /*             original data set.  E.g., if BI(3) = 4 and B(3) = 220 then the */
  /*             element 200 in the original sample was assigned to cluster 4 */
  /*             and therefore I(220) = 4 */
  /*  */
  /*  not all of original samples are included in a bootstrap data set, so */
  /*  initialize I with NaN so that these excluded elements remain NaN after */
  /*  the loop below is finished.   Note: if indices are at some point converted to  */
  /*  integers, the NaNs will become zeros.  */


  /*  main loop */
  for (int count = 0; count < n; count++) {
    i[b[count]] = bi[count];
  }
  return;
}

/* End of code generation (cluster_util_bootpartition2partition.c) */
