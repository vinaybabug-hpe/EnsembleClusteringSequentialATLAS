/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * my_count.c
 *
 * Code generation for function 'my_count'
 *
 */

/* Include files */
#include <iostream>
#include <unordered_set>
#include <set>


/* Function Definitions */
void indices_count(int *X, int nSamples, int *valueList, int nClusters)
{

  /* MY_COUNT Computes the counts of each unique value in X */
  /*    INPUT ARGS */
  /*        X           :  column or row vector of values */
  /*        valueList   :  list of unique values to be counted */
  /*                       OPTIONAL: if not provided, counts are based on data in X */
  /*  */
  /*    OUTPUT ARGS */
  /*        counts      :  count of each value */
  /*        valueList   :  list of values associated with each count */
  /*                       will be in same order as valueList if valueList given */
  /*                       as argument, otherwise will be in sorted order. */
  /*  Ensure X is a column vector */
  /*  if size(X,2)>1 */
  /*      X=X'; */
  /*  end */
  /*  If value list not provided, base list on unique values...in sorted order */
  /*  if ~exist('valueList','var') */

	std::unordered_set<int> s(X, X + nSamples);
	std::set<int> nValues(s.begin(), s.end());

	for(int iValue = 0; iValue < nClusters; iValue++){
		int label = *std::next(nValues.begin(), iValue);
		int counter = 0;
		// doing count if
		for(int count = 0; count < nSamples; count++){
			if(label == X[count]){
				counter++;
			}
		}
		valueList[iValue] = counter;
	}

	s.clear();
	nValues.clear();
	return;
}

/* End of code generation (my_count.c) */
