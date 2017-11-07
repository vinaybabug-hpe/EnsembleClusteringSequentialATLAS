/*
 ============================================================================
 Name        : cluster_data_outliers.h
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
 


#ifndef CLUSTER_DATA_OUTLIERS_H
#define CLUSTER_DATA_OUTLIERS_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>



/* Function Declarations */
float* cluster_data_outliers(int nrows, int ncols, float* X, char* distFun, float cutoff, int *newNElements);

#endif

/* End of code generation (cluster_data_outliers.h) */
