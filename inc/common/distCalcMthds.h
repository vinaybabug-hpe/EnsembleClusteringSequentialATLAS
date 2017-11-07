/*
 ============================================================================
 Name        : distCalcMthds.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Feb 26, 2017
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

#ifndef DISTCALCMTHDS_H_
#define DISTCALCMTHDS_H_

float* distancematrix (int nrows, int ncolumns, float** centroids, float** data,
  int** mask, float weights[], int *active,char *is_centroid, char dist, int row,int transpose);

float** distancematrix (int nrows, int ncolumns, float** data,
  int** mask, float weights[], char dist, int transpose);

float* distancematrix (int nrows, int ncolumns, float** data,
  int** mask, float weights[], char dist, int row,int transpose);

float distancematrix (int n, float** data1, float** data2, int** mask1, int** mask2,
		  const float weight[], int index1, int index2,char dist, int transpose);


#endif /* DISTCALCMTHDS_H_ */
