/*
 ============================================================================
 Name        : zscore.cpp
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


/* Include files */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "common/wrapper.h"
#include "common/wrapperFuncs.h"





/* Function Definitions */
void zscore(int n, float* z)
{
  float sigma;
  int k;
  float mu;


//  sigma = x[0];


  mu = 0;
  for(int count = 0; count < n; count++){
	  mu += z[count];
  }

  mu /= n;

  // setup arguments
  float init = 0;


  // compute norm
  for(int count = 0; count < n; count++){
	  float r = z[count] - mu;
	  sigma += r*r;
    }

  if (sigma == 0.0) {
    sigma = 1.0;
  }

  for (k = 0; k < n; k++) {
    z[k] =(z[k] - mu)/sigma;
  }

}


