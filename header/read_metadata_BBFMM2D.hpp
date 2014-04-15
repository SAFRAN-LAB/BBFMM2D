/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	read_metadata_BBFMM2D.hpp
Read Number of unknowns, number of terms in structure,number of measurements,number of measurement sets from file
*/
#ifndef __read_metadata_BBFMM2D_hpp__
#define __read_metadata_BBFMM2D_hpp__

#include"environment.hpp"
#include "Eigen/Dense"


using namespace Eigen;
using namespace std;


void read_Metadata_BBFMM2D (const string& filenameMetadata, unsigned long& N, unsigned& m);


#endif //(__read_metadata_BBFMM2D_hpp__)
