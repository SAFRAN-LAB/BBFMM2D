/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file read_Location_Htranpose_binary.hpp
 read location and Htranspose from binary files
*/
#ifndef __read_Location_Htranpose_binary_hpp__
#define __read_Location_Htranpose_binary_hpp__

#include"H2_2D_Node.hpp"
#include"environment.hpp"

using namespace Eigen;
using namespace std;


void read_Location_Htranpose_binary(const string& filenameLocation, unsigned long N, vector<Point>& location, const string& filenameHtranspose,unsigned m, MatrixXd& Htranspose);

#endif //(__read_Location_Htranpose_binary_hpp__)
