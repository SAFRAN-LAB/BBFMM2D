//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	read_location_H.hpp
//
#ifndef __read_location_H_hpp__
#define __read_location_H_hpp__

#include"iostream"
#include <sstream>
#include<fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <string>
#include <stdlib.h>
#include"cmath"
#include"Eigen/Dense"


using namespace Eigen;
using namespace std;


void read_Location_and_Measurement_operator (const string& filename, unsigned long& N, VectorXd* location, unsigned& m, MatrixXd& H);

void read_Measurement_operator(const string& s, unsigned long row, MatrixXd& H, unsigned m);

#endif //(__read_location_H_hpp__)
