/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file	read_Location_Htranpose_binary.cpp
 source file to read location and Htranpose matrix from binary files.
*/

#include"read_Location_Htranpose_binary.hpp"

using namespace Eigen;
using namespace std;

void read_Location_Htranpose_binary(const string& filenameLocation, unsigned long N, vector<Point>& location, const string& filenameHtranspose,unsigned m, MatrixXd& Htranspose) {
    ifstream fin;

    /* Read location */
	fin.open(filenameLocation.c_str(),ios::binary);
	
	if (!fin.good()){
		cerr << "Failed to open file " << filenameLocation << endl;
		throw runtime_error("Failed to open file!");
	}
    for (unsigned long i = 0; i < N; i++) {
        double elem =   0.;
        Point new_Point;
        fin.read((char *)&elem, sizeof(double));
        new_Point.x =   elem;
        fin.read((char *)&elem, sizeof(double));
        new_Point.y =   elem;
        location.push_back(new_Point);

    } 
	fin.close();
    
    /* Read Htranspose */
    fin.open(filenameHtranspose.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenameHtranspose << endl;
		throw runtime_error("Failed to open file!");
	}
    MatrixXd H(m,N);
    fin.read((char*) H.data(), N*m*sizeof(double));
    Htranspose  =   H.transpose();
	fin.close();
}