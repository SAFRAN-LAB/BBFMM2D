/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file	read_Location_Charges_binary.cpp
 source file to read location and charges from binary files.
*/

#include"read_Location_Charges_binary.hpp"

using namespace std;
using namespace Eigen;

void read_Location_Charges_binary(const string& filenameLocation, unsigned long N, vector<Point>& location, const string& filenameCharges,unsigned m, double* charges) {
    ifstream fin;

    /* Read location */
	fin.open(filenameLocation.c_str(),ios::binary);
	
	if (!fin.good()){
		cerr << "Failed to open file " << filenameLocation << endl;
		throw runtime_error("Failed to open file!");
	}
    MatrixXd loc(N,2);
    fin.read((char*) loc.data(), N*2*sizeof(double));
    for (unsigned long i = 0; i < N; i++) {
        Point new_Point(loc(i,0),loc(i,1));
        location.push_back(new_Point);
    } 
	fin.close();
    
    /* Read charges */
    fin.open(filenameCharges.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenameCharges << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) charges, N*m*sizeof(double));
	fin.close();
}