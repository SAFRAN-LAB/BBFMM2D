/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file	write_Into_Binary_File.cpp
 source file to read location and charges from binary files.
*/

#include"write_Into_Binary_File.hpp"

using namespace std;

void write_Into_Binary_File(const string& filename, double* data, int numOfElems) {
    
    ofstream outdata;
	outdata.open(filename.c_str(),ios::binary);
    
    outdata.write((char *)data, numOfElems*sizeof(double));
	
	outdata.close();
}