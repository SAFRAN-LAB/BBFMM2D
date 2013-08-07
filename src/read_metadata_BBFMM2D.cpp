/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file read_metadata_BBFMM2D.cpp
 */

#include"read_metadata_BBFMM2D.hpp"



void read_Metadata_BBFMM2D (const string& filenameMetadata, unsigned long& N, unsigned& m) {
    ifstream fin;
	fin.open(filenameMetadata.c_str());
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameMetadata << endl;
		throw runtime_error("Failed to open file!");
	}
    string line;
    getline(fin,line);
    line.erase(remove(line.begin(), line.end(), ' '),
               line.end());
    stringstream ss;
    ss << line;
    char comma;
    ss >> N >> comma >> m;
    fin.close();
}

