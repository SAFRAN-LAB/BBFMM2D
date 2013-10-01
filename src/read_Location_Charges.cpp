/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file	read_Location_Charges.cpp
*/

#include"read_Location_Charges.hpp"

using namespace Eigen;


void read_Charges(const string& s, unsigned long row, double*& charges, const unsigned m, const unsigned long N);

void read_Location_Charges (const string& filename, unsigned long N, vector<Point>& location, unsigned m, double*& charges) {
    ifstream fin;
	fin.open(filename.c_str());
	
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    
    string line;
        
    // read location and measurement operator
    unsigned long row = 0;
    while(getline(fin, line)){
        line.erase(remove(line.begin(), line.end(), ' '),
                   line.end());
        int i=0;
        while (line[i]!=',') {
            i++;
        }
        Point new_Point((double)atof(&line[1]), (double)atof(&line[i+1]));
        location.push_back(new_Point);
 
        while (line[i]!=')') {
            i++;
        }
        if (line.length()!=(unsigned)i+2) {
            read_Charges(line.substr(i+1), row, charges, m, N);
        }
        row++;
    }
    fin.close();
}

void read_Charges(const string& s, unsigned long row, double*& charges, const unsigned m, const unsigned long N) {
    if (!s.empty()) {
        unsigned k = 0;
        const char* start_pt = NULL;
        for (unsigned i = 0; i < s.length();)
            if ( s[i]==',' || s[i]=='(' ) {
                start_pt=&s[++i];
                while(s[i]!=',' && s[i]!=')') {
                    i++;
                }
                if(start_pt!=&s[i]) {
                    charges[k*N+row]=(double)atof(start_pt);
                }
                k++;
            }
            else {
                i++;
            }
        if(k!=m)
            throw runtime_error("Number of measurement is not consistent with input");
    }
}


