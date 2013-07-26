//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	read_location_H.cpp
//

#include"read_location_H.hpp"

void read_Location_and_Measurement_operator (const string& filename, unsigned long& N, VectorXd* location, unsigned& m, MatrixXd& Htranspose) {
    ifstream fin;
	fin.open(filename.c_str());
	
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    // read first line in the file : N, m
    string line;
    getline(fin,line);
    line.erase(remove(line.begin(), line.end(), ' '),
               line.end());
    stringstream ss;
    ss << line;
    char nonuse;
    ss >> N >> nonuse >> m;
    Htranspose  =   MatrixXd::Zero(N,m);
    location[0]	=	VectorXd::Random(N);
	location[1]	=	VectorXd::Random(N);
    
    // read location and measurement operator
    unsigned long row = 0;
    while(getline(fin, line)){
        line.erase(remove(line.begin(), line.end(), ' '),
                   line.end());
        int i=0;
        while (line[i]!=',') {
            i++;
        }
        location[0](row) = (double)atof(&line[1]);
        location[1](row) = (double)atof(&line[i+1]);
 
        while (line[i]!=')') {
            i++;
        }
        if (line.length()!=(unsigned)i+2) {
            read_Measurement_operator(line.substr(i+1), row, Htranspose, m);
        }
        row++;
    }
    fin.close();
}

void read_Measurement_operator(const string& s, unsigned long row, MatrixXd& Htranspose, unsigned m) {
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
                    Htranspose(row,k)=(double)atof(start_pt);
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


