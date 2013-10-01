/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file write_Into_Binary_File.hpp
 write data into binary file
*/
#ifndef __write_Into_Binary_File_hpp__
#define __write_Into_Binary_File_hpp__

#include"environment.hpp"

using namespace std;


void write_Into_Binary_File(const string& filename, double* outdata, int numOfElems);

#endif //(__write_Into_Binary_File_hpp__)
