//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran</author>
//	
//	kernel.hpp
//
#ifndef __kernel_hpp__
#define __kernel_hpp__

#include"cmath"
#include"Eigen/Dense"

using namespace Eigen;

//	Evaluate the kernel
void kernel2D(const unsigned long M, const VectorXd* x, const unsigned long N, const VectorXd* y, MatrixXd& Kernel);

#endif
