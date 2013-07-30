//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	H2_2D_tree.hpp
//
#ifndef __kernel_types_hpp__
#define __kerbel_types_hpp__

#include"environment.hpp"
#include"kernel_base.hpp"

using namespace Eigen;
using namespace std;


class kernel_Logarithm: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};

class kernel_OneOverR2: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};

class kernel_Gaussian: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};

class kernel_Quadric: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};

class kernel_InverseQuadric: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};

class kernel_ThinPlateSpline: public kernel_base {
public:
    virtual double kernel_func(Point r0, Point r1);
};


#endif //(__kerbel_types_hpp__)
