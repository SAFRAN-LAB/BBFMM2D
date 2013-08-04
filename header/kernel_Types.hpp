/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!	\file kernel_Types.hpp
  Defines differnt types of kernel
*/
#ifndef __kernel_Types_hpp__
#define __kernel_Types_hpp__

#include"environment.hpp"
#include"kernel_Base.hpp"

using namespace Eigen;
using namespace std;


/*! Logarithm kernel */
class kernel_Logarithm: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};

/*! OneOverR2 kernel */
class kernel_OneOverR2: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};

/*! Gaussian kernel */
class kernel_Gaussian: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};

/*! Quadric kernel */
class kernel_Quadric: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};

/*! InverseQuadric kernel */
class kernel_InverseQuadric: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};

/*! ThinPlateSpline kernel */
class kernel_ThinPlateSpline: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1);
};


#endif //(__kerbel_Types_hpp__)
