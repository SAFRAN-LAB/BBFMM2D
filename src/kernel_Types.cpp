//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	kernel_Types.cpp
//

#include"kernel_Types.hpp"


double kernel_Logarithm::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    if (rSquare == 0){
        return 0;
    }
    else{
        return 0.5*log(rSquare);
    }
}

double kernel_OneOverR2::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    if (rSquare == 0){
        return 0;
    }
    else{
        return 1.0/rSquare;
    }
}

double kernel_Gaussian::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    return exp(-rSquare);
}

double kernel_Quadric::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    return 1.0+rSquare;
}

double kernel_InverseQuadric::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    return 1.0/(1.0+rSquare);
}

double kernel_ThinPlateSpline::kernel_Func(Point r0, Point r1){
    double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
    if (rSquare == 0){
        return 0;
    }
    else{
        return 0.5*rSquare*log(rSquare);
    }
}

