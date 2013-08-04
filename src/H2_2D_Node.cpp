/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	H2_2D_Node.cpp
*/
#include "H2_2D_Node.hpp"
Point Point::operator+ (Point const& rhs ) const {
    Point result;
    result.x    =   this->x+rhs.x;
    result.y    =   this->y+rhs.y;
    return result;
}

Point Point::operator* (double rhs ) const {
    Point result;
    result.x    =   this->x*rhs;
    result.y    =   this->y*rhs;
    return result;
}


H2_2D_Node::H2_2D_Node(unsigned short nLevel, unsigned short nodeNumber){
//	Set parent NULL
	parent	=	NULL;
//	Set children NULL
	for(unsigned short k=0; k<4; ++k){
		child[k]	=	NULL;
	}
//	Set neighbors NULL
	for(unsigned short k=0; k<8; ++k){
		neighbor[k]	=	NULL;
	}
//	Set interactions NULL
    for(unsigned short k=0; k<27; ++k){
		interaction[k]	=	NULL;
	}
	nNeighbor	=	0;
	nInteraction	=	0;

	isLeaf		=	false;
	isEmpty		=	false;
	chargeComputed	=	false;

	this->nLevel	=	nLevel;
	this->nodeNumber=	nodeNumber;
}

H2_2D_Node::~H2_2D_Node(){
    for (unsigned short k=0; k<4; ++k) {
        if (child[k]!=NULL) {
            delete child[k];
            child[k]    =   NULL;
        }
    }
}
