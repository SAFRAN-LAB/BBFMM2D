//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	H2_2D_node.cpp
//
#include "./H2_2D_node.hpp"

H2_2D_node::H2_2D_node(unsigned short nlevel, unsigned short nodenumber){
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
	for(unsigned short k=0; k<40; ++k){
		interaction[k]	=	NULL;
	}
	nneighbor	=	0;
	ninteraction	=	0;

	isleaf		=	false;
	isempty		=	false;
	charge_computed	=	false;

	this->nlevel	=	nlevel;
	this->nodenumber=	nodenumber;
}

H2_2D_node::~H2_2D_node(){
//	Delete the children
//	delete [] child;
}
