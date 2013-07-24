//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran</author>
//	
//	H2_2D_node.hpp
//
#ifndef __H2_2D_node_hpp__
#define __H2_2D_node_hpp__

#include"Eigen/Dense"

using namespace Eigen;

class H2_2D_node{
public:
	H2_2D_node* child[4];		//	Each node has 4 children;
	H2_2D_node* parent;		//	Each node has 1 parent;
	H2_2D_node* neighbor[8];	//	Each node has atmost 8(=3^2-1^2) neighbors;
	H2_2D_node* interaction[27];	//	Each node has atmost 27(=6^2-3^2) wellseparated interactions. 40(=7^2-3^2) is the number of distinct well-separated interactions possible;

	unsigned short nneighbor;	//	Number of neighbors for each node;
	unsigned short ninteraction;	//	Number of nodes in the interaction list;
	unsigned short nlevel;		//	Level of the node. Root is at level 0;
	unsigned short nodenumber;	//	The number of the node for identification purposes. Each node is assigned an integer from 0 to 3. The left bottom is assigned a value 0, the right bottom is assigned a value 1, the left top takes the value 2 and the right top takes the value 3;

	double center[2], radius[2];	//	Center and radius of the node;

	unsigned long N;		//	Number of points inside the node;

	VectorXi index;			//	Index of the charges in the original vector;
	VectorXd location[2];		//	Location of the points inside the node;
	VectorXd scaledcnode[2];	//  Scaled Chebyshev nodes along both the directions;

	MatrixXd charge;		//	Value of the charge for the points inside the node;
	MatrixXd potential;		//	Value of the potential for the points inside the node;

	MatrixXd nodecharge;		//	The charges on the chebyshev nodes inside this cluster;
	MatrixXd nodepotential;		//	The potential on the chebyshev nodes inside this cluster;

	MatrixXd R;			//	Transfer matrix from chebnode to cluster;

	bool isleaf;			//	TRUE-leaf node; FALSE-not a leaf node;
	bool isempty;			//	Checks if there are charges inside a cluster;
	bool charge_computed;		//	Checks if charge is already in the cluster;
	H2_2D_node(unsigned short nlevel, unsigned short nodenumber);

	~H2_2D_node();
};

#endif //__H2_2D_node_hpp_
