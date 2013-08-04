/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file H2_2D_Node.hpp
 Node in the FMM tree
 */
#ifndef __H2_2D_Node_hpp__
#define __H2_2D_Node_hpp__

#include"environment.hpp"
#include"Eigen/Dense"

using namespace Eigen;
using namespace std;

/*! Point in 2 dimesion */
struct Point {
    double x;
    double y;
    Point operator+ (Point const& rhs ) const;
    Point operator* (double rhs ) const;
};

/*! Node in FMM tree */
class H2_2D_Node{
public:
	H2_2D_Node* child[4];           /**<	Each node has 4 children; */
	H2_2D_Node* parent;             /**<	Each node has 1 parent;   */
	H2_2D_Node* neighbor[8];        /**<	Each node has atmost 8(=3^2-1^2) neighbors;*/
	H2_2D_Node* interaction[27];    /**<	Each node has atmost 27(=6^2-3^2) wellseparated interactions. */
	unsigned short nNeighbor;       /**<	Number of neighbors for each node;*/
	unsigned short nInteraction;	/**<	Number of nodes in the interaction list;*/
	unsigned short nLevel;          /**<	Level of the node. Root is at level 0; */
	unsigned short nodeNumber;      /**<	The number of the node for identification purposes. Each node is assigned an integer from 0 to 3. The left bottom is assigned a value 0, the right bottom is assigned a value 1, the left top takes the value 2 and the right top takes the value 3;*/

    Point center;                   /**<    Center of the node;*/
    Point radius;                   /**<    Radius of the node;*/
	unsigned long N;                /**<	Number of points inside the node;*/

	VectorXi index;                 /**<	Index of the charges in the original vector;*/
    vector<Point> location;         /**<	Location of the points inside the node;*/
    vector<Point> scaledCnode;      /**<    Scaled Chebyshev nodes along both the directions;*/

	MatrixXd charge;                /**<	Value of the charge for the points inside the node;*/
	MatrixXd potential;             /**<	Value of the potential for the points inside the node;*/

	MatrixXd nodeCharge;            /**<	The charges on the chebyshev nodes inside this cluster;*/
	MatrixXd nodePotential;         /**<	The potential on the chebyshev nodes inside this cluster;*/

	MatrixXd R;                     /**<	Transfer matrix from chebnode to cluster;*/

	bool isLeaf;                    /**<	TRUE-leaf node; FALSE-not a leaf node;*/
	bool isEmpty;                   /**<	Checks if there are charges inside a cluster;*/
	bool chargeComputed;            /**<	Checks if charge is already in the cluster;*/
	H2_2D_Node(unsigned short nLevel, unsigned short nodeNumber);

	~H2_2D_Node();
};

#endif //__H2_2D_Node_hpp_
