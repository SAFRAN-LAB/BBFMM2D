/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file H2_2D_Tree.hpp
 The FMM tree
 */
#ifndef __H2_2D_Tree_hpp__
#define __H2_2D_Tree_hpp__


#include"environment.hpp"
#include"H2_2D_Node.hpp"

using namespace Eigen;
using namespace std;


const double PI	=	4.0*atan(1.0);

class H2_2D_Tree{
    friend class kernel_Base;
public:
	H2_2D_Tree(const unsigned short nChebNodes, const MatrixXd& charge, const vector<Point>& location);
    unsigned short nChebNodes;	/**<	Number of chebyshev nodes along one direction;*/
    ~H2_2D_Tree();
    
private:
    
    H2_2D_Node* root;           /**<	Root node of tree;*/
    unsigned long N;            /**<	Total number of particles;*/
	unsigned m;                 /**<	Number of sets of particles;*/
	unsigned short rank;        /**<	Rank of interaction: rank = chebnodes^2;*/
	VectorXd cNode;             /**<	Standardized chebyshev nodes in unit square [-1,1];*/
	MatrixXd TNode;             /**<	First rank chebyshev polynomials evaluated at chebyshev nodes;*/

	unsigned maxLevels;         /**<	Counts maximum number of levels in tree;*/

	MatrixXd chargeTree;		/**<	All charges;*/
	vector<Point> locationTree;	/**<	All locations;*/

	MatrixXd R[4];              /**<	Translation matrix from cluster to its parent;*/
	unsigned short nLevels;		/**<	Total number of levels in tree;*/
    Point center;               /**<    Center of tree;*/
    Point radius;               /**<    Radius of tree;*/
/*!	Assigns children;*/
	void assign_Children(H2_2D_Node*& node);
    
/*!  Build tree(assign siblings and cousins);*/
    void build_Tree(H2_2D_Node*& node);


/*!	Obtains charge to node when needed;*/
	void get_Charge(H2_2D_Node*& node);


/*!	Obtains standard Chebyshev nodes in interval [-1,1];*/
	void get_Standard_Chebyshev_Nodes(const unsigned short nChebNodes, VectorXd& cNode);

/*!	Obtains standard Chebyshev polynomials evaluated at given set of Points;*/
	void get_Standard_Chebyshev_Polynomials(const unsigned short nChebPoly, const unsigned long N, const VectorXd& x, MatrixXd& T);

/*!	Obtains center and radius of node location;*/
	void get_Center_Radius(const vector<Point>& location, Point& center, Point& radius);
    
/*! Obtain the maximum and minimum of x coordinate and y coornidate of a vector of points; */
    void max_And_Min_Coordinates(const vector<Point>& vec, double& maxX, double& maxY, double& minX, double& minY);

/*!	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Points in children;*/
	void get_Transfer_From_Parent_To_Children(const unsigned short nChebNodes, const vector<Point>& location, const Point& center, const Point& radius, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd& R);

/*!	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Chebyshev nodes of children;*/
	void get_Transfer_From_Parent_CNode_To_Children_CNode(const unsigned short nChebNodes, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd& Transfer);

/*!	Evaluates transfer from four children to parent;*/
	void get_Transfer(const unsigned short nChebNodes, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd* R);

/*!	Evaluates 'nChebNodes' standardized chebyshev nodes in any interval;*/
	void get_Scaled_ChebNode(const unsigned short& nChebNodes, const VectorXd& cNode, const Point& center, const Point& radius, vector<Point>& chebNode);


/*!	Displays desired information;*/
	void display(H2_2D_Node* node);

/*!	Assigns siblings;*/
	void assign_Siblings(H2_2D_Node*& node);

/*!	Assign cousins to children of node;*/
	void assign_Cousin(H2_2D_Node*& node, unsigned short neighborNumber);
};


#endif //(__H2_2D_Tree_hpp__)
