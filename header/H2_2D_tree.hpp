//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	H2_2D_tree.hpp
//
#ifndef __H2_2D_tree_hpp__
#define __H2_2D_tree_hpp__


#include"environment.hpp"
#include"H2_2D_node.hpp"

using namespace Eigen;
using namespace std;


const double PI	=	4.0*atan(1.0);

class H2_2D_tree{
    friend class kernel_base;
public:
	H2_2D_tree(const unsigned short nchebnodes, const MatrixXd& charge, const vector<Point>& location);
    unsigned short nchebnodes;	//	Number of chebyshev nodes along one direction;
    ~H2_2D_tree();
    
private:
    
    H2_2D_node* root;		//	Root node of tree;
    unsigned long N;		//	Total number of particles;
	unsigned m;             //	Number of sets of particles;
	unsigned short rank;		//	Rank of interaction: rank = chebnodes^2;
	VectorXd cnode;			//	Standardized chebyshev nodes in unit square [-1,1];
	MatrixXd Tnode;			//	First rank chebyshev polynomials evaluated at chebyshev nodes;

	unsigned maxlevels;		//	Counts maximum number of levels in tree;

	MatrixXd charge_tree;		//	All charges;
	MatrixXd location_tree[2];	//	All locations;

	MatrixXd R[4];			//	Translation matrix from cluster to its parent;
	unsigned short nlevels;		//	Total number of levels in tree;
    Point center;           // Center of tree;
    Point radius;           // Radius of tree;
//	Assigns children;
	void assignchildren(H2_2D_node*& node);
    
//  Build tree(assign siblings and cousins);
    void build_tree(H2_2D_node*& node);


//	Obtains charge to node when needed;
	void getcharge(H2_2D_node*& node);


//	Obtains standard Chebyshev nodes in interval [-1,1];
	void get_standard_Chebyshev_nodes(const unsigned short nchebnodes, VectorXd& cnode);

//	Obtains standard Chebyshev polynomials evaluated at given set of Points;
	void get_standard_Chebyshev_polynomials(const unsigned short nchebpoly, const unsigned long N, const VectorXd& x, MatrixXd& T);

//	Obtains center and radius of node location;
	void get_center_radius(const vector<Point>& location, Point& center, Point& radius);
    
    void max_And_Min_Coordinates(const vector<Point>& vec, double& max_x, double& max_y, double& min_x, double& min_y);

//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Points in children;
	void get_transfer_from_parent_to_children(const unsigned short nchebnodes, const vector<Point>& location, const Point& center, const Point& radius, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd& R);

//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Chebyshev nodes of children;
	void get_transfer_from_parent_cnode_to_children_cnode(const unsigned short nchebnodes, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd& Transfer);

//	Evaluates transfer from four children to parent;
	void gettransfer(const unsigned short nchebnodes, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd* R);

//	Evaluates 'nchebnodes' standardized chebyshev nodes in any interval;
	void getscaledchebnode(const unsigned short& nchebnodes, const VectorXd& cnode, const Point& center, const Point& radius, vector<Point>& chebnode);


//	Displays desired information;
	void display(H2_2D_node* node);

//	Assigns siblings;
	void assignsiblings(H2_2D_node*& node);

//	Assign cousins to children of node;
	void assigncousin(H2_2D_node*& node, unsigned short neighbor_number);

//	Deletes the tree;
	void deletetree(H2_2D_node*& node);
    
    
};


#endif //(__H2_2D_tree_hpp__)
