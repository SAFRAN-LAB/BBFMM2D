//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran</author>
//	
//	H2_2D_tree.hpp
//
#ifndef __H2_2D_tree_hpp__
#define __H2_2D_tree_hpp__

#include"iostream"
#include <sstream>
#include<fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <string>
#include <stdlib.h>
#include"cmath"
#include"Eigen/Dense"
#include"./H2_2D_node.hpp"
#include"./kernel.hpp"

using namespace Eigen;
using namespace std;


const double PI	=	4.0*atan(1.0);

class H2_2D_tree{
    friend class kernel_base;
public:
	H2_2D_tree(const unsigned short nchebnodes, const MatrixXd& charge, const VectorXd* location);
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
	double center[2], radius[2];	//	Center and radius of tree;

//	Assigns children;
	void assignchildren(H2_2D_node*& node);
    
//  Build tree(assign siblings and cousins);
    void build_tree(H2_2D_node*& node);


//	Obtains charge to node when needed;
	void getcharge(H2_2D_node*& node);


//	Obtains standard Chebyshev nodes in interval [-1,1];
	void get_standard_Chebyshev_nodes(const unsigned short nchebnodes, VectorXd& cnode);

//	Obtains standard Chebyshev polynomials evaluated at given set of points;
	void get_standard_Chebyshev_polynomials(const unsigned short nchebpoly, const unsigned long N, const VectorXd& x, MatrixXd& T);

//	Obtains center and radius of node location;
	void get_center_radius(const VectorXd* location, double* center, double* radius);

//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to points in children;
	void get_transfer_from_parent_to_children(const unsigned short nchebnodes, const VectorXd* location, const double* center, const double* radius, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd& R);

//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Chebyshev nodes of children;
	void get_transfer_from_parent_cnode_to_children_cnode(const unsigned short nchebnodes, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd& Transfer);

//	Evaluates transfer from four children to parent;
	void gettransfer(const unsigned short nchebnodes, const VectorXd& cnode, const MatrixXd& Tnode, MatrixXd* R);

//	Evaluates 'nchebnodes' standardized chebyshev nodes in any interval;
	void getscaledchebnode(const unsigned short& nchebnodes, const VectorXd& cnode, const double& center, const double& radius, VectorXd& chebnode);


//	Displays desired information;
	void display(H2_2D_node* node);

//	Assigns siblings;
	void assignsiblings(H2_2D_node*& node);

//	Assign cousins to children of node;
	void assigncousin(H2_2D_node*& node, unsigned short neighbor_number);

//	Deletes the tree;
	void deletetree(H2_2D_node*& node);
    
    
};


class kernel_base {
public:
    void calculatepotential(H2_2D_tree& tree, MatrixXd& potential);
    void calculatepotential(H2_2D_node*& node, MatrixXd& potential,H2_2D_tree& tree);
    void set_Tree_Potential_Zero(H2_2D_node* node);
    //	Obtains Chebyshev node potential from well separated clusters;
	void calculate_nodepotential_from_wellseparated_clusters(H2_2D_node*& node, unsigned short rank,unsigned short nchebnodes);
    //	Tranfers potential from node to final potential matrix when needed;
	void tranfer_potential_to_potential_tree(H2_2D_node*& node, MatrixXd& potential);
    //	Evaluate kernel at Chebyshev nodes;
	void kernelcheb2D(const unsigned short& M, const VectorXd* x, const unsigned short& N, const VectorXd* y, MatrixXd& K);
    //	Evaluate the kernel;
    void kernel2D(const unsigned long M, const VectorXd* x, const unsigned long N, const VectorXd* y, MatrixXd& Kernel);// new changed
    //	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;
	void transfer_nodepotential_to_child(H2_2D_node*& node,MatrixXd R[]);

    virtual double kernel_func(double) = 0;
};

// Child class of H2_2D_tree;
class kernel_Logarithm: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

class kernel_OneOverR2: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

class kernel_Gaussian: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

class kernel_Quadric: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

class kernel_InverseQuadric: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

class kernel_ThinPlateSpline: public kernel_base {
public:
    virtual double kernel_func(double R_square);
};

/* file format:
    location   Measurement operator
   (0.5, 0.5) (y,x,,)
   (1.5, 0.5)
    ...
 */
void read_Location_and_Measurement_operator (const string& filename, unsigned long& N, VectorXd* location, unsigned& m, MatrixXd& H);

void read_Measurement_operator(const string& s, unsigned long row, MatrixXd& H, unsigned m);

#endif //(__H2_2D_tree_hpp__)
