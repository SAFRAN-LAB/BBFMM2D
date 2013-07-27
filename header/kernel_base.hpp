//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	H2_2D_tree.hpp
//
#ifndef __kernel_base_hpp__
#define __kerbel_base_hpp__

#include"environment.hpp"
#include"H2_2D_tree.hpp"

using namespace Eigen;
using namespace std;


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
    void kernel2D(const unsigned long M, const VectorXd* x, const unsigned long N, const VectorXd* y, MatrixXd& Kernel);
    
    //	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;
	void transfer_nodepotential_to_child(H2_2D_node*& node,MatrixXd R[]);

    virtual double kernel_func(double) = 0;
};


#endif //(__kerbel_base_hpp__)
