//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	kernel_base.cpp
//

#include"kernel_base.hpp"
void kernel_base::calculatepotential(H2_2D_node*& node, MatrixXd& potential,H2_2D_tree& tree){
    if(!node->isempty){
		if(node->isleaf){
			MatrixXd temp_K;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					kernel2D(node->N , node->location, node->neighbor[k]->N, node->neighbor[k]->location, temp_K);
                    //			Potential from neighbors
					tree.getcharge(node->neighbor[k]);
					node->potential+=temp_K*node->neighbor[k]->charge;
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodepotential;
            //			Self potential
			kernel2D(node->N , node->location, node->N , node->location, temp_K);
			node->potential+=temp_K*node->charge;
            
			tranfer_potential_to_potential_tree(node, potential);
		}
		else{
			bool compute_potential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(node->neighbor[k]->isleaf){
						MatrixXd temp_K;
						kernel2D(node->N, node->location, node->neighbor[k]->N, node->neighbor[k]->location, temp_K);
						tree.getcharge(node->neighbor[k]);
						node->potential+=temp_K*node->neighbor[k]->charge;
						compute_potential	=	true;
					}
				}
			}
			calculate_nodepotential_from_wellseparated_clusters(node,tree.rank,tree.nchebnodes);
			transfer_nodepotential_to_child(node,tree.R);
			if(compute_potential){
				tranfer_potential_to_potential_tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculatepotential(node->child[k], potential,tree);
			}
		}
	}
    
}

void kernel_base::set_Tree_Potential_Zero(H2_2D_node* node){
    if (node) {
        node->potential     =   MatrixXd::Zero(node->potential.rows(),node->potential.cols());
        node->nodepotential =   MatrixXd::Zero(node->nodepotential.rows(),node->potential.cols());
        for (unsigned short k=0; k<4; ++k) {
            set_Tree_Potential_Zero(node->child[k]);
        }
    }
}

//	Calculates potential;
void kernel_base::calculatepotential(H2_2D_tree& tree, MatrixXd& potential){
    potential		=	MatrixXd(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
    std::cout << "Calculating potential..." << std::endl;
    calculatepotential(tree.root,potential,tree);
    std::cout << "Calculated potential." << std::endl;
}

//	Obtains Chebyshev node potential from well separated clusters;
void kernel_base::calculate_nodepotential_from_wellseparated_clusters(H2_2D_node*& node, unsigned short rank,unsigned short nchebnodes){
	MatrixXd K(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		for(unsigned short i=0; i<node->child[k]->ninteraction; ++i){
            if (!node->child[k]->interaction[i]->isempty && !node->child[k]->isempty) {
                kernelcheb2D(nchebnodes,node->child[k]->scaledcnode,nchebnodes,node->child[k]->interaction[i]->scaledcnode,K);
                node->child[k]->nodepotential	=	node->child[k]->nodepotential+K*node->child[k]->interaction[i]->nodecharge;

            }
		}
	}
}

//	Tranfers potential from node to final potential matrix when needed;
void kernel_base::tranfer_potential_to_potential_tree(H2_2D_node*& node, MatrixXd& potential){
	for(unsigned long k=0; k<node->N; ++k){
		potential.row(node->index(k))+=node->potential.row(k);
	}
}

//	Evaluate kernel at Chebyshev nodes;
void kernel_base::kernelcheb2D(const unsigned short& M, const vector<Point>& x_Vec, const unsigned short& N, const vector<Point>& y_Vec, MatrixXd& K){
	vector<Point> xnew;
	vector<Point> ynew;
	K		=	MatrixXd(M*M,N*N);
	for(unsigned short j=0;j<M;++j){
		for(unsigned short i=0;i<M;++i){
            Point new_Point;
            new_Point.x     =   x_Vec[i].x;
            new_Point.y     =   x_Vec[j].y;
			xnew.push_back(new_Point);
		}
	}
	for(unsigned short j=0;j<N;++j){
		for(unsigned short i=0;i<N;++i){
            Point new_Point;
            new_Point.x     =   y_Vec[i].x;
            new_Point.y     =   y_Vec[j].y;
            ynew.push_back(new_Point);
		}
	}
	kernel2D(M*M, xnew, N*N, ynew, K);
}

//	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;
void kernel_base::transfer_nodepotential_to_child(H2_2D_node*& node, MatrixXd R[]){
	for(unsigned short k=0;k<4;++k){
		node->child[k]->nodepotential	=	node->child[k]->nodepotential+R[k]*node->nodepotential;
	}
}


void kernel_base::kernel2D(const unsigned long M, const vector<Point>& x, const unsigned long N, const vector<Point>& y, MatrixXd& Kernel) {
	Kernel	=	MatrixXd::Zero(M,N);
	for(unsigned long i=0;i<M;++i){
		for(unsigned long j=0;j<N;++j){
            Kernel(i,j) =   kernel_func(x[i],y[j]);
        }
    }
}



