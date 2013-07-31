//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	kernel_Base.cpp
//

#include"kernel_Base.hpp"
void kernel_Base::calculate_Potential(H2_2D_Node*& node, MatrixXd& potential,H2_2D_Tree& tree){
    if(!node->isEmpty){
		if(node->isLeaf){
			MatrixXd tempK;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					kernel_2D(node->N , node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
                    //			Potential from neighbors
					tree.get_Charge(node->neighbor[k]);
					node->potential+=tempK*node->neighbor[k]->charge;
				}
			}
            //			Potential from Chebyshev nodes
			node->potential+=node->R*node->nodePotential;
            //			Self potential
			kernel_2D(node->N , node->location, node->N , node->location, tempK);
			node->potential+=tempK*node->charge;
            
			tranfer_Potential_To_Potential_Tree(node, potential);
		}
		else{
			bool computePotential	=	false;
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(node->neighbor[k]->isLeaf){
						MatrixXd tempK;
						kernel_2D(node->N, node->location, node->neighbor[k]->N, node->neighbor[k]->location, tempK);
						tree.get_Charge(node->neighbor[k]);
						node->potential+=tempK*node->neighbor[k]->charge;
						computePotential	=	true;
					}
				}
			}
			calculate_NodePotential_From_Wellseparated_Clusters(node,tree.rank,tree.nChebNodes);
			transfer_NodePotential_To_Child(node,tree.R);
			if(computePotential){
				tranfer_Potential_To_Potential_Tree(node, potential);
			}
			for(unsigned short k=0;k<4;++k){
				calculate_Potential(node->child[k], potential,tree);
			}
		}
	}
    
}

void kernel_Base::set_Tree_Potential_Zero(H2_2D_Node* node){
    if (node) {
        node->potential     =   MatrixXd::Zero(node->potential.rows(),node->potential.cols());
        node->nodePotential =   MatrixXd::Zero(node->nodePotential.rows(),node->potential.cols());
        for (unsigned short k=0; k<4; ++k) {
            set_Tree_Potential_Zero(node->child[k]);
        }
    }
}

//	Calculates potential;
void kernel_Base::calculate_Potential(H2_2D_Tree& tree, MatrixXd& potential){
    potential		=	MatrixXd(tree.N,tree.m);
    set_Tree_Potential_Zero(tree.root);
    std::cout << "Calculating potential..." << std::endl;
    calculate_Potential(tree.root,potential,tree);
    std::cout << "Calculated potential." << std::endl;
}

//	Obtains Chebyshev node potential from well separated clusters;
void kernel_Base::calculate_NodePotential_From_Wellseparated_Clusters(H2_2D_Node*& node, unsigned short rank,unsigned short nChebNodes){
	MatrixXd K(rank, rank);
	for(unsigned short k=0; k<4; ++k){
		for(unsigned short i=0; i<node->child[k]->nInteraction; ++i){
            if (!node->child[k]->interaction[i]->isEmpty && !node->child[k]->isEmpty) {
                kernel_Cheb_2D(nChebNodes,node->child[k]->scaledCnode,nChebNodes,node->child[k]->interaction[i]->scaledCnode,K);
                node->child[k]->nodePotential	=	node->child[k]->nodePotential+K*node->child[k]->interaction[i]->nodeCharge;

            }
		}
	}
}

//	Tranfers potential from node to final potential matrix when needed;
void kernel_Base::tranfer_Potential_To_Potential_Tree(H2_2D_Node*& node, MatrixXd& potential){
	for(unsigned long k=0; k<node->N; ++k){
		potential.row(node->index(k))+=node->potential.row(k);
	}
}

//	Evaluate kernel at Chebyshev nodes;
void kernel_Base::kernel_Cheb_2D(const unsigned short& M, const vector<Point>& xVec, const unsigned short& N, const vector<Point>& yVec, MatrixXd& K){
	vector<Point> xNew;
	vector<Point> yNew;
	K		=	MatrixXd(M*M,N*N);
	for(unsigned short j=0;j<M;++j){
		for(unsigned short i=0;i<M;++i){
            Point newPoint;
            newPoint.x     =   xVec[i].x;
            newPoint.y     =   xVec[j].y;
			xNew.push_back(newPoint);
		}
	}
	for(unsigned short j=0;j<N;++j){
		for(unsigned short i=0;i<N;++i){
            Point newPoint;
            newPoint.x     =   yVec[i].x;
            newPoint.y     =   yVec[j].y;
            yNew.push_back(newPoint);
		}
	}
	kernel_2D(M*M, xNew, N*N, yNew, K);
}

//	Tranfers potential from Chebyshev node of parent to Chebyshev node of children;
void kernel_Base::transfer_NodePotential_To_Child(H2_2D_Node*& node, MatrixXd R[]){
	for(unsigned short k=0;k<4;++k){
		node->child[k]->nodePotential	=	node->child[k]->nodePotential+R[k]*node->nodePotential;
	}
}


void kernel_Base::kernel_2D(const unsigned long M, const vector<Point>& x, const unsigned long N, const vector<Point>& y, MatrixXd& kernel) {
	kernel	=	MatrixXd::Zero(M,N);
	for(unsigned long i=0;i<M;++i){
		for(unsigned long j=0;j<N;++j){
            kernel(i,j) =   kernel_Func(x[i],y[j]);
        }
    }
}



