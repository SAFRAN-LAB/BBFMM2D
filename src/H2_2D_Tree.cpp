/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	H2_2D_Tree.cpp
*/

#include"H2_2D_Tree.hpp"
H2_2D_Tree::H2_2D_Tree(const unsigned short nChebNodes, const MatrixXd& charge, const vector<Point>& location){
	this->nChebNodes    =	nChebNodes;
	this->rank          =	nChebNodes*nChebNodes;
	this->N             =	charge.rows();
	this->m             =	charge.cols();
	this->maxLevels		=	0;
	this->chargeTree	=	charge;
    this->locationTree  =   location;
        
    //	Get Chebyshev nodes
	cNode               =	VectorXd(nChebNodes);
	get_Standard_Chebyshev_Nodes(nChebNodes,cNode);
    //	Get Chebyshev polynomials evaluated at Chebyshev nodes
	TNode               =	MatrixXd(nChebNodes,nChebNodes);
	get_Standard_Chebyshev_Polynomials(nChebNodes,nChebNodes,cNode,TNode);
    //	Gets transfer matrices
	get_Transfer(nChebNodes,cNode,TNode,R);
	nLevels             =	0;
	get_Center_Radius(location, center, radius);
    //	Create root
	root                =	new H2_2D_Node(0, 0);
	root->nNeighbor		=	0;
	root->nInteraction	=	0;
	root->N             =	N;
    root->center        =   center;
    root->radius        =   radius;
    
	root->index.setLinSpaced(N,0,N-1);
    std::cout << "Assigning children..." << std::endl;
	assign_Children(root);
	std::cout << "Assigned children." << std::endl;
    build_Tree(root);
    std::cout << "Maximum levels is: " << this->maxLevels << std::endl;
}


/*void H2_2D_Tree::delete_Tree(H2_2D_Node*& node){
	if(node!=NULL){
		for(unsigned short k=0; k<4; ++k){
			delete_Tree(node->child[k]);
		}
		delete node;
		node = NULL;
	}
}*/

//	Assigns children;
void H2_2D_Tree::assign_Children(H2_2D_Node*& node){
	if(node->N==0){
		node->isLeaf	=	true;
		node->isEmpty	=	true;
	}
	else{
		node->potential		=	MatrixXd::Zero(node->N,m);
		node->nodePotential	=	MatrixXd::Zero(rank,m);
		node->nodeCharge	=	MatrixXd::Zero(rank,m);
        
        get_Scaled_ChebNode(nChebNodes, cNode, node->center, node->radius, node->scaledCnode);
        for(unsigned long k=0;k<node->N;++k){
            node->location.push_back(locationTree[node->index(k)]);
        }

		if(node->N<= (unsigned long) 4*rank){
			node->isLeaf        =	true;
			get_Transfer_From_Parent_To_Children(nChebNodes,node->location,node->center,node->radius,cNode,TNode,node->R);
			get_Charge(node);
			node->nodeCharge    =	node->nodeCharge+node->R.transpose()*node->charge;
			if(this->maxLevels<node->nLevel){
				this->maxLevels	=	node->nLevel;
			}
		}
		else{
			for(unsigned short k=0; k<4; ++k){
				node->child[k]              =	new H2_2D_Node(node->nLevel+1, k);
				node->child[k]->parent		=	node;
				node->child[k]->center.x	=	node->center.x+((k%2)-0.5)*node->radius.x;
				node->child[k]->center.y	=	node->center.y+((k/2)-0.5)*node->radius.y;
				node->child[k]->radius      =	node->radius*0.5;
                node->child[k]->N           =	0;
			}
			for(unsigned long k=0; k<node->N; ++k){
				if(locationTree[node->index(k)].x<node->center.x){
					if(locationTree[node->index(k)].y<node->center.y){
						node->child[0]->index.conservativeResize(node->child[0]->N+1);
						node->child[0]->index(node->child[0]->N)	=	node->index(k);
						++node->child[0]->N;
					}
					else{
						node->child[2]->index.conservativeResize(node->child[2]->N+1);
						node->child[2]->index(node->child[2]->N)	=	node->index(k);
						++node->child[2]->N;
					}
				}
				else{
					if(locationTree[node->index(k)].y<node->center.y){
						node->child[1]->index.conservativeResize(node->child[1]->N+1);
						node->child[1]->index(node->child[1]->N)	=	node->index(k);
						++node->child[1]->N;
					}
					else{
						node->child[3]->index.conservativeResize(node->child[3]->N+1);
						node->child[3]->index(node->child[3]->N)	=	node->index(k);
						++node->child[3]->N;
					}
				}
			}

			for(unsigned short k=0; k<4; ++k){
				assign_Children(node->child[k]);
				node->nodeCharge=	node->nodeCharge+R[k].transpose()*(node->child[k]->nodeCharge);
			}
		}
	}
}

//	Calculates potential;
void H2_2D_Tree::build_Tree(H2_2D_Node*& node){
	if(!node->isEmpty){
		if(!node->isLeaf){
			assign_Siblings(node);
			for(unsigned short k=0;k<8;++k){
				if(node->neighbor[k]!=NULL){
					if(!node->neighbor[k]->isLeaf){
                        assign_Cousin(node,k);
					}
				}
			}
            for(unsigned short k=0;k<4;++k) {
                build_Tree(node->child[k]);
            }
        }
    }
}



//	Obtains charge to node when needed;
void H2_2D_Tree::get_Charge(H2_2D_Node*& node){
	if(node->chargeComputed==true){
		return;
	}
	else{
		node->chargeComputed	=	true;
		node->charge		=	MatrixXd::Zero(node->N,m);
		for(unsigned long k=0;k<node->N;++k){
			node->charge.row(k)	=	chargeTree.row(node->index(k));
		}
	}
}


//	Obtains standard Chebyshev nodes in interval [-1,1];
void H2_2D_Tree::get_Standard_Chebyshev_Nodes(const unsigned short nChebNodes, VectorXd& cNode){
	cNode	=	VectorXd(nChebNodes);
	for(unsigned short k=0;k<nChebNodes;++k){
		cNode(k)	=	-cos((k+0.5)*PI/nChebNodes);
	}
}

//	Obtains standard Chebyshev polynomials evaluated at given set of Points;
void H2_2D_Tree::get_Standard_Chebyshev_Polynomials(const unsigned short nChebPoly, const unsigned long N, const VectorXd& x, MatrixXd& T){
	T	=	MatrixXd(N,nChebPoly);
	T.col(0)=	VectorXd::Ones(N);
	if(nChebPoly>1){
		T.col(1)	=	x;
		for(unsigned short k=2; k<nChebPoly; ++k){
			T.col(k)=	2.0*x.cwiseProduct(T.col(k-1))-T.col(k-2);
		}
	}
}

//	Obtains center and radius of node location;
void H2_2D_Tree::get_Center_Radius(const vector<Point>& location, Point& center, Point& radius){
	double maxX;
	double maxY;
	double minX;
	double minY;
    max_And_Min_Coordinates(location, maxX, maxY, minX, minY);
	center.x	=	0.5*(maxX+minX);
	center.y	=	0.5*(maxY+minY);
	radius.x	=	0.5*(maxX-minX);
	radius.y	=	0.5*(maxY-minY);
}



void H2_2D_Tree::max_And_Min_Coordinates(const vector<Point>& vec, double& maxX, double& maxY, double& minX, double& minY) {
    maxX   =   vec[0].x;
    maxY   =   vec[0].y;
    minX   =   maxX;
    minY   =   maxY;
    for (unsigned int i = 0; i < vec.size(); i++) {
        if (vec[i].x>maxX) {
            maxX   =   vec[i].x;
        }
        if (vec[i].y>maxY) {
            maxY   =   vec[i].y;
        }
        if (vec[i].x<minX) {
            minX   =   vec[i].x;
        }
        if (vec[i].y<minY) {
            minY   =   vec[i].y;
        }
    }
}


//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Points in children;
void H2_2D_Tree::get_Transfer_From_Parent_To_Children(const unsigned short nChebNodes, const vector<Point>& location, const Point& center, const Point& radius, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd& R){
	unsigned long N		=	location.size();
    VectorXd standlocation[2];
    standlocation[0].resize(N);
    standlocation[1].resize(N);
    for (unsigned long i = 0; i < N; i++) {
        standlocation[0](i) =   (location[i].x-center.x)/radius.x;
        standlocation[1](i) =   (location[i].y-center.y)/radius.y;
    }
    
	MatrixXd Transfer[2];
	for(unsigned short k=0;k<2;++k){
		get_Standard_Chebyshev_Polynomials(nChebNodes, N, standlocation[k], Transfer[k]);
		Transfer[k]	=	(2.0*Transfer[k]*TNode.transpose()-MatrixXd::Ones(N,nChebNodes))/nChebNodes;
	}
	unsigned short rank	=	nChebNodes*nChebNodes;
	R			=	MatrixXd(N,rank);
	for(unsigned short k=0;k<N;++k){
		for(unsigned short i=0; i<nChebNodes; ++i){
			for(unsigned short j=0;j<nChebNodes;++j){
				R(k,i+nChebNodes*j)=Transfer[0](k,i)*Transfer[1](k,j);
			}
		}
	}
}


//	Obtains interpolation operator, which interpolates information from Chebyshev nodes of parent to Chebyshev nodes of children;
void H2_2D_Tree::get_Transfer_From_Parent_CNode_To_Children_CNode(const unsigned short nChebNodes, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd& Transfer){
	VectorXd childcNode(2*nChebNodes);
	childcNode.segment(0,nChebNodes)		=	0.5*(cNode-VectorXd::Ones(nChebNodes));
	childcNode.segment(nChebNodes,nChebNodes)	=	0.5*(cNode+VectorXd::Ones(nChebNodes));
	get_Standard_Chebyshev_Polynomials(nChebNodes, 2*nChebNodes, childcNode, Transfer);
	Transfer	=	(2.0*Transfer*TNode.transpose()-MatrixXd::Ones(2*nChebNodes,nChebNodes))/nChebNodes;
}

//	Evaluates transfer from four children to parent;
void H2_2D_Tree::get_Transfer(const unsigned short nChebNodes, const VectorXd& cNode, const MatrixXd& TNode, MatrixXd* R){
	MatrixXd S(2*nChebNodes,nChebNodes);
	get_Transfer_From_Parent_CNode_To_Children_CNode(nChebNodes, cNode, TNode, S);
	MatrixXd Transfer[2];
	Transfer[0]	=	S.block(0,0,nChebNodes,nChebNodes);
	Transfer[1]	=	S.block(nChebNodes,0,nChebNodes,nChebNodes);
	unsigned short rank	=	nChebNodes*nChebNodes;
	for(unsigned short k=0; k<4; ++k){
		R[k]	=	MatrixXd(rank,rank);
	}
	for(unsigned short i=0;i<nChebNodes;++i){
		for(unsigned short j=0;j<nChebNodes;++j){
			for(unsigned short k=0;k<nChebNodes;++k){
				for(unsigned short l=0;l<nChebNodes;++l){
					R[0](i*nChebNodes+j,k*nChebNodes+l)	=	Transfer[0](i,k)*Transfer[0](j,l);
					R[1](i*nChebNodes+j,k*nChebNodes+l)	=	Transfer[0](i,k)*Transfer[1](j,l);
					R[2](i*nChebNodes+j,k*nChebNodes+l)	=	Transfer[1](i,k)*Transfer[0](j,l);
					R[3](i*nChebNodes+j,k*nChebNodes+l)	=	Transfer[1](i,k)*Transfer[1](j,l);
				}
			}
		}
	}
}

//	Evaluates 'nChebNodes' standardized chebyshev nodes in any interval;
void H2_2D_Tree::get_Scaled_ChebNode(const unsigned short& nChebNodes, const VectorXd& cNode, const Point& center, const Point& radius, vector<Point>& chebNode){
	for(unsigned short k=0;k<nChebNodes;++k){
		chebNode.push_back(center+radius*cNode(k));
	}
}


//	Displays desired information;
void H2_2D_Tree:: display(H2_2D_Node* node){
	std::cout << node->nLevel << ": " << node->N << std::endl << node->index << std::endl;
	if(!node->isLeaf){
		for(unsigned short k=0; k<4; ++k){
			display(node->child[k]);
		}
	}
}

//	Assign siblings to children of a the node
void H2_2D_Tree:: assign_Siblings(H2_2D_Node*& node){
//	Assign siblings to child[0]
	node->child[0]->neighbor[3]	=	node->child[1];
	node->child[0]->neighbor[5]	=	node->child[2];
	node->child[0]->neighbor[4]	=	node->child[3];

//	Assign siblings to child[1]
	node->child[1]->neighbor[7]	=	node->child[0];
	node->child[1]->neighbor[6]	=	node->child[2];
	node->child[1]->neighbor[5]	=	node->child[3];

//	Assign siblings to child[2]
	node->child[2]->neighbor[1]	=	node->child[0];
	node->child[2]->neighbor[2]	=	node->child[1];
	node->child[2]->neighbor[3]	=	node->child[3];

//	Assign siblings to child[3]
	node->child[3]->neighbor[0]	=	node->child[0];
	node->child[3]->neighbor[1]	=	node->child[1];
	node->child[3]->neighbor[7]	=	node->child[2];

	for(unsigned short k=0;k<4;++k){
		node->child[k]->nNeighbor+=3;
	}
}

//	Assign cousins to children of the node
void H2_2D_Tree:: assign_Cousin(H2_2D_Node*& node, unsigned short neighborNumber){
//	Assigning children of neighbor 0
	if(neighborNumber==0){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[0]->child[2];
		node->child[0]->neighbor[0]					=	node->neighbor[0]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[1];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[0];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[0]->child[3];

//	Update neighbor count.

		node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+1;
	}
//	Assigning children of neighbor 1
	else if(neighborNumber==1){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[1]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[1]->child[1];
		node->child[0]->neighbor[1]					=	node->neighbor[1]->child[2];
		node->child[0]->neighbor[2]					=	node->neighbor[1]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[1]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[1]->child[1];
		node->child[1]->neighbor[0]					=	node->neighbor[1]->child[2];
		node->child[1]->neighbor[1]					=	node->neighbor[1]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[0];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[1]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[1]->child[3];

//	Update neighbor count.

		node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+2;
		node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+2;
	}
//	Assigning children of neighbor 2
	else if (neighborNumber==2){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[2];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[1];
		node->child[1]->neighbor[2]					=	node->neighbor[2]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[0];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[2]->child[3];

//	Update neighbor count.

		node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+1;
	}
//	Assigning children of neighbor 3
	else if(neighborNumber==3){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[2];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->neighbor[3]					=	node->neighbor[3]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[3]->child[1];
		node->child[1]->neighbor[4]					=	node->neighbor[3]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[0];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->neighbor[2]					=	node->neighbor[3]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[3]->child[1];
		node->child[3]->neighbor[3]					=	node->neighbor[3]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[3]->child[3];

//	Update neighbor count.

		node->child[1]->nNeighbor					=	node->child[1]->nNeighbor+2;
		node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+2;
	}
//	Assigning children of neighbor 4
	else if(neighborNumber==4){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[2];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[1];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[0];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->neighbor[4]					=	node->neighbor[4]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[4]->child[3];

//	Update neighbor count.

		node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+1;
	}
//	Assigning children of neighbor 5
	else if(neighborNumber==5){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[2];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[1];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->neighbor[5]					=	node->neighbor[5]->child[0];
		node->child[2]->neighbor[4]					=	node->neighbor[5]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[5]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->neighbor[6]					=	node->neighbor[5]->child[0];
		node->child[3]->neighbor[5]					=	node->neighbor[5]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[5]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[5]->child[3];

//	Update neighbor count.

		node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+2;
		node->child[3]->nNeighbor					=	node->child[3]->nNeighbor+2;
}
//	Assigning children of neighbor 6
	else if (neighborNumber==6){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[0];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[2];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[1];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[0];
		node->child[2]->neighbor[6]					=	node->neighbor[6]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[2];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[6]->child[3];

//	Update neighbor count.

		node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+1;
	}
//	Assigning children of neighbor 7
	else if (neighborNumber==7){
//	Assigning the cousins to child0. One neighbor and three well-separated cousins.

		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[7]->child[0];
		node->child[0]->neighbor[7]					=	node->neighbor[7]->child[1];
		node->child[0]->interaction[node->child[0]->nInteraction++]	=	node->neighbor[7]->child[2];
		node->child[0]->neighbor[6]					=	node->neighbor[7]->child[3];

//	Assigning the cousins to child1. Four well-separated cousins.

		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[0];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[1];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[2];
		node->child[1]->interaction[node->child[1]->nInteraction++]	=	node->neighbor[7]->child[3];

//	Assigning the cousins to child2. Four well-separated cousins.

		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[7]->child[0];
		node->child[2]->neighbor[0]					=	node->neighbor[7]->child[1];
		node->child[2]->interaction[node->child[2]->nInteraction++]	=	node->neighbor[7]->child[2];
		node->child[2]->neighbor[7]					=	node->neighbor[7]->child[3];

//	Assigning the cousins to child3. Four well-separated cousins.

		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[0];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[1];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[2];
		node->child[3]->interaction[node->child[3]->nInteraction++]	=	node->neighbor[7]->child[3];

//	Update neighbor count.

		node->child[0]->nNeighbor					=	node->child[0]->nNeighbor+2;
		node->child[2]->nNeighbor					=	node->child[2]->nNeighbor+2;
	}
}

H2_2D_Tree::~H2_2D_Tree() {
    if(root!=NULL){
		delete root;
		root = NULL;
	}
}



