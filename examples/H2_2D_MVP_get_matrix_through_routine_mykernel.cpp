//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//
/* Input type : Through matrix generating routine;
   Types of kernel: kernel defined by user */

#include"iostream"
#include<fstream>
#include"Eigen/Dense"
#include"cmath"
#include"ctime"
#include"./H2_2D_tree.hpp"
#include"kernel_types.hpp"

using namespace std;
using namespace Eigen;

void get_Location(unsigned long& N, VectorXd* location){
	N           =	5000;
	location[0]	=	VectorXd::Random(N);
	location[1]	=	VectorXd::Random(N);
}


void get_charges(const unsigned long N, unsigned& m, MatrixXd& Htranspose){
	m               =	10;
	Htranspose		=	MatrixXd::Random(N,m);
}

class mykernel: public kernel_base {
public:
    virtual double kernel_func(double R_square){
        //implement your own kernel here
        return 1.0 + R_square;
    }
};


int main(){
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
	unsigned long N;      // Number of charges;
	VectorXd location[2]; // Locations of the charges;
    get_Location(N,location);

	unsigned m;           // Number of sets of charges;
	MatrixXd Htranspose;  // All the different sets of charges;
    get_charges(N,m,Htranspose);


	cout << endl << "Number of charges:"    << N << endl;
	cout << endl << "Number of sets of charges:" << m << endl;
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /****************      Building fmm tree     **************/
    
	clock_t start_build	=	clock();
	unsigned short nchebnodes	=	6;                 // Number of Chebyshev nodes( >= 3) per dimension;
    H2_2D_tree Atree(nchebnodes, Htranspose, location);// Build the fmm tree;
    clock_t end_build	=	clock();
    
    double FMM_total_time_build	=	double(end_build-start_build)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(build tree) is: " << FMM_total_time_build << endl;

    /****************    Calculating potential   *************/
    
    clock_t start_A	=	clock();
	MatrixXd potentialA(N,m);
    /* Options of kernel: 
            LOGARITHM:          kernel_Logarithm
            ONEOVERR2:          kernel_OneOverR2
            GAUSSIAN:           kernel_Gaussian
            QUADRIC:            kernel_Quadric
            INVERSEQUADRIC:     kernel_InverseQuadric
            THINPLATESPLINE:    kernel_ThinPlateSpline 
     */
    mykernel A;
    A.calculatepotential(Atree,potentialA);
    clock_t end_A	=	clock();
    
    double FMM_total_time_A	=	double(end_A-start_A)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(calculate A) is: " << FMM_total_time_A << endl;


    /****     If you want to use more than one kernels    ****/
    
    /*clock_t start_B	=	clock();
    MatrixXd potentialB(N,m);
    kernel_Gaussian B;
    B.calculatepotential(Atree,potentialB);
    clock_t end_B	=	clock();
    
    double FMM_total_time_B	=	double(end_B-start_B)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(calculate B) is: " << FMM_total_time_B << endl;*/

    
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/

    cout << "Starting Exact computating..." << endl;
	clock_t start	=	clock();
	MatrixXd Q;
	A.kernel2D(N, location, N, location, Q);// Make sure the type of A here 
                                            // corresponds to the kernel used
                                            // to generate Q. 
	clock_t end	=	clock();

	double Exact_Assembly_time	=	double(end-start)/double(CLOCKS_PER_SEC);
	start	=	clock();
	MatrixXd potential_exact	=	Q*Htranspose;
	end	=	clock();
	double Exact_Computation_time	=	double(end-start)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for exact matrix vector product is: " << Exact_Assembly_time+Exact_Computation_time << endl;
    

	MatrixXd error              =	potentialA-potential_exact;
	double absolute_Error		=	(error).norm();
	double potential_Norm		=	(potential_exact).norm();
	cout << endl << "Relative difference in the solution is: " << absolute_Error/potential_Norm << endl;
    
    return 0;
}
