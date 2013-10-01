/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file H2_2D_MVP_get_input_through_routine_mykernel.cpp
  Input type : Through matrix generating routine;
  Types of kernel: kernel defined by user 
*/

#include"BBFMM2D.hpp"

using namespace std;
using namespace Eigen;


/*! Get the location */
void get_Location(unsigned long& N, vector<Point>& location){
	N           =	5000;
    for (unsigned long i = 0; i < N; i++) {
        double x = rand() % 10 - 5;
        double y = rand() % 10 - 5;
        Point newPoint(x,y);
        location.push_back(newPoint);
    }
}

/*! Get charges */
void get_Charges(const unsigned long N, unsigned& m, double*& charges){
	m               =	10;
    charges =   new double[N*m];
    for (unsigned i = 0; i<m; i++) {
        for (unsigned long int j=0; j < N; j++) {
            charges[i*N+j]  =   rand() % 10 - 5;
        }
    }
}

/* Define user's own kernel */
class myKernel: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
        double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
        return exp(-pow(pow(rSquare,0.5)/900.0,0.5));
    }
};


int main(){
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
	unsigned long N;        // Number of charges;
	vector<Point> location; // Locations of the charges;
    get_Location(N,location);

	unsigned m;             // Number of sets of charges;
	double* charges;        // All the different sets of charges;
    get_Charges(N,m,charges);


	cout << endl << "Number of charges:"    << N << endl;
	cout << endl << "Number of sets of charges:" << m << endl;
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /****************      Building fmm tree     **************/
    
	clock_t startBuild	=	clock();
	unsigned short nchebnodes	=	6;    // Number of Chebyshev nodes( >= 3)
                                          // per dimension;
    cout << "Number of Chebyshev Nodes: " << nChebNodes << endl;

    H2_2D_Tree Atree(nchebnodes, charges, location, N, m);// Build the fmm tree;
    clock_t endBuild	=	clock();
    
    double FMMTotalTimeBuild	=	double(endBuild-startBuild)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(build tree) is: " << FMMTotalTimeBuild << endl;

    /****************    Calculating potential   *************/
    
    clock_t startA	=	clock();
	double* potentialA;
    potentialA = new double[N*m];
    /* Other options of kernel:
     LOGARITHM:          kernel_Logarithm
     ONEOVERR2:          kernel_OneOverR2
     GAUSSIAN:           kernel_Gaussian
     QUADRIC:            kernel_Quadric
     INVERSEQUADRIC:     kernel_InverseQuadric
     THINPLATESPLINE:    kernel_ThinPlateSpline
     */
    myKernel A;
    A.calculate_Potential(Atree,potentialA);
    clock_t endA	=	clock();
    
    double FMMTotalTimeA	=	double(endA-startA)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(calculating potential) is: " << FMMTotalTimeA << endl;
    
    
    /****     If you want to use more than one kernels    ****/
    
    /*clock_t startB	=	clock();
     double* potentialB;
     potentialB = new double[N*m];
     
     kernel_Gaussian B;
     B.calculate_Potential(Atree,potentialB);
     clock_t endB	=	clock();
     
     double FMMTotalTimeB	=	double(endB-startB)/double(CLOCKS_PER_SEC);
     cout << endl << "Total time taken for FMM(calculate B) is: " << FMMTotalTimeB << endl;*/
    
    
    
    /****           write data into binary file          ****/
    string outputfilename = "../output/potential.bin";
    write_Into_Binary_File(outputfilename, potentialA, N*m);
    
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/
    
    cout << "Starting Exact computating..." << endl;
	clock_t start	=	clock();
	MatrixXd Q;
	A.kernel_2D(N, location, N, location, Q);// Make sure the type of A here
                                             // corresponds to the kernel used
                                             // to generate Q.
	clock_t end	=	clock();
    
	double exactAssemblyTime	=	double(end-start)/double(CLOCKS_PER_SEC);
	start	=	clock();
    MatrixXd charges_           =   Map<MatrixXd>(charges, N, m);
	MatrixXd potentialExact     =	Q*charges_;
	end	=	clock();
	double exactComputationTime	=	double(end-start)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for exact matrix vector product is: " << exactAssemblyTime+exactComputationTime << endl;
    MatrixXd potentialA_   =   Map<MatrixXd>(potentialA, N, m);
    
	MatrixXd error              =	potentialA_-potentialExact;
	double absoluteError		=	(error).norm();
	double potentialNorm		=	(potentialExact).norm();
	cout << endl << "Relative difference in the solution is: " << absoluteError/potentialNorm << endl;
    
    /********       Clean Up        *******/
    delete []charges;
    delete []potentialA;
    //delete []potentialB;
    
    return 0;
}
