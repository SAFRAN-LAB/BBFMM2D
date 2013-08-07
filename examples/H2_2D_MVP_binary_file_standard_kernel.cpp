/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file H2_2D_MVP_binary_file_standard_kernel.cpp
   Input type : Binary file;
   Types of kernel: standard kernels 
 */

#include"environment.hpp"
#include"BBFMM2D.hpp"

using namespace std;
using namespace Eigen;


int main(){
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
	unsigned long N;        // Number of charges;
	vector<Point> location; // Locations of the charges;
	unsigned m;             // Number of sets of charges;
	MatrixXd Htranspose;    // All the different sets of charges;
    
    string filenameLocation     =   "../input/test_Location.bin";
    string filenameH            =   "../input/test_H.bin";
    string filenameMetadata     =   "../input/metadata.txt";
    
    read_Metadata_BBFMM2D (filenameMetadata, N, m);
    read_Location_Htranpose_binary (filenameLocation, N, location, filenameH, m, Htranspose);

	cout << endl << "Number of charges:"    << N << endl;
	cout << endl << "Number of sets of charges:" << m << endl;
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /****************      Building fmm tree     **************/
    
	clock_t startBuild	=	clock();
	unsigned short nChebNodes	=	8;                 // Number of Chebyshev nodes( >= 3)
                                                       // per dimension;
    H2_2D_Tree Atree(nChebNodes, Htranspose, location);// Build the fmm tree;
    clock_t endBuild	=	clock();
    
    double FMMTotalTimeBuild	=	double(endBuild-startBuild)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(build tree) is: " << FMMTotalTimeBuild << endl;
    
    /****************    Calculating potential   *************/
    
    clock_t startA	=	clock();
	MatrixXd potentialA(N,m);
    /* Options of kernel:
     LOGARITHM:          kernel_Logarithm
     ONEOVERR2:          kernel_OneOverR2
     GAUSSIAN:           kernel_Gaussian
     QUADRIC:            kernel_Quadric
     INVERSEQUADRIC:     kernel_InverseQuadric
     THINPLATESPLINE:    kernel_ThinPlateSpline
     */
    kernel_Gaussian A;
    A.calculate_Potential(Atree,potentialA);
    clock_t endA	=	clock();
    
    double FMMTotalTimeA	=	double(endA-startA)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for FMM(calculating potential) is: " << FMMTotalTimeA << endl;
    
    
    /****     If you want to use more than one kernels    ****/
    
    /*clock_t startB	=	clock();
     MatrixXd potentialB(N,m);
     kernel_Gaussian B;
     B.calculate_Potential(Atree,potentialB);
     clock_t endB	=	clock();
     
     double FMMTotalTimeB	=	double(endB-startB)/double(CLOCKS_PER_SEC);
     cout << endl << "Total time taken for FMM(calculate B) is: " << FMMTotalTimeB << endl;*/
    
   
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
	MatrixXd potentialExact	=	Q*Htranspose;
	end	=	clock();
	double exactComputationTime	=	double(end-start)/double(CLOCKS_PER_SEC);
	cout << endl << "Total time taken for exact matrix vector product is: " << exactAssemblyTime+exactComputationTime << endl;
    

	MatrixXd error              =	potentialA-potentialExact;
	double absoluteError		=	(error).norm();
	double potentialNorm		=	(potentialExact).norm();
	cout << endl << "Relative difference in the solution is: " << absoluteError/potentialNorm << endl;
    
    return 0;
}
