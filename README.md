
<style TYPE="text/css">
code.has-jax {font: inherit; font-size: 100%; background: inherit; border: inherit;}
</style>
<script type="text/x-mathjax-config">

MathJax.Hub.Config({
    tex2jax: {
        inlineMath: [['$','$'], ['\\(','\\)']],
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'] // removed 'code' entry
    }
});
MathJax.Hub.Queue(function() {
    var all = MathJax.Hub.getAllJax(), i;
    for(i = 0; i < all.length; i += 1) {
        all[i].SourceElement().parentNode.className += ' has-jax';
    }
});
</script>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
#BBFMM2D  

This is the first public release of the BBFMM2D library.  
Date: July 24th, 2013

**Version 3.1 - First external release.**

%% Copyleft 2013: Sivaram Ambikasaran, Ruoxi Wang and Eric Darve         
%% Developed by Sivaram Ambikasaran, Ruoxi Wang         
%% Contact: <siva.1985@gmail.com> (Sivaram), <ruoxi@stanford.edu> (Ruoxi)     
%%    
%% This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license.      
%% The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not %% distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>  

###1. INTRODUCTION
BBFMM2D is an open source package of the <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">Black-box Fast Multipole Method</a> in 2 dimensions.   
The Black-box Fast Multipole Method is an $\mathcal{O}(N)$ fast multipole method, which is a technique to calculate sums of the form $f(x_i) = \displaystyle \sum_{j=1}^N K(x_i,y_j) \sigma_j, \,\,\, \forall i \in\{1,2,\ldots,N\}$.  
BBFMM2D provides an $\mathcal{O}(Nlog(N))$ or  $\mathcal{O}(N)$ solution to matrix-matrix product $Q \times H$, where $Q$ is a covariance matrix of size $N\times N$ with a kernel, and $N$ is the number of unknown values at points $(x,y)$ in a 2D domain. $H$ is a $N \times m$ matrix with $N >> m$. 
This fast multipole method is applicable for a wide range of non-oscillatory kernels, which can be evaluated for any pair of points. This Black-Box Fast Multipole Method relies on Chebyshev interplation to construct low-rank approximations for well-separated clusters. The main advantage of this fast multipole method is that it requires minimal pre-computation time and is applicable for a wide range of kernels.

###2. DIRECTORIES AND FILES


	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.  
	./src/			:	Source code in C++  
	./header/		:	Relevant header files  
	./exec/			:	Executables for BBFMM2D  
	./input/		:	The input file.  
	./README.md		:	This file  
	./License.md	:	License file  
	./Makefile		:	Makefile
	
###3. TUTORIAL
####3.1 To Get Started  
1. To use BBFMM2D, you need to have <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a> library.

2. Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page>.

3. Create a directory named `Codes/` inside the main Eigen folder and copy the directory `BBFMM2D/` into the directory `Codes/`.

4. Open the Makefile, which is inside the folder BBFMM2D. Ensure that you have included the path to Eigen in the line containing `CFLAGS`. For instance, in the above setting, the path `"-I ./../../"` should be included in the Makefile.  

5. To check whether things are set up correctly, you can perform the following action: Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make binary_file_mykernel
		cd exec/
		./H2_2D_MVP_binary_file_mykernel

####3.2 FMM with standard kernel
#####3.2.1 Basic usage

The basic usage of BBFMM2D is as follows: 

	#include"BBFMM2D.hpp"  
	...
	{
	unsigned long N;        // Number of charges;
	unsigned m;             // Number of sets of charges;
	vector<Point> location; // Locations of the charges;
	double* charges;        // All the different sets of charges;
	unsigned short nChebNodes	=	9; // Number of Chebyshev nodes per dimension; 
	…
	H2_2D_Tree Atree(nChebNodes, charges, location, N, m);// Build the fmm tree;
	
	/* The following can be repeated with different kernels */
	double* potential;
    potential = new double[N*m];
    
    kernel_Gaussian A;
    A.calculate_Potential(Atree,potential);
    ...
    }
    
This example first build a fmm tree with this line:  
`H2_2D_Tree Atree(nChebNodes, charges, location, N, m);`  
where H2_2D_Tree is a class of fmm tree, the constructor takes 5 arguments:  

* nChebNodes(unsigned short):   
	Number of Chebyshev nodes per dimension. It should take value $\ge$ 3, and we recommend to take value from 3 to 10. (Larger number of Chebyshev nodes would give better result but with much more time)
* charges(double*):   
	All the different sets of charges. This pointer should point to an array with size $N \times m$, and the data should be stored in column-wise. ( i.e. first set of charges, followd by second set of charges, etc)
* location(vector<Point>):  
	Locations of the charges in $2D$ domain. Here Point is a structure type with $x$ and $y$ coordinate defined.  
* N(unsigned long):  
	Number of charges.  
* m(unsigned):  
	Number of sets of charges.  
	
Once the tree is created, you can compute the matrix-matrix product with as many kernels as you want.(Please see **3.2.3**) The code shows an example using Gaussian kernel:  

	kernel_Gaussian A;
    A.calculate_Potential(Atree,potential);
The result is computed via `calculate_Potential()`, which is a method of class `kernel_Gaussian`. The first argument of `calculate_Potential()` is the fmm tree that we just created; and the second argument potential is a pointer to the result, and the result is stored column-wise in `potential`.  

#####3.2.2 Options of kernels

We have provided several standard kernels:  
To demonstrate the deltails of each kernel, we denote the element of the covariance matrix $Q_{ij} = k(r )$, where $r$ is the distance between point i and point j, is described by a kernel function.

Options of kernels:  

* LOGARITHM kernel:           
	usage: kernel_Logarithm  
	kernel function: $k(r) = 0.5 \times log(r^2)\, (r\neq 0);\, k(r )= 0 \,(r=0).$  
	
	
* ONEOVERR2 kernel:  
	usage: kernel_OneOverR  
	kernel function:  $k(r ) = 1 / r^2 \,(r \neq 0);\, k(r )= 0 \,(r=0)$.  
	
* GAUSSIAN kernel:  
	usage: kernel_Gaussian  
	kernel function: $k(r ) = exp(-r^2)$.  
	
* QUADRIC kernel:  
	usage: kernel_Quadric  
	kernel function: $ k(r ) = 1 + r^2$.  

* INVERSEQUADRIC kernel:  
	usage: kernel_InverseQuadric  
	kernel function: $k(r ) = 1 / (1+r^2)$.
	
* THINPLATESPLINE kernel:  
	usage:  kernel_ThinPlateSpline  
	kernel function: $k(r ) =  0.5 \times r^2 \times log(r^2 )\, (r \neq 0);\, k(r )=0\,(r=0)$.
    		
If you want to define your own kernel, please see **3.3**.

#####3.2.3 Usage of multiple kernels

	#include"BBFMM2D.hpp"  
	...
	{
	…
	H2_2D_Tree Atree(nChebNodes, charges, location, N, m);// Build the fmm tree;
	
	/* The following can be repeated with different kernels */
	...
    kernel_Gaussian A;
    A.calculate_Potential(Atree,potentialA);
    ...   
    kernel_Quadric B;
    A.calculate_Potential(Atree,potentialB);
	...
    }
The basic usage is already domonstrated in **3.2.1**. Once you have built the fmm tree, you can use different kernels to compute the matrix-matrix multiplication without rebuilding the tree. You can choose kernel type from standard kernels given by us ( see **3.2.2** ), or you can define your own kernel ( see **3.3** )
	 
####3.3 FMM with user defined kernel

#####3.3.1 Basic usage
The basic usage is almost the same as **3.2** except that you have to define your own routine of computing kernel. One example code is as follows:  

	#include"BBFMM2D.hpp"  
	class myKernel: public kernel_Base {
	public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
        double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
        return exp(-pow(pow(rSquare,0.5)/900.0;8,0.5));
    	}
	};
	...
	{
	…
	H2_2D_Tree Atree(nChebNodes, charges, location, N, m);// Build the fmm tree;
	
	/* The following can be repeated with different kernels */
	...
    myKernel A;
    A.calculate_Potential(Atree,potential);
    ...
    }

You can define your own kernel inside `kernel_Func(Point r0, Point r1)`, it takes two Points as input and returns a double value ( $Q_{ij}$ ). 

#####3.3.2 Usage of multiple kernels

You can also define and use multiple kernels in one file, but make sure to have different class names. 
e.g.  

	class myKernelA: public kernel_Base{
		public:
    	virtual double kernel_Func(Point r0, Point r1) {
    	...
    	}
	}
	class myKernelB: public kernel_Base{
		public:
    	virtual double kernel_Func(Point r0, Point r1) {
    	...
    	}
	}
	…
	
	{
	…
	H2_2D_Tree Atree(nChebNodes, charges, location, N, m);// Build the fmm tree;
	
	/* The following can be repeated with different kernels */
	...
    myKernelA A;
    A.calculate_Potential(Atree,potentialA);
    ...
    myKernelB B;
    B.calculate_Potential(Atree,potentialB);
    …  
    kernel_Gaussian C;
    C.calculate_Potential(Atree,potentialC);
    ...
    }

###4. ROUTINES FOR INPUTING AND OUTPUTING DATA:  
We have provided several routines for reading data from text file and binary file, and writing data into binary file.	

####4.1 Reading meta data from text file
	
	void read_Metadata_BBFMM2D (const string& filenameMetadata, unsigned long& N, unsigned& m);
The first argument, filenameMetadata is the filename for your meta data (the number of locations, number of sets of charges). The number of locations is stored in N, and number of sets of charges is stored in m.

**File format:**  
 
`Number of charges, Number of sets of charges`

For example:

 	5000, 10
 	
####4.2 Reading from text file 
The prototype of function to read input from text file is:  

	void read_Location_Charges (const string& filename, unsigned long N, vector<Point>& location, unsigned m, double*& charges);

The first argument is the filename of your text file, the second argument N and forth argument m are the number of locations and number of sets of charges respectively.
This function stores location in `location` and stores charges column-wise in `charges`.

**File format:**  
For each row, it should start with locations, and followed by a row in $H^T$ ( if we do $QH^T$ multiplication, $H^T$ is the R.H.S.). Here note that elements should be separated using ','. If some element is 0, you can leave it as empty instead of 0. If all the elements in a row is 0, nothing need to be typed after the location.(spaces are allowed)

The row should look like this:  
  
`(location[0],location[1]) (elem1,elem2,elem3,elem4,…,elemn)`

For example:

	(-0.999984,-0.676221)   	(0.480685,0.869803,-0.188232,0.548587,-0.771039,0.73709,0.126494)    
	(0.869386,-0.5408)  
	(0.0655345,0.891162) (-0.193033,,-0.0287383,,-0.520512,0.33891,)  
	(0.342299,-0.246828) (0.0732668,,,,,,0.0951028)  
	(-0.984604,-0.44417) (,0.782447,-0.867924,0.485731,-0.729282,-0.481031,0.541473)  


####4.3 Reading from binary file  

	void read_Location_Charges_binary(const string& filenameLocation, unsigned long N, vector<Point>& location, const string& filenameHtranspose,unsigned m, double* charges);

The first argument filenameLocation and the forth argument filenameHtranspose are binary file names for location and $H^T$ respectively. N is the number of locations and m is the number of sets of charges. The data of locations is stored in `location` and the data of charges is stored in `charges` column-wise.  

**File format:** 
 
1. Binary file for charges: 

	Elements of charges is stored this way:  
	It should be stored column-wise, i.e.   
	first set of charges, followed by second set of charges, etc.

2. Binary file for Location:

	Elements are stored this way(row-wise):
		
		loc0.x loc0.y  
		loc1.x loc1.y  
		…

####4.4 Writing into binary file  
	
	void write_Into_Binary_File(const string& filename, double* outdata, int numOfElems);  
This first argument is the filename for your output data. The second argument is a pointer to the output data, and the last argument is the number of elements in the array of your output data.  


###5. EXAMPLES

We have provided several examples for BBFMM2D. Go to examples/, read through the files both must be self explanatory for the most part.
You can use our examples with your own input.
####5.1 Chage input of examples

1. If you want to generate input through your own routine, and use the standard kernels:

    Go to `/examples/`, open `"H2_2D_MVP_get_input_through_routine_standard_kernel.cpp"`.        
    * To generate input through routine:   
    
        Change `void get_Location()` and `void get_Charges()`  
    * To use standard kernels:   
      
      Choose the kernel type in `main()`, options of kernels are in **3.2.2**
  

	
2. If you want to read input from text file, and use standard kernels:

    Go to the folder `/input/`, put your input file inside of this folder. 
     
    Go to the folder `/examples/`, open `"H2_2D_MVP_textfile_standard_kernel.cpp"`. 
    * To change input filename: 
     
      `string filenameInput = "../input/test_input.txt"`;  
    * To use standard kernels:  
    
      The same step as described in 1.


3. If you want to generate input through your own routine, and use your own kernel:

    Go to `/examples/`, open `"H2_2D_MVP_get_input_through_routine_myKernel.cpp"`.    
    * To define your own kernel: 
     
      Modify `class myKernel`. 
    * To generate your input:  
    
      The same step as described in 1.

4. If you want to read input from text file, and use your own kernel:  

	Go to `/examples/`, open `"H2_2D_MVP_textfile_myKernel.cpp"`.      
    * To define your own kernel:  
    
  	    Modify `class myKernel`. 
    * To change input filename: 
     
  	    The same step as described in 2. 
  	     
5. If you want to read input from binary file, and use standard kernel:
	
	Go to `/examples/`, open `"H2_2D_MVP_binary_file_standard_kernel.cpp"`. 			     
    * To change input filename:  
      
		`string filenameLocation     =   "../input/test_Location.bin";`  
    	`string filenameCharges			  =   "../input/test_Charges.bin";`  
    	`string filenameMetadata     =   "../input/metadata.txt";`  
    change them into your input filenames. 
    * To use standard kernels: 
     
  	    The same step as described in 1. ;  

6. If you want to read input from binary file, and use your own kernel:  

	Go to `/examples/`, open `"H2_2D_MVP_binary_file_mykernel.cpp"`.
	* To change the input filename:  
	
	  	The same step as described in 5.  
	* To define your own kernel:  
	
	   Modify `class myKernel`.
	 
When using our examples, make sure that the input file format are the same as described in  **4.**  	

####5.2 Run examples  

Here we give an example:  
If you want to use `"H2_2D_MVP_textfile_standard_kernel.cpp"`  

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make textfile_standard_kernel

2. Make sure you have changed or cleaned the .o files from previous compilation. To clean the irrelevant files, key in:

		make clean

3. To tar the file, key in

		make tar

4. Read through the makefile for other options.

To run other .cpp files:  

1) `H2_2D_MVP_get_input_through_routine_myKernel.cpp`     
   key in: 

      	make get_input_through_routine_myKernel   
   
2) `H2_2D_MVP_get_input_through_routine_standard_kernel.cpp`     
   key in:      

   		make get_input_through_routine_standard_kernel
   
3) `H2_2D_MVP_textfile_myKernel.cpp`    
   key in:  
   
   		make textfile_myKernel 
4) `H2_2D_MVP_binary_file_mykernel.cpp`  
   key in:
   
   		make binary_file_mykernel
5) `H2_2D_MVP_binary_file_standard_kernel.cpp`  
   key in:
   
   		make binary_file_standard_kernel
   		
