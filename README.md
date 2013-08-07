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


###DIRECTORIES AND FILES:


	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.  
	./src/			:	Source code in C++  
	./header/		:	Relevant header files  
	./exec/			:	Executables for BBFMM2D  
	./input/		:	The input file.  
	./README.md		:	This file  
	./License.md	:	License file  
	./Makefile		:	Makefile

###SETTING THINGS UP:

1. To run this package, you need to have **Eigen**.

2. Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page>

3. Create a directory named `Codes/` inside the main Eigen folder and copy the directory `BBFMM2D/` into the directory `Codes/`.

4. Open the Makefile, which is inside the folder BBFMM2D. Ensure that you have included the path to Eigen in the line containing `CFLAGS`. For instance, in the above setting, the path `"-I ./../../"` should be included in the Makefile.

5. Once you have this set up, you should be able to run the code. Check whether the code runs by performing the following action. Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make get_matrix_through_routine_standard_kernel
		cd exec/
		./H2_2D_MVP_get_matrix_through_routine_standard_kernel

The code should now run.

	
###CHANGING THE INPUTS:

The files you have control over are the files inside the directory ./examples/, read through the files both must be self explanatory for the most part.

1. If you want to generate matrix through your own routine, and use the standard kernels:

    Go to `/examples/`, open `"H2_2D_MVP_get_matrix_through_routine_standard_kernel.cpp"`.        
    * To generate matrix through routine:   
    
        Change `void get_Location(unsigned long& N, vector<Point>& location)` and `void get_Charges(const unsigned long N, unsigned& m, MatrixXd& Htranspose)`  
    * To use standard kernels:   
      
      Select the kernel type in `main()`.  
      Options of kernels:  
      
  	    	LOGARITHM:          kernel_Logarithm  
  	   		ONEOVERR2:          kernel_OneOverR2  
  			GAUSSIAN:           kernel_Gaussian  
  			QUADRIC:            kernel_Quadric  
    		INVERSEQUADRIC:     kernel_InverseQuadric  
    		THINPLATESPLINE:    kernel_ThinPlateSpline 
  

	
2. If you want to read matrix from text file, and use standard kernels:

    Go to the folder `/input/`, put your input file inside of this folder. 
     
    Go to the folder `/examples/`, open `"H2_2D_MVP_input_from_file_standard_kernel.cpp"`. 
    * To change input filename: 
     
      `string filenameInput = "../input/test_input.txt"`;  
    * To use standard kernels:  
    
      The same step as described in 1.


3. If you want to generate matrix through your own routine, and use your own kernel:

    Go to `/examples/`, open `"H2_2D_MVP_get_matrix_through_routine_myKernel.cpp"`.    
    * To define your own kernel: 
     
      Modify `class myKernel`. 
    * To generate your matrix:  
    
      The same step as described in 1.

4. If you want to read matrix from text file, and use your own kernel:  

	Go to `/examples/`, open `"H2_2D_MVP_input_from_file_myKernel.cpp"`.      
    * To define your own kernel:  
    
  	    Modify `class myKernel`. 
    * To change input filename: 
     
  	    The same step as described in 2. 
  	     
5. If you want to read matrix from binary file, and use standard kernel:
	
	Go to `/examples/`, open `"H2_2D_MVP_binary_file_standard_kernel.cpp"`.  
	* To change the input filename:  
	
		string filenameLocation     = "../input/test_Location.bin";
    	string filenameHtranspose   =   "../input/test_H.bin";
    	string filenameMetadata     =   "../input/metadata.txt";
    * To use standard kernels:  
    
    	The same step as described in 1.  

6. If you want to read matrix from binary file, and use your own kernel:  

	Go to `/examples/`, open `"H2_2D_MVP_binary_file_mykernel.cpp"`.
	* To change the input filename:  
	
	  	The same step as described in 5.  
	* To define your own kernel:  
	
	   Modify `class myKernel`.
	 
    	


###INPUT FILES  

Go to `/input/`, you should put your own input file in the input folder.

####TEXT FILES

The file format is described as follows:

The first row should be like this: 
 
`Number of charges, Number of sets of charges`

For example:

 	5000, 10

For the rest of the rows, it should start with locations, followed by a row in Htranspose matrix(elements should be separated using ','). If some element is 0, you can leave it as empty instead of 0. If all the elements in a row is 0, nothing need to be typed after the location.(spaces are allowed)

The row should look like this:  
  
`(location[0],location[1]) (elem1,elem2,elem3,elem4,…,elemn)`

For example:

	(-0.999984,-0.676221)   	(0.480685,0.869803,-0.188232,0.548587,-0.771039,0.73709,0.126494)    
	(0.869386,-0.5408)  
	(0.0655345,0.891162) (-0.193033,,-0.0287383,,-0.520512,0.33891,)  
	(0.342299,-0.246828) (0.0732668,,,,,,0.0951028)  
	(-0.984604,-0.44417) (,0.782447,-0.867924,0.485731,-0.729282,-0.481031,0.541473)  

####BINARY FILES

You should have 3 files:  

1. A binary file for H matrix: 

	Elements of H is stored in binary file row-wise.

2. A binary file for Location:

	Elements are stored this way(row-wise):
		
		loc0.x loc0.y  
		loc1.x loc1.y  
		…
3. A text file for metadata: 
 
   The file format is like this:  
   
   `Number of charges, Number of sets of charges`
   
   For example:

 		5000, 10
   
###RUNNING THE CODE:  

Here we give an example:  
If you want to use `"H2_2D_MVP_input_from_file_standard_kernel.cpp"`  

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make input_from_file_standard_kernel

2. Make sure you have changed or cleaned the .o files from previous compilation. To clean the irrelevant files, key in:

		make clean

3. To tar the file, key in

		make tar

4. Read through the makefile for other options.

To run other .cpp files:  

1) `H2_2D_MVP_get_matrix_through_routine_myKernel.cpp`     
   key in: 

      	make get_matrix_through_routine_myKernel   
   
2) `H2_2D_MVP_get_matrix_through_routine_standard_kernel.cpp`     
   key in:      

   		make get_matrix_through_routine_standard_kernel
   
3) `H2_2D_MVP_input_from_file_myKernel.cpp`    
   key in:  
   
   		make input_from_file_myKernel 
4) `H2_2D_MVP_binary_file_mykernel.cpp`  
   key in:
   
   		make binary_file_mykernel
5) `H2_2D_MVP_binary_file_standard_kernel.cpp`  
   key in:
   
   		make binary_file_standard_kernel
   
