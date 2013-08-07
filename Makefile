CC=g++
CFLAGS=-c -Wall -DNDEBUG -O3 -ffast-math -ffinite-math-only -I header/ -I ./../../
LDFLAGS=

SOURCES=src/H2_2D_Tree.cpp src/H2_2D_Node.cpp src/read_Location_H.cpp src/kernel_Base.cpp src/kernel_Types.cpp

SOURCESBINARY=src/read_metadata_BBFMM2D.cpp src/read_Location_Htranpose_binary.cpp

SOURCESA=examples/H2_2D_MVP_input_from_file_standard_kernel.cpp
SOURCESB=examples/H2_2D_MVP_get_matrix_through_routine_standard_kernel.cpp 
SOURCESC=examples/H2_2D_MVP_input_from_file_mykernel.cpp
SOURCESD=examples/H2_2D_MVP_get_matrix_through_routine_mykernel.cpp
SOURCESE=examples/H2_2D_MVP_binary_file_standard_kernel.cpp
SOURCESF=examples/H2_2D_MVP_binary_file_mykernel.cpp


OBJECTSA=$(SOURCES:.cpp=.o) $(SOURCESA:.cpp=.o) 
OBJECTSB=$(SOURCES:.cpp=.o) $(SOURCESB:.cpp=.o) 
OBJECTSC=$(SOURCES:.cpp=.o) $(SOURCESC:.cpp=.o)
OBJECTSD=$(SOURCES:.cpp=.o) $(SOURCESD:.cpp=.o)
OBJECTSE=$(SOURCES:.cpp=.o) $(SOURCESE:.cpp=.o) $(SOURCESBINARY:.cpp=.o)
OBJECTSF=$(SOURCES:.cpp=.o) $(SOURCESF:.cpp=.o) $(SOURCESBINARY:.cpp=.o)


EXECUTABLEA=./exec/H2_2D_MVP_input_from_file_standard_kernel
EXECUTABLEB=./exec/H2_2D_MVP_get_matrix_through_routine_standard_kernel
EXECUTABLEC=./exec/H2_2D_MVP_input_from_file_mykernel
EXECUTABLED=./exec/H2_2D_MVP_get_matrix_through_routine_mykernel
EXECUTABLEE=./exec/H2_2D_MVP_binary_file_standard_kernel
EXECUTABLEF=./exec/H2_2D_MVP_binary_file_mykernel


input_from_file_standard_kernel: $(SOURCES) $(SOURCESA) $(EXECUTABLEA)

$(EXECUTABLEA): $(OBJECTSA)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSA) -o $@


get_matrix_through_routine_standard_kernel: $(SOURCES) $(SOURCESB) $(EXECUTABLEB)

$(EXECUTABLEB): $(OBJECTSB)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSB) -o $@


input_from_file_mykernel: $(SOURCES) $(SOURCESC) $(EXECUTABLEC)

$(EXECUTABLEC): $(OBJECTSC)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSC) -o $@


get_matrix_through_routine_mykernel: $(SOURCES) $(SOURCESD) $(EXECUTABLED)

$(EXECUTABLED): $(OBJECTSD)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSD) -o $@

binary_file_standard_kernel: $(SOURCES) $(SOURCESE) $(SOURCESBINARY) $(EXECUTABLEE)

$(EXECUTABLEE): $(OBJECTSE)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSE) -o $@

binary_file_mykernel: $(SOURCES) $(SOURCESF) $(SOURCESBINARY) $(EXECUTABLEF)

$(EXECUTABLEF): $(OBJECTSF)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSF) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(INDEX) $< -o $@

clean:
	rm -rf *.out examples/*.o src/*.o exec/*

tar:
	tar -zcvf BBFMM2D.tar.gz ./exec ./src ./header ./examples ./Makefile ./README ./LICENSE