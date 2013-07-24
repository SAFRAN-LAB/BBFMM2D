CC=g++
CFLAGS=-c -Wall -DNDEBUG -O3 -ffast-math -ffinite-math-only -I header/ -I ./../
LDFLAGS=

SOURCES=exec/H2_2D_tree.cpp exec/H2_2D_node.cpp #src/kernel.cpp

SOURCESA=src/H2_2D_MVP_input_from_file.cpp


OBJECTSA=$(SOURCES:.cpp=.o) $(SOURCESA:.cpp=.o) 

EXECUTABLEA=./release/H2_2D_MVP_input_from_file


all: $(SOURCES) $(SOURCESA) $(EXECUTABLEA)
	
#$(EXECUTABLEA): $(OBJECTSA)
	#$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSA) -o $@
$(EXECUTABLEB): $(OBJECTSB)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSB) -o $@
$(EXECUTABLEC): $(OBJECTSC)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSC) -o $@
.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(INDEX) $< -o $@

clean:
	rm -rf *.out src/*.o exec/*.o release/*

tar:
	tar -zcvf H2_2D_MVP.tar.gz H2_2D_MVP_input.cpp H2_2D_tree.cpp H2_2D_node.cpp kernel.cpp H2_2D_tree.hpp H2_2D_node.hpp kernel.hpp makefile.mk

taruser:
	tar -zcvf H2_2D_MVP.tar.gz H2_2D_MVP_input.cpp kernel.cpp H2_2D_tree.hpp H2_2D_node.hpp kernel.hpp makefile.mk H2_2D_MVP_input.o H2_2D_tree.o H2_2D_node.o kernel.o
