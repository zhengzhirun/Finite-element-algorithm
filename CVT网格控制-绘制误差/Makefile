all: CVT_Error

CVT_Error: CVT_Error.o Matrix.o Possion.o
	g++ -std=c++0x -o CVT_Error CVT_Error.o Possion.o Matrix.o 

CVT_Error.o: CVT_Error.cpp Possion.h 
	g++ -std=c++0x -c CVT_Error.cpp

Possion.o: Possion.cpp Possion.h
	g++ -std=c++0x -c Possion.cpp

Matrix.o: Matrix.cpp Matrix.h
	g++ -std=c++0x -c Matrix.cpp
