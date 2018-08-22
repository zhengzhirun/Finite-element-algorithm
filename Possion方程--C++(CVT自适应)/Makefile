all: CVT_Adaptive_finite

CVT_Adaptive_finite: CVT_Adaptive_finite.o Possion.o Matrix.o density.o
	g++ -std=c++0x -o CVT_Adaptive_finite CVT_Adaptive_finite.o Possion.o Matrix.o density.o

CVT_Adaptive_finite.o: CVT_Adaptive_finite.cpp Possion.h density.h
	g++ -std=c++0x -c CVT_Adaptive_finite.cpp

Possion.o: Possion.cpp Possion.h
	g++ -std=c++0x -c Possion.cpp

Matrix.o: Matrix.cpp Matrix.h
	g++ -std=c++0x -c Matrix.cpp

density.o: density.cpp density.h
	g++ -std=c++0x -c density.cpp
