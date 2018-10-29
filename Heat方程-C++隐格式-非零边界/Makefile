c++ = g++

VERSION = -std=c++0x

all: Test

Test: Test.o Point.o Mesh.o Parabolic.o TmpEle.o PDE.o
	$(c++) $(VERSION) -o Test Test.o Point.o Mesh.o Parabolic.o TmpEle.o PDE.o

Test.o: Test.cpp
	$(c++) $(VERSION) -c Test.cpp

Point.o: Point.cpp Point.h
	$(c++) $(VERSION) -c Point.cpp

Mesh.o: Mesh.cpp Mesh.h
	$(c++) $(VERSION) -c Mesh.cpp

Parabolic.o: Parabolic.cpp Parabolic.h
	$(c++) $(VERSION) -c Parabolic.cpp

TmpEle.o: TmpEle.cpp TmpEle.h
	$(c++) $(VERSION) -c TmpEle.cpp

PDE.o: PDE.cpp PDE.h
	$(c++) $(VERSION) -c PDE.cpp
