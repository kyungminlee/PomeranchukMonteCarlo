theta: theta.cc theta1.f
	g++ -std=c++11 -c theta.cc
	gfortran -c theta1.f
	g++ -std=c++11 theta.o theta1.o 
