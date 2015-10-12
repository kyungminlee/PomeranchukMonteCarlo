HEADERS = \
	montecarlo.h \
	pome/basics.h \
	pome/geometry.h \
	pome/geometry_impl.h \
	pome/measurement.h \
	pome/naive_theta.h \
	pome/particle.h \
	pome/pomeranchuk.h \
	pome/theta.h \
	pome/wavefunction.h \

SOURCES = main.cc

LDFLAGS = -L$(HOME)/.brew/opt/openblas/lib -lopenblas
CXXFLAGS = -g -Wall -std=c++11 -I$(HOME)/.brew/opt/openblas/include -I$(HOME)/.brew/include/eigen3 -I../../koretools/include -I$(HOME)/.brew/include
main: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
