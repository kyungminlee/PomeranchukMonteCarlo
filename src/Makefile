include make.inc

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

main: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
