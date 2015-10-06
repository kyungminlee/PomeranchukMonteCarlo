#define _USE_MATH_DEFINES

#include <cmath>

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>


#include <tclap/CmdLine.h>

#include "pome/basics.h"
#include "pome/pomeranchuk.h"
//#include <boost/math/special_functions/erf.hpp>

using namespace pome;


int main(int argc, char** argv)
{
  try {
    //TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    //TCLAP::ValueArg<std::string> nameArg("n", "name", "Name to print",
    //  true, "homer", "string");
    //TCLAP::ValueArg<Real> l1Arg("l1", "l1", "Length1", true, 1.0, )
    //cmd.add(nameArg);
    //cmd.parse(argc, argv);
    //std::string name = nameArg.getValue();
    
    std::mt19937 gen;

    Real tolerance = 1E-12;
    Integer n_electron = 37;
    Integer n_flux = 2 * n_electron;

    Real l1 = ::sqrt(2 * M_PI * n_flux);
    Real l2 = l1;
    Real angle = M_PI_2;
    Geometry geometry(l1, l2, angle);
    //Real d0 = l1 / n_flux;

    // distance from corners
    auto distance_squared = [l1, l2](Complex r) -> Real {
      Real x = r.real();
      Real y = r.imag();
      Real r1_sq = std::norm(Complex(x     , y     ));
      Real r2_sq = std::norm(Complex(l1 - x, y     ));
      Real r3_sq = std::norm(Complex(x     , l2 - y));
      Real r4_sq = std::norm(Complex(l1 - x, l2 - y));
      return std::min({ r1_sq, r2_sq, r3_sq, r4_sq });
    };
    
    //std::ofstream outfile("out.txt");
 
    // == SELECT OUT DVECTORS

    std::vector<Complex> displacements(geometry.displacements().cbegin(), 
                                       geometry.displacements().cend());
    std::sort(displacements.begin(), displacements.end(), [&distance_squared](Complex c1, Complex c2) { return distance_squared(c1) < distance_squared(c2); });
    if (n_electron > displacements.size() ) {
      std::cout << "n_electron = " << n_electron;
      std::cout << " is too large." << std::endl;
      exit(1);
    } else if (n_electron == displacements.size()) {
      std::cout << "n_electron = " << n_electron;
      std::cout << " is the same as the number of displacements.";
      std::cout << " Doing nothing." << std::endl;
    }
    else {
      std::cout << "n_electron = " << n_electron;
      std::cout << " is smaller than the number of displacements.";
      std::cout << " Truncating." << std::endl;
      std::cout << "r_{n}^2 - r_{n-1}^2 = ";
      std::cout << (distance_squared(displacements[n_electron])
                    - distance_squared(displacements[n_electron - 1]));
      if (n_electron > 1) {
        std::cout << "(cf. r_{n-1}^2 - r_{n-2}^2 = ";
        std::cout << (distance_squared(displacements[n_electron - 1])
                      - distance_squared(displacements[n_electron - 2]));
        std::cout << ")";
      }
      std::cout << std::endl;
      displacements.resize(n_electron);
    }
    //outfile << std::setprecision(16);

    // == Generate Wavefunction
    CflWavefunction wavefunction(geometry, n_electron, displacements.data());
    CflCache cache(n_electron);
    ParticleCoordinates coord(n_electron);

    // initialize monte carlo engine
    // initialize 
    ParticleMoveOpt opt{ 0.1 };
    Measurement measurement(geometry, wavefunction);


    auto const & ks = wavefunction.momentums();
    for (Integer ie = 0; ie < n_electron; ++ie) {
      Complex k1 = ks(ie);
      for (Integer je = 0; je < n_electron; ++je) {
        Complex k2 = ks(je);
        Complex q = k1 - k2;
        if (q.real() >= 0 && q.imag() >= 0) { // TODO ???
          measurement.add_momentum(ie, je, q);
        }
      }
    }


    std::uniform_real_distribution<Real> dist_01(0.0, 1.0);

    for (Integer i = 0; i < n_electron; ++i) {
      Real x = dist_01(gen), y = dist_01(gen);

      Complex z(l1 * x + l2 * y * ::cos(angle),
                l2 * y * sin(angle));
      coord(i) = z;
    }

    wavefunction.generate_theta(coord, cache);
    wavefunction.generate_slater(coord, cache);


    // TODO
    // print wavefunction
    // print coord
    // print cache
    // print measurement
    // print opt

    MonteCarloEngine<PomeranchukProblem> engine(wavefunction,
                                                coord,
                                                cache,
                                                measurement,
                                                opt,
                                                gen,
                                                true,
                                                tolerance
                                                );

    engine.run(100, 1);

    auto stat = engine.get_stat();
    std::cout << "ACCEPT : " << std::get<0>(stat) << std::endl;
    std::cout << "REJECT : " << std::get<1>(stat) << std::endl;
    

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


  return 0;
}
