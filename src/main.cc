#include <cmath>

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include <tclap/CmdLine.h>
#include <yaml-cpp/yaml.h>

#include "pome/basics.h"
#include "pome/pomeranchuk.h"
#include "montecarlo.h"

using namespace pome;



struct PomeranchukBuilder
{
  std::vector<Complex> displacements;
  ParticleMoveOpt opt;
  ParticleCoordinates coordinates;

  std::unique_ptr<CflWavefunction> p_wavefunction;
  std::unique_ptr<CflCache> p_cache;
  std::unique_ptr<Measurement> p_measurement;

public:
  PomeranchukBuilder()
  {
  }

  void read_configuration(const char* filename)
  {
    YAML::Node config = YAML::LoadFile(filename);
    Integer n_electron = config["n_electron"].as<Integer>();
    auto node_displacement = config["displacements"];
    for (auto disp : node_displacements) {
      double dx = disp[0].as<Real>();
      double dy = disp[1].as<Real>();
    }
  }

  void read_displacements(const char* filename);
  void init();

  void init_wavefunction();
  void init_cache();
  void init_measurement();


};

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

    // == SELECT OUT DVECTORS

#if 1
    std::vector<Complex> displacements(geometry.displacements().cbegin(),
                                       geometry.displacements().cend());

    {
      auto offset = (geometry.tau() + 1.0) * 0.5 * l1;
      for (size_t i = 0; i < displacements.size(); ++i)
      {
        auto d = displacements[i];
        displacements[i] = geometry.mod_z(d + offset) - offset;
      }
    }

    std::sort(displacements.begin(), displacements.end(), [&distance_squared](Complex c1, Complex c2) { return std::norm(c1) < std::norm(c2); });

    if (n_electron > displacements.size()) {
      std::cout << "n_electron = " << n_electron;
      std::cout << " is too large." << std::endl;
      exit(1);
    }
    else if (n_electron == displacements.size()) {
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
#endif

    // == Generate Wavefunction
    CflWavefunction wavefunction(geometry, n_electron, displacements.data());
    CflCache cache(n_electron);
    ParticleCoordinates coord(n_electron);

    // initialize monte carlo engine
    // initialize 
    ParticleMoveOpt opt{ 3.0 * ::sqrt(2.0) * geometry.kappa() };
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

    Integer nel_sqrt = static_cast<Integer>(std::round(::sqrt(n_electron*1.0)));
    //Real nel_sqrt2 = n_electron / nel_sqrt;

    for (Integer i = 0; i < n_electron; ++i) {
      Real x = dist_01(gen), y = dist_01(gen);
      Complex z(l1 * x + l2 * y * ::cos(angle),
                l2 * y * sin(angle));
      coord(i) = z * geometry.kappa();
      //coord(i) = Complex(1.0, 3.0) * ((M_PI * i) / n_electron);
      //coord(i) = displacements[i];
      coord(i) = geometry.mod_kappa_z(coord(i));
    }

#if 0
    std::cout << coord.coordinates << std::endl;
    std::cout << std::endl;
    {
      std::ofstream coordfile("coord.txt");
      coordfile << std::setprecision(16);
      coordfile << std::scientific;
      for (int i = 0; i < n_electron; ++i) {
        coordfile << std::setw(12) << i;
        auto t = coord(i);
        coordfile << std::setw(30) << t.real();
        coordfile << std::setw(30) << t.imag();
        coordfile << std::endl;
      }
    }
#endif
    wavefunction.generate_theta(coord, cache);
#if 0
    {
      std::ofstream thetafile("theta.txt");
      thetafile << std::setprecision(16);
      thetafile << std::scientific;
      for (int i = 0; i < n_electron; ++i) {
        for (int k = 0; k < n_electron; ++k) {
          for (int j = 0; j < n_electron; ++j) {
            thetafile << std::setw(12) << i;
            thetafile << std::setw(12) << k;
            thetafile << std::setw(12) << j;
            auto t = cache.theta(i, k, j);
            thetafile << std::setw(30) << t.real();
            thetafile << std::setw(30) << t.imag();
            thetafile << std::endl;
          }
        }
      }
    }
#endif
    wavefunction.generate_slater(coord, cache);
#if 0
    {
      std::ofstream slaterfile("slater.txt");
      slaterfile << std::setprecision(16);
      slaterfile << std::scientific;

      for (int i = 0; i < n_electron; ++i) {
        for (int j = 0; j < n_electron; ++j) {
          slaterfile << std::setw(12) << i;
          slaterfile << std::setw(12) << j;
          slaterfile << std::setw(30) << cache.slater(i, j).real();
          slaterfile << std::setw(30) << cache.slater(i, j).imag();
          slaterfile << std::endl;
        }
      }
      slaterfile.close();
    }
#endif
    //exit(0);
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

    std::ofstream coordlog("coordlog.txt");
    engine.run(1000, 0);

    for (int i_run = 0; i_run < 20000; ++i_run) {
      engine.step(false);
#if 1
      coordlog << std::setprecision(16);
      for (auto k : coord.coordinates) {
        coordlog << k.real() << "\t" << k.imag() << "\t";
      }
      coordlog << std::endl;
#endif
      if(true){
        auto stat = engine.get_stat();
        std::cout << "ACCEPT : " << std::get<0>(stat) << "\t";
        std::cout << "REJECT : " << std::get<1>(stat) << std::endl;
      }

    }
    coordlog.close();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


  return 0;
}
