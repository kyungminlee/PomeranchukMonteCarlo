#ifndef POMERANCHUK_VERSION
#define POMERANCHUK_VERSION "dev"
#endif

#include <cmath>

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>

#include <tclap/CmdLine.h>
#include <yaml-cpp/yaml.h>

#include "pome/basics.h"
#include "pome/pomeranchuk.h"
#include "montecarlo.h"

using namespace pome;

struct PomeParam {
  std::vector<Complex> displacements;
  Real aspect_ratio;
  Real angle;
  ParticleMoveOpt move_opt;
  Real tolerance;

  // derived parameters
  Integer n_electron;
  Integer n_flux;
  Real ell1, ell2;

  // Engine Option
  Integer n_thermalize;
  Integer n_measure_every;
  Integer n_run;
  Integer random_seed;
  bool keep_track;

  void show(std::ostream& os, const std::string& prefix)
  {
    auto prec = os.precision();
    auto ba = os.boolalpha;
    os << std::setprecision(16);
    os << std::boolalpha;
    os << prefix << "---" << std::endl;
#define _POME_SHOW(X) {os << prefix << #X << " : " << (X) << std::endl; }
    _POME_SHOW(aspect_ratio);
    _POME_SHOW(angle);
    _POME_SHOW(move_opt);
    _POME_SHOW(tolerance);
    _POME_SHOW(n_electron);
    _POME_SHOW(n_flux);
    _POME_SHOW(ell1);
    _POME_SHOW(ell2);
    _POME_SHOW(n_thermalize);
    _POME_SHOW(n_measure_every);
    _POME_SHOW(random_seed);
    _POME_SHOW(keep_track);
#undef _POME_SHOW
    os << prefix << "displacements :" << std::endl;
    for (auto d : displacements) {
      os << prefix << "  - [" << d.real() << ", " << d.imag() << "]" << std::endl;
    }
    os << prefix << "..." << std::endl;
    os.precision(prec);
  }
};


struct PomeranchukReader
{
  PomeranchukReader() = default;

  std::unique_ptr<PomeParam> read_yaml(const char* filename)
  {
    std::unique_ptr<PomeParam> pp = std::make_unique<PomeParam>();
    try {
      PomeParam& pome_param = *pp;

      YAML::Node config = YAML::LoadFile(filename);
      { // displacements
        auto node_displacement = config["displacements"];
        for (auto disp : node_displacement) {
          auto dx = disp[0].as<Real>();
          auto dy = disp[1].as<Real>();
          pome_param.displacements.emplace_back(dx, dy);
        }
      }
      {  // aspect_ratio
        pome_param.aspect_ratio = config["aspect_ratio"].as<Real>();
      }
      { // angle
        auto angle = config["angle"].as<Real>(); // NOTE: IN DEGREES
        pome_param.angle = angle * M_PI / 180.0;
      }


      { // move_opt
        auto node_move_opt = config["move_opt"];
        auto s = node_move_opt["strategy"].as<std::string>();
        pome_param.move_opt.radius = node_move_opt["radius"].as<Real>();
        if (s == "circle") {
          pome_param.move_opt.strategy = ParticleMoveStrategy::CIRCLE;
        }
        else if (s == "square") {
          pome_param.move_opt.strategy = ParticleMoveStrategy::SQUARE;
        }
        else {
          throw YAML::RepresentationException(YAML::Mark::null_mark(),
                                              "strategy not supported.");
        }
      }

      { // monte carlo
        pome_param.n_thermalize = config["n_thermalize"].as<Integer>();
        pome_param.n_measure_every = config["n_measure_every"].as<Integer>();
        pome_param.n_run = config["n_run"].as<Integer>();
        pome_param.random_seed = config["random_seed"].as<uint64_t>();
        pome_param.keep_track = config["keep_track"].as<bool>();
        //pome_param.keep_track = true;
      }

      { // derived
        pome_param.n_electron = pome_param.displacements.size();
        pome_param.n_flux = 2 * pome_param.n_electron;
        pome_param.ell1 = pome_param.ell2 = ::sqrt(2.0 * M_PI * pome_param.n_flux);
        pome_param.tolerance = config["tolerance"].as<Real>();
        if (pome_param.tolerance < 0) {
          throw YAML::RepresentationException(YAML::Mark::null_mark(),
                                              "tolerance should be no less than zero.");
        }
      }


    }
    catch (YAML::Exception & e)
    {
      std::cerr << "Reading file " << filename << " failed (";
      std::cerr << e.what() << ")" << std::endl;
      std::cerr << typeid(e).name() << std::endl;
      pp = nullptr;
    }
    return pp;
  }
};


int main(int argc, char** argv)
{
  std::unique_ptr<PomeParam> param = nullptr;
  std::string output_filename;
  std::string state_filename;
  try {
    TCLAP::CmdLine cmd("PomeranchukMonteCarlo", ' ', POMERANCHUK_VERSION);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg(
      "input_filename",         // name
      "Parameter Filename",     // desc
      true,                     // req
      "test.yaml",              // value
      "string"                  // typeDesc
      );
    TCLAP::ValueArg<std::string> output_filename_arg(
      "o",                      // flag
      "out",                    // name
      "Result Filename",        // desc
      true,                     // req
      "output.txt",             // value
      "string"                  // typeDesc
      );
    TCLAP::ValueArg<std::string> state_filename_arg(
      "s",                      // flag
      "state",                    // name
      "State Filename",        // desc
      false,                     // req
      "state.dat",             // value
      "string"                  // typeDesc
      );

    cmd.add(input_filename_arg);
    cmd.add(output_filename_arg);
    cmd.add(state_filename_arg);
    cmd.parse(argc, argv);
    auto input_filename = input_filename_arg.getValue();
    output_filename = output_filename_arg.getValue();
    state_filename = state_filename_arg.getValue();
    PomeranchukReader reader;
    param = reader.read_yaml(input_filename.c_str());
  }
  catch (TCLAP::ArgException &e) {
    std::cerr << "Arg ERROR" << std::endl;
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  if (!param) { std::cout << "Failed to read file. Quitting." << std::endl; }
  param->show(std::cout, "");

  using RNG = std::mt19937_64;

  RNG gen;

  Geometry geometry(param->ell1, param->ell2, param->angle);
  CflWavefunction wavefunction(geometry, param->n_electron, param->displacements.data());
  CflCache cache(param->n_electron);
  ParticleCoordinates coord(param->n_electron);
  param->move_opt.radius *= geometry.kappa(); // TODO
                                              
  //std::ofstream log_stream("log.txt");
  //log_stream << std::setprecision(16);

  std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app;

  std::ofstream result_stream(output_filename.c_str(), mode);
  result_stream.precision(16);
  result_stream << "# <Parameters>" << std::endl;
  param->show(result_stream, "#   ");
  result_stream << "#" << std::endl;
  result_stream << "# <Data>" << std::endl;
  TableFormatStream result_table_stream(result_stream);
  { // Construct Columns
    result_table_stream.add_column("n_measure");
    for (Integer ie = 0; ie < param->n_electron; ++ie) {
      {
        std::stringstream ss;  ss << "dx(" << ie << ")";
        result_table_stream.add_column(ss.str());
      }
      {
        std::stringstream ss;  ss << "dy(" << ie << ")";
        result_table_stream.add_column(ss.str());
      }
    }
  }
  result_table_stream.print_header();

  Measurement measurement(geometry, wavefunction, result_table_stream);

  bool file_loaded = false;
  {
    std::ifstream state_file(state_filename.c_str());
    try {
      if (state_file)
      {
        Real kz_r, kz_i;
        for (Integer ie = 0; ie < param->n_electron; ++ie) {
          state_file >> kz_r >> kz_i;
          coord.coordinates(ie) = Complex(kz_r, kz_i);
        }
        state_file >> gen;
        std::cout << "Continuuing..." << std::endl;
        file_loaded = true;
      }
    }
    catch (std::invalid_argument& e) {
      std::cout << "Failed reading state file." << std::endl;
      // pass
    }
  }

  wavefunction.generate_theta(coord, cache);
  wavefunction.generate_slater(coord, cache);
  MonteCarloEngine<PomeranchukProblem, RNG> engine(wavefunction,
                                              coord,
                                              cache,
                                              measurement,
                                              param->move_opt,
                                              gen,
                                              param->keep_track,
                                              param->tolerance
                                              );
  if (!file_loaded) {
    std::uniform_real_distribution<Real> dist_01(0.0, 1.0);
    for (Integer i = 0; i < param->n_electron; ++i) {
      Real x = dist_01(gen), y = dist_01(gen);
      Complex z(param->ell1 * x + param->ell2 * y * ::cos(param->angle),
                param->ell2 * y * sin(param->angle));
      coord(i) = z * geometry.kappa();
      coord(i) = geometry.mod_kappa_z(coord(i));
    }

    std::cout << "Thermalizing ( " << param->n_thermalize << " steps ) ..." << std::endl;
    engine.run_notrack(param->n_thermalize);
  }

  std::ofstream state_file(state_filename.c_str(), std::ofstream::trunc | std::ofstream::out);
  state_file.precision(16);
  std::cout << std::fixed << std::setprecision(1);
  for (Integer i_run = 0; i_run < param->n_run; ++i_run) {
    auto start_time = std::chrono::steady_clock::now();
    engine.run(param->n_measure_every);

    { // save state
      state_file.seekp(std::ofstream::beg);
      for (Integer ie = 0; ie < param->n_electron; ++ie) {
        auto kz = coord.coordinates(ie);
        state_file << kz.real() << "\t" << kz.imag() << std::endl;
      }
      state_file << gen << std::endl;
      state_file.flush();
    }


    measurement.measure(coord, cache);
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() * 1.0 / param->n_measure_every;
    //auto time_left = (end_time - start_time) * (param->n_measure_every * (param->n_run - i_run - 1)*1.0);
    auto stat = engine.get_stat();
    auto nac = std::get<0>(stat);
    auto nrj = std::get<1>(stat);
    std::cout << "Acceptance = " << nac << " / " << (nac + nrj);
    std::cout << " = " << std::setprecision(1) << (nac * 100.0 / (nac + nrj)) << "%";
    std::cout << " ( " << std::setprecision(3) << duration << " ms/step";
    std::cout << " : " << std::setprecision(1) << (1000.0 / duration) << " step/sec";
    //std::cout << " ,  ETA = " << std::setprecision(0) << time_left << " sec";
    std::cout << " )" << std::endl;
  }

  return 0;
}




#if 0
int test(int argc, char** argv)
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

    std::sort(displacements.begin(), displacements.end(), [](Complex c1, Complex c2) { return std::norm(c1) < std::norm(c2); });

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
      std::cout << (std::norm(displacements[n_electron])
                    - std::norm(displacements[n_electron - 1]));
      if (n_electron > 1) {
        std::cout << "(cf. r_{n-1}^2 - r_{n-2}^2 = ";
        std::cout << (std::norm(displacements[n_electron - 1])
                      - std::norm(displacements[n_electron - 2]));
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
#if 0
    std::ofstream foo;
    Measurement measurement(geometry, wavefunction, false, foo);

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
#endif

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


#endif