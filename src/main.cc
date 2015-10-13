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
#include <memory>

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
    os << std::setprecision(std::numeric_limits<Real>::max_digits10);
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
    std::unique_ptr<PomeParam> pp(new PomeParam);
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
  auto kDigits = std::numeric_limits<Real>::max_digits10;
  using RNG = std::mt19937_64;

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

  RNG gen(param->random_seed);

  Geometry geometry(param->ell1, param->ell2, param->angle);
  CflWavefunction wavefunction(geometry, param->n_electron, param->displacements.data());
  CflCache cache(param->n_electron);
  ParticleCoordinates coord(param->n_electron);
  param->move_opt.radius *= geometry.kappa(); // TODO
                                              
  std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app;

  std::ofstream result_stream(output_filename.c_str(), mode);
  result_stream.precision(kDigits);
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

  MonteCarloEngine<PomeranchukProblem, RNG> 
    engine(wavefunction, coord, cache, measurement,
           param->move_opt, gen, param->keep_track);

  auto file_loaded = false;
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

  if (!file_loaded) {
    std::uniform_real_distribution<Real> dist_01(0.0, 1.0);
    for (Integer i = 0; i < param->n_electron; ++i) {
      auto x = dist_01(gen), y = dist_01(gen);
      Complex z(param->ell1 * x + param->ell2 * y * ::cos(param->angle),
                param->ell2 * y * sin(param->angle));
      coord(i) = z * geometry.kappa();
      coord(i) = geometry.mod_kappa_z(coord(i));
    }
    wavefunction.generate_theta(coord, cache);
    wavefunction.generate_slater(coord, cache);

    std::cout << "Thermalizing ( " << param->n_thermalize << " steps ) ..." << std::endl;
    engine.run_notrack(param->n_thermalize);
  } else {
    wavefunction.generate_theta(coord, cache);
    wavefunction.generate_slater(coord, cache);
  }


  std::ofstream state_file(state_filename.c_str(), std::ofstream::trunc | std::ofstream::out);
  state_file.precision(kDigits);

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

