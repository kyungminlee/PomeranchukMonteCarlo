#pragma once

#include "../tableformatstream.h"

namespace pome
{
class Measurement
{

private:
  const Geometry& geometry_;
  const CflWavefunction & wavefunction_;
  TableFormatStream & output_stream_;

  //std::vector<Integer> rows_, cols_;         // 
  //std::vector<Complex> momentums_;           // q[i] = k[ rows[i] ] - k[ cols[i] ]

  Integer n_measure_;
  //std::vector<Complex> structure_factors_;   // S(q[i]) 
  //std::vector< std::vector<Complex> > structure_factor_history_;

public:

  Measurement(Geometry const & geometry,
              CflWavefunction const & wavefunction,
              TableFormatStream & os)
    : geometry_(geometry)
    , wavefunction_(wavefunction)
    , output_stream_(os)
    , n_measure_(0)
  {
  }

  void measure(const ParticleCoordinates& coord,
               const CflCache& cache)
  {
    output_stream_ << ++n_measure_;
    for (auto kz : coord.coordinates) {
      Complex z = kz * geometry_.inverse_kappa();
      output_stream_ << z.real() << z.imag();
    }
    output_stream_.linebreak();
  }

#if 0
  void reset()
  {
    rows_.clear();
    cols_.clear();
    momentums_.clear();
    structure_factors_.clear();
  }
#endif

#if 0
  void add_momentum(Integer row, Integer col, Complex q)
  {
    rows_.push_back(row);
    cols_.push_back(col);
    momentums_.push_back(q);
    structure_factors_.push_back(0.0);
  }

  void clear_structure_factors()
  {
    std::fill(structure_factors_.begin(), structure_factors_.end(), Complex(0.0));
  }

  void measure_structure_factor(kore::array::Array<Complex const ,1> const & coordinates)
  {
    ++n_measure_;
    kore::array::Array<Complex, 1> zs = coordinates.clone();
    zs /= geometry_.kappa();

    assert(rows_.size() == momentums_.size());
    assert(cols_.size() == momentums_.size());
    assert(structure_factors_.size() == momentums_.size());

    auto ne = wavefunction_.n_electron();
    for (Integer i = 0, n = momentums_.size(); i < n; ++i) {
      //auto r = rows_[i]; auto c = cols_[i]; 
      auto ki = momentums_[rows_[i]];
      auto kj = momentums_[cols_[i]];
      Complex str_fac(0.0, 0.0);
      for (Integer i_elec = 0; i_elec < ne; ++i_elec) {
        for (Integer j_elec = 0; j_elec < ne; ++j_elec) {
          Complex zr = zs(i_elec) - zs(j_elec); // already divided by kappa
          Real prod1 = ki.real() * zr.real(), prod2 = q.imag() * zr.imag();
          auto skp = Complex(2.0 * ::cos(prod1) * ::cos(prod2),
                            2.0 * ::sin(prod1) * ::sin(prod2));
          str_fac += skp;
        } // for j_elec
      } // for i_elec
      structure_factors_[i] += str_fac;
      if (keep_track_) { output_stream_ << str_fac.real() << str_fac.imag(); }
    } // i
    if (keep_track_) { output_stream_.linebreak(); }
  }

#endif

  Integer n_measure() const { return n_measure_; }
  //auto rows() const -> decltype(rows_) const & { return rows_; }
  //auto cols() const -> decltype(cols_) const & { return cols_; }
  //auto momentums() const -> decltype(momentums_) const & { return momentums_; }
  //auto structure_factors() const -> decltype(structure_factors_) const & { return structure_factors_; }
  //auto structure_factor_history() const -> decltype(structure_factor_history_) const & { return structure_factor_history_; }
};


}
