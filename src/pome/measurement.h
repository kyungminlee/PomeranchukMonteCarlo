#pragma once

namespace pome
{
class Measurement
{

private:
  const Geometry& geometry_;
  const CflWavefunction & wavefunction_;

  std::vector<Integer> rows_, cols_;         // 
  std::vector<Complex> momentums_;           // q[i] = k[ rows[i] ] - k[ cols[i] ]
  std::vector<Complex> structure_factors_;   // S(q[i]) 

public:

  Measurement(const Geometry & geometry, const CflWavefunction & wavefunction)
    : geometry_(geometry)
    , wavefunction_(wavefunction)
  {
    //TODO
  }

  void measure(const ParticleCoordinates& coord, const CflCache& cache)
  {
    measure_structure_factor(coord.coordinates);
  }


  void reset()
  {
    rows_.clear();
    cols_.clear();
    momentums_.clear();
    structure_factors_.clear();
  }


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
    //static const Complex C_ZERO(0.0, 0.0);
    kore::array::Array<Complex, 1> zs = coordinates.clone();
    zs /= geometry_.kappa();

    assert(rows_.size() == momentums_.size());
    assert(cols_.size() == momentums_.size());

    auto ne = wavefunction_.n_electron();
    for (Integer i = 0, n = momentums_.size(); i < n; ++i) {
      //auto r = rows_[i]; auto c = cols_[i]; 
      auto q = momentums_[i];
      structure_factors_[i] = 0.0;

      for (Integer i_elec = 0; i_elec < ne; ++i_elec) {
        for (Integer j_elec = 0; j_elec < ne; ++j_elec) {
          Complex zr = (zs(i_elec) - zs(j_elec)); // already divided by kappa
          Real prod1 = q.real() * zr.real(), prod2 = q.imag() * zr.imag();
          structure_factors_[i] += Complex(2.0 * ::cos(prod1) * ::cos(prod2),
                                           2.0 * ::sin(prod1) * ::sin(prod2));
        } // for j_elec
      } // for i_elec
    } // i
  }



  auto rows() const -> decltype(rows_) const & { return rows_; }
  auto cols() const -> decltype(cols_) const & { return cols_; }
  auto momentums() const -> decltype(momentums_) const & { return momentums_; }
  auto structure_factors() const -> decltype(structure_factors_) const & { return structure_factors_; }


};


}
