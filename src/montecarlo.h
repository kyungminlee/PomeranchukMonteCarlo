#pragma once

#include <random>

class MonteCarloStat;
class MonteCarloManager;

class MonteCarloStat
{
public:
  explicit MonteCarloStat(bool keep_track)
    : keep_track_(keep_track)
    , n_accept_(0), n_reject_(0)
  { }

  void accept() {
    n_accept_++;
    if (keep_track_) { history_.push_back(true); }
  }

  void reject() {
    n_reject_++;
    if (keep_track_) { history_.push_back(false); }
  }

  void keep_track(bool kt) { keep_track_ = kt; }
  bool keep_track() const { return keep_track_; }

  int n_accept() const { return n_accept_; }
  int n_reject() const { return n_reject_; }
  int n_total() const { return n_accept_ + n_reject_; }
  const std::vector<bool>& history() const { return history_; }

private:
  bool keep_track_;
  int n_accept_, n_reject_;
  std::vector<bool> history_;
};

template <typename _ProblemType>
class MonteCarloEngine
{
public:
  using ProblemType = _ProblemType;
  using SystemType = typename ProblemType::SystemType;
  using StateType = typename ProblemType::StateType;
  using UpdateType = typename ProblemType::UpdateType;
  using UpdateOptType = typename ProblemType::UpdateOptType;
  using CacheType = typename ProblemType::CacheType;
  using MeasurementType = typename ProblemType::MeasurementType;
  using RealType = typename ProblemType::RealType;

  MonteCarloEngine(
    const SystemType& system,
    StateType& state,
    CacheType& cache,  // state and cache changes all the time
    MeasurementType & measurement,
    const UpdateOptType& update_opt,
    std::mt19937 & gen,
    bool keep_track = false,
    RealType tolerance = 1E-15)
    : system_(system)
    , state_(state)
    , cache_(cache)
    , update_opt_(update_opt)
    , measurement_(measurement)
    , stat_(keep_track)
    , tolerance_(tolerance)
    , random_generator_(gen)
    , real_dist_(0.0, 1.0)
  {
    assert(tolerance > 0);
  }

  void step(bool measure)
  {
    UpdateType update = system_.candidate(state_, cache_, update_opt_, random_generator_);
    //std::cout << update.i_electron << "\t" << update.z_new << std::endl; 
    system_.update(update, state_, cache_);
    RealType p = system_.acceptance(update, state_, system_, cache_);
    if (p >= 1.0) {  system_.commit(update, state_, cache_);      stat_.accept(); }
    else {
      RealType p2 = real_dist_(random_generator_);
      if (p2 <= p) { system_.commit(update, state_, cache_);      stat_.accept(); }
      else {         system_.revert(update, state_, cache_);      stat_.reject(); }
    }
    if (measure) { measurement_.measure(state_, cache_); }
  }
  
  void run(size_t n_step, size_t measure_every)
  {
    if (measure_every == 0) {
      for (size_t i_step = 0; i_step < n_step; ++i_step) { step(false); }
    } // if measure_every == 0 
    else {
      for (size_t i_step = 0; i_step < n_step; ++i_step) {
        if ((i_step + 1) % measure_every == 0) { step(true); }
        else { step(false); }
      } // for i_step
    } // else measure_every == 0
  }

  std::tuple<int, int> get_stat() const {
    return std::make_tuple(stat_.n_accept(), stat_.n_reject());
  }

private:
  const SystemType& system_;
  StateType & state_;
  CacheType & cache_;
  const UpdateOptType & update_opt_;
  MeasurementType & measurement_;
  MonteCarloStat stat_;

  RealType tolerance_;
  std::mt19937 & random_generator_;
  std::uniform_real_distribution<double> real_dist_;
};





#if 0
void step()
{
  Complex old_log_determinant = wavefunction_cache_.log_determinant;

  //// ==== NEXT_MOVE
  Integer i_elec = particle_dist_(random_generator_);
  //Complex z_new = real_dist_(random_generator_) + geometry_.tau() * real_dist_(random_generator_);
  Complex z_diff = real_dist_(random_generator_) + geometry_.tau() * real_dist_(random_generator_);
  Complex z_new = wavefunction_cache_.coordinates(i_elec) + z_diff;
  // TODO use different distribution

  z_new = geometry_.mod_kappa_z(z_new);

  //// ==== UPDATE
  wavefunction_.update_theta(i_elec, z_new, wavefunction_cache_);
  wavefunction_.generate_slater(wavefunction_cache_);

  Complex log_determinant = wavefunction_cache_.log_determinant;
  Real old_log_probability = old_log_determinant.real();
  Real log_probability = log_determinant.real();

  if (log_probability > old_log_probability) {
    // accept
    wavefunction_.accept_theta(i_elec, z_new, wavefunction_cache_);
  }
  else {
    Real diff_log_probability = log_probability - old_log_probability;
    Real dice = ::log(real_dist_(random_generator_));
    if (dice > diff_log_probability) {
      wavefunction_.accept_theta(i_elec, z_new, wavefunction_cache_);
    }
    else {
      wavefunction_.revert_theta(i_elec, wavefunction_cache_);
    }
  }
  // compute Sq
}
#endif
