Concepts

~~~
class Wavefunction {
  virtual Complex operator()(ComplexVector z) const = 0;
};

template <typename Wavefunction>
class WavefunctionCache {
  Complex WavefunctionCache(ComplexVector z); // initializer
  Complex ratio(int i, Complex zi) const;
  Real    ratio_mag(int i, Complex zi) const;
  Complex ratio_phase(int i, Complex zi) const;
  
  Complex update(int i, Complex zi);
  void    recompute();
  void    
  
  
};
~~~
