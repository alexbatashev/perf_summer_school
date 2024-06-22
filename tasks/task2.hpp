/*
    fftbench
    Standard C++ version
    10 May 2003

    Written by Scott Robert Ladd.
    No rights reserved. This is public domain software, for use by anyone.

    A number-crunching benchmark that can be used as a fitness test for
    evolving optimal compiler options via genetic algorithm.

    A simplified version of the FFT component from my standard libcoyote
    library, this benchmark test FFT performance on complex<double> values.

    Note that the code herein is design for the purpose of testing
    computational performance; error handling and other such "niceties"
    is virtually non-existent.

    Actual benchmark results can be found at:
            http://www.coyotegulch.com

    Please do not use this information or algorithm in any way that might
    upset the balance of the universe or to produce imaginary friends for
    imaginary numbers.
*/

#include <complex>

template <class T> class polynomial {
public:
  // creation constructor, uninitialized (coefficients will be unknown values)
  polynomial(size_t degree);

  // creation constructor, with initialization from array
  polynomial(size_t degree, const T *coefficients);

  // creation constructor, with single-value initialization constructor
  polynomial(size_t degree, const T &value);

  // copy constructor
  polynomial(const polynomial<T> &source);

  // destructor
  virtual ~polynomial();

  // assignment
  polynomial<T> &operator=(const polynomial<T> &source);

  // initialize to specific value
  void initialize(const T &value = T(0));

  // increase polynomial length
  polynomial<T> &stretch(size_t degrees);

  // interrogate for degree
  size_t degree() const { return m_degree; }

  // get specific coefficient
  T get(size_t term) const;

  T &operator[](size_t term);

  // evaluate for a specific value
  T operator()(const T &x) const;

  // unitary operators
  polynomial<T> operator-() const;
  polynomial<T> operator+() const;

  // binary mathematical operators
  polynomial<T> operator+(const polynomial<T> &poly) const;
  polynomial<T> operator-(const polynomial<T> &poly) const;

  // constant required by FFT routines
  static const std::complex<T> PI2I;

  // returns largest power of two that holds n
  static size_t log2(size_t n);

  // reverses a sequence of bits
  static size_t flip_bits(size_t k, size_t bits);

  // stretches the length of a polynomial to a power of two
  size_t stretch_fft();

  // performs a reverse-bit copy of a polynomial<T> into a new polynomial<
  // complex<T> >
  static polynomial<std::complex<T>> bit_reverse(const polynomial<T> &poly);
  static polynomial<std::complex<T>>
  bit_reverse(const polynomial<std::complex<T>> &poly);

  // Fast Fourier Transform of polynomial<T> to polynomial< complex<T> >
  static polynomial<std::complex<T>> fft(const polynomial<T> &poly);

  // inverse FFT of polynomial< complex<T> > to polynomial<T>
  static polynomial<std::complex<T>>
  inverse_fft(const polynomial<std::complex<T>> &poly);

  // multiplication via FFT
  polynomial<T> operator*(const polynomial<T> &poly) const;

  // shorthand mathematical operators
  polynomial<T> &operator+=(const polynomial<T> &poly);
  polynomial<T> &operator-=(const polynomial<T> &poly);
  polynomial<T> &operator*=(const polynomial<T> &poly);

  const T* data() const { return m_coeff; };

protected:
  // coefficients
  T *m_coeff;

  // number of terms
  size_t m_degree;

  // acquire resources
  void acquire();

  // release resources
  void release();

  // deep copy
  void deep_copy(const T *source);
};