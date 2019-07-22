#ifndef DNPSOUP_EULER_H
#define DNPSOUP_EULER_H

namespace dnpsoup {
  class Quaternion;

  // ZYZ rotations
  // Active rotations
  class Euler {
  public:
      Euler() : m_alpha(0.0), m_beta(0.0), m_gamma(0.0) {}
      Euler(double a, double b, double g) : m_alpha(a), m_beta(b), m_gamma(g) {}
      Euler(const Quaternion &);
      Euler(const Euler &) = default;
      Euler(Euler &&) noexcept = default;
      Euler& operator=(const Euler &) = default;
      Euler& operator=(Euler &&) noexcept = default;
      ~Euler() {}

      double alpha() const { return m_alpha; }
      double beta() const { return m_beta; }
      double gamma() const { return m_gamma; }

      Euler& alpha(double a) { m_alpha = a; return *this; }
      Euler& beta(double b) { m_beta = b; return *this; }
      Euler& gamma(double g) { m_gamma = g; return *this; }
  private:
      double m_alpha;
      double m_beta;
      double m_gamma;
  };  // class Euler

  // using quaternion to calculate
  Euler operator*(const Euler &e1, const Euler &e2);
} // namespace dnpsoup

#endif
