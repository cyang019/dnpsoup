#ifndef DNPSOUP_SHIFT_H
#define DNPSOUP_SHIFT_H

#include "Matrix.h"
#include "constants.h"
#include "common.h"
#include "errors.h"
#include "spin_physics_components/spin.h"


namespace dnpsoup {
  // abstract base class
  class Shift {
  public:
    Shift(size_t n);
    Shift(size_t n, size_t nbefore, size_t nafter);

    // active rotation
    virtual matrix::Matrix<cxdbl> genMatrix(double gyro, double bz, 
        double sxx, double syy, double szz,
        double alpha, double beta, double gamma) const;

    /// @param sxx: in Hz
    /// @param syy: in Hz
    /// @param szz: in Hz
    /// @param alpha: angle in rads
    /// @param beta: angle in rads
    /// @param gamma: angle in rads
    double calcCoeffX(double sxx, double syy, double szz,
        double alpha, double beta, double gamma) const; 

    /// @param sxx: in Hz
    /// @param syy: in Hz
    /// @param szz: in Hz
    /// @param alpha: angle in rads
    /// @param beta: angle in rads
    /// @param gamma: angle in rads
    double calcCoeffY(double sxx, double syy, double szz,
        double alpha, double beta, double gamma) const; 

    /// @param sxx: in Hz
    /// @param syy: in Hz
    /// @param szz: in Hz
    /// @param beta: angle in rads
    /// @param gamma: angle in rads
    double calcCoeffZ(double sxx, double syy, double szz,
        double beta, double gamma) const; 

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_ix;
    matrix::Matrix<cxdbl> m_iy;
    matrix::Matrix<cxdbl> m_iz;

    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;
  };  // class Shift
}   // namespace dnpsoup

#endif
