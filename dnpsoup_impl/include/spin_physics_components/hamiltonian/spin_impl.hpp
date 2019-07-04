#ifndef DNPSOUP_SPIN_IMPL_HPP
#define DNPSOUP_SPIN_IMPL_HPP

namespace dnpsoup {
  template<SpinType T>
  matrix::Matrix<cxdbl> spin(size_t n)
  {
    if constexpr(T == SpinType::Identity){
      auto mat = matrix::identity<cxdbl>(n);
      return mat;
    }

    auto mat = matrix::zeros<cxdbl>(n, n);
    if(n == 0) return mat;

    const double j = (static_cast<double>(n) - 1.0) * 0.5;
    if constexpr (T == SpinType::X){    // spin x
      /// \f$1/2(I_+)\f$ 
      double m = j - 1.0;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i, i+1) = 0.5 * std::sqrt(j * (j + 1.0) - m * (m + 1.0)); 
        m -= 1.0;
      }

      /// \f$1/2(I_-)\f$
      m = j;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i+1, i) = 0.5 * std::sqrt(j * (j + 1.0) - m * (m - 1.0));
        m -= 1.0;
      }
    } else if constexpr (T == SpinType::Y){   // spin y
      /// \f$-i/2(I_+)\f$
      double m = j - 1.0;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i, i+1) = cxdbl(0,0.5) * std::sqrt(j*(j+1) - m*(m+1));
        m -= 1.0;
      }

      /// \f$i/2(I_-)\f$
      m = j;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i+1, i) = cxdbl(0,0.5) * std::sqrt(j*(j+1) - m*(m-1));
        m -= 1.0;
      }
    } else if constexpr (T == SpinType::Z){   // spin z
      double m = j;
      for(size_t i = 0; i < n; ++i){
        mat(i, i) = m;
        m -= 1.0;
      }
    } else if constexpr (T == SpinType::Plus){    // spin +
      double m = j - 1.0;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i, i+1) = std::sqrt(j * (j + 1.0) - m * (m + 1.0));
        m -= 1.0;
      }
    } else if constexpr (T == SpinType::Minus){   // spin -
      double m = j;
      for(size_t i = 0; i < n - 1; ++i){
        mat(i+1, i) = std::sqrt(j * (j + 1.0) - m * (m - 1.0));
        m -= 1.0;
      }
    } else {
      throw NotImplementedError("Operator symbol not recognized.");
    }

    return mat;
  }
} // namespace dnpsoup

#endif
