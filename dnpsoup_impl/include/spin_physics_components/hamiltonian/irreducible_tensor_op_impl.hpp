namespace dnpsoup {
  //Ernst, Richard R., Geoffrey Bodenhausen, and Alexander Wokaun. 
  //Principles of nuclear magnetic resonance in one and two dimensions. 1987.
  // Chapter 2

  /// single spin tensors
  template<int l, int m>
  inline
  MatrixCxDbl tensor(size_t n)
  {
    if constexpr(l < 0){
      throw IndexError("Tensor rank cannot be negative.");
    }

    if constexpr(l > 1){
      throw IndexError("Tensor rank for a single spin cannot be larger than 2.");
    }

    if constexpr(m > l || m < -l){
      throw IndexError("m cannot be larger than l for tensor ranks.");
    }

    if constexpr(l == 0){
      auto res = identity<cxdbl>(n) * (0.5 * sqrt(2));
      return res;
    }
    else {  // l == 1
      if constexpr(m == 1) 
      { return -spin<P>(n); }
      else if constexpr(m == 0)
      { return sqrt(2) * spin<Z>(n); }
      else { // m == -1
        { return spin<M>(n); }
      }
    }
  }

  /// two spin tensors
  template<int l, int m>
  inline
  MatrixCxDbl tensor(size_t n1, size_t n2)
  {
    if constexpr(l < 0){
      throw IndexError("Tensor rank cannot be negative.");
    }

    if constexpr(l > 2){
      throw IndexError("Tensor rank for 2-spin system cannot be larger than 2.");
    }

    if constexpr(m > l || m < -l){
      throw IndexError("m cannot be larger than l for tensor ranks.");
    }

    if constexpr(l == 0){
      auto res = (-1.0/sqrt(3.0)) * (kron(spin<P>(n1), spin<M>(n2)) + kron(spin<M>(n1), spin<P>(n2))) 
        + (-2.0/sqrt(3.0)) * kron(spin<Z>(n1), spin<Z>(n2));
      return res;
    } 
    else if constexpr(l == 1){
      if constexpr(m == 0){
        auto res = (0.5 * sqrt(2.0)) * (kron(spin<M>(n1), spin<P>(n2)) - kron(spin<P>(n1), spin<M>(n2)));
        return res;
      } else if constexpr(m == 1){
        auto res = kron(spin<Z>(n1), spin<P>(n2)) - kron(spin<P>(n1), spin<Z>(n2));
        return res;
      } else {  // m == -1
        auto res = kron(spin<Z>(n1), spin<M>(n2)) - kron(spin<M>(n1), spin<Z>(n2));
        return res;
      }
    } else {  // l == 2
      if constexpr(m == 0){
        auto res = sqrt(2.0/3.0) * (2.0 * kron(spin<Z>(n1), spin<Z>(n2)) 
            - kron(spin<X>(n1), spin<X>(n2))
            - kron(spin<Y>(n1), spin<Y>(n2)));
        return res;
      } else if constexpr(m == 1){
        auto res = -kron(spin<P>(n1), spin<Z>(n2)) - kron(spin<Z>(n1), spin<P>(n2));
        return res;
      } else if constexpr(m == -1){
        auto res = kron(spin<M>(n1), spin<Z>(n2)) + kron(spin<Z>(n1), spin<M>(n2));
        return res;
      } else if constexpr(m == 2){
        auto res = kron(spin<P>(n1), spin<P>(n2));
        return res;
      } else {  // m == -2
        return kron(spin<M>(n1), spin<M>(n2));
      }
    }
  }
} // namespace dnpsoup
