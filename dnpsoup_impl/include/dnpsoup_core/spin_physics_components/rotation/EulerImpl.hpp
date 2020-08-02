namespace dnpsoup {
  template<typename T>
  Euler<T>::Euler(const Euler<T> &e)
    : m_alpha(e.m_alpha), m_beta(e.m_beta), m_gamma(e.m_gamma)
  {}

  template<typename T>
  Euler<T>::Euler(Euler<T> &&e) noexcept
    : m_alpha(std::move(e.m_alpha)), 
    m_beta(std::move(e.m_beta)), 
    m_gamma(std::move(e.m_gamma))
  {}

  template<typename T>
  Euler<T>& Euler<T>::operator=(const Euler<T> &rhs)
  {
    m_alpha = rhs.m_alpha;
    m_beta = rhs.m_beta;
    m_gamma = rhs.m_gamma;
    return *this;
  }

  template<typename T>
  Euler<T>& Euler<T>::operator=(Euler<T> &&rhs) noexcept
  {
    m_alpha = std::move(rhs.m_alpha);
    m_beta = std::move(rhs.m_beta);
    m_gamma = std::move(rhs.m_gamma);
    return *this;
  }

  template<typename R>
  inline
  Euler<R> toEuler(const Quaternion &q)
  {
      static_assert(is_rotation_type<R>::value, 
          "Euler can either be ActiveRotation or PassiveRotation");
      
      double mat02 = 2.0 * (q.i() * q.k() + q.r() * q.j());
      double mat12 = 2.0 * (q.j() * q.k() - q.r() * q.i());
      double mat20 = 2.0 * (q.i() * q.k() - q.r() * q.j());
      double mat21 = 2.0 * (q.j() * q.k() + q.r() * q.i());
      double mat22 = q.r() * q.r() - q.i() * q.i() - q.j() * q.j() + q.k() * q.k();
  
      double alpha, beta, gamma;
      if constexpr(std::is_same<R, ActiveRotation>::value){
        alpha = dnpsoup::atan(mat12, mat02);
        beta = dnpsoup::atan(std::sqrt(1.0 - mat22 * mat22), mat22);
        gamma = dnpsoup::atan(mat21, -mat20);
      } else {
        gamma = -dnpsoup::atan(mat12, mat02);
        beta = -dnpsoup::atan(std::sqrt(1.0 - mat22 * mat22), mat22);
        alpha = -dnpsoup::atan(mat21, -mat20);
      }

      /// range 0 ~ pi
      //alpha = alpha > -eps ? alpha : -alpha;
      //beta = beta > -eps ? beta : -beta;
      //gamma = gamma > -eps ? gamma : -gamma;
      return Euler(alpha, beta, gamma);
  }

  template<typename R>
  inline
  Quaternion toQuaternion(const Euler<R> &e)
  {
    static_assert(is_rotation_type<R>::value, 
        "Euler can either be ActiveRotation or PassiveRotation");
    double r, i, j, k;
    if constexpr (std::is_same<R, ActiveRotation>::value){
      r = cos(0.5 * e.beta()) * cos(0.5 * (e.alpha() + e.gamma()));
      i = -sin(0.5 * e.beta()) * sin(0.5 * (e.alpha() - e.gamma()));
      j = sin(0.5 * e.beta()) * cos(0.5 * (e.alpha() - e.gamma()));
      k = cos(0.5 * e.beta()) * sin(0.5 * (e.alpha() + e.gamma()));
    } else {
      r = cos(0.5 * e.beta()) * cos(0.5 * (e.alpha() + e.gamma()));
      i = sin(0.5 * e.beta()) * sin(0.5 * (e.alpha() - e.gamma()));
      j = -sin(0.5 * e.beta()) * cos(0.5 * (e.alpha() - e.gamma()));
      k = -cos(0.5 * e.beta()) * sin(0.5 * (e.alpha() + e.gamma()));
    }
    return Quaternion(r,i,j,k);
  }

  template<typename T>
  inline
  Euler<T> operator*(const Euler<T> &e1, const Euler<T> &e2)
  {
    Quaternion q1 = toQuaternion(e1);
    Quaternion q2 = toQuaternion(e2);
    Quaternion res = q1 * q2;
    return toEuler<T>(res);
  }

  template<typename T>
  inline
  std::ostream& operator<<(std::ostream &os, const Euler<T> &euler)
  {
		std::streamsize ss = std::cout.precision();
		os << std::setprecision(14)
       << euler.alpha() << ", " << euler.beta() << ", " << euler.gamma();
		os << std::setprecision(ss);
    return os;
  }
} // namespace dnpsoup

