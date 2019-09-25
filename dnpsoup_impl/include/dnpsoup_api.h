#ifndef DNPSOUP_DNPSOUP_API
#define DNPSOUP_DNPSOUP_API


extern "C" {
  void eigenValues(const char *json_str, 
      double *values);    ///< values were stored row major

  double calculateIntensity(const char *json_str);

  double powderIntensity(const char *json_str,
      int ncores);

  void fieldProfile(const char *json_str, double *values, 
      int ncores);
} // extern C


#endif
