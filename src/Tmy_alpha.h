#include <iostream>
#include "T_alpha_container.h"
#include "Tmy_double.h"
#include "global.h"

using namespace std;

#ifndef Included_Tmy_alpha_H

#define Included_Tmy_alpha_H

struct Treturn_is_pass
{
  bool is_pass;
  Tmy_double alpha_i;
  Tmy_double alpha_j;
  Tmy_double new_alpha_i;
  Tmy_double new_alpha_j;
};

class Tmy_alpha
{
private:
  Tconfig *_config;

  vector<Tmy_double> calculateBoundaries(int i, int j, T_alpha_container alpha);
  vector<Tmy_double> limit_alpha(Tmy_double alpha_a, Tmy_double alpha_b, Tmy_double Low, Tmy_double High, int flag);
  vector<Tmy_double> calculateNewAlpha(int i, int j, Tmy_double delta, Tmy_double Low, Tmy_double High, T_alpha_container alpha);

public:
  Tmy_alpha(Tconfig *v_config);
  ~Tmy_alpha();
  void init(int jml_data, T_alpha_container& alpha, T_alpha_container& alpha_v1, T_alpha_container& alpha_v2);
  Treturn_is_pass is_pass(int i, int j, Tmy_double delta, T_alpha_container alpha);

};


#endif