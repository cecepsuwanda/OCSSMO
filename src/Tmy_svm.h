#include <string>
#include <iterator>
#include <map>
#include "Tdataframe.h"
#include "T_alpha_container.h"
#include "T_grad_container.h"
#include "Tmy_alpha.h"
#include "Tmy_kernel.h"
#include "Tmy_G.h"
#include "Tmy_double.h"

using namespace std;

#ifndef Included_Tmy_svm_H

#define Included_Tmy_svm_H

struct Treturn_train
{
   int jml_iterasi;
   int n_kkt;
   Tmy_double jml_alpha;
   Tmy_double jml_alpha_v1;
   Tmy_double jml_alpha_v2;
   int n_all_sv;
   int n_sv;
   Tmy_double rho_v1;
   Tmy_double rho_v2;
   bool is_optimum;
};

class Tmy_svm
{

private:
   Tconfig *_config;

   T_alpha_container _alpha;
   T_alpha_container _alpha_v1;
   T_alpha_container _alpha_v2;

   T_grad_container _grad;
   T_grad_container _grad_v1;
   T_grad_container _grad_v2;

   Tmy_alpha *_my_alpha;
   Tmy_kernel *_my_kernel;
   Tmy_G _my_G;

   Treturn_update_rho _rho;
   map<int, vector<string>> _model;
   vector<Tmy_double> _alpha_sv;

public:
   Tmy_svm(Tconfig *v_config);
   ~Tmy_svm();
   Treturn_train train(Tdataframe &df);
   vector<string> test(Tdataframe &df);
   int examineExample();
   int examineExample(int idx_b);

   int take_step(int idx_b, int idx_a);
};

#endif
