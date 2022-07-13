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
   int n_all_sv;
   int n_sv;

   int n_kkt_v1;
   int n_all_sv_v1;
   int n_sv_v1;

   int n_kkt_v2;
   int n_all_sv_v2;
   int n_sv_v2;

   Tmy_double jml_alpha;
   Tmy_double jml_alpha_v1;
   Tmy_double jml_alpha_v2;

   Tmy_double rho_v1;
   Tmy_double rho_v2;

   bool is_optimum;

   friend ostream &operator<<(ostream &out, const Treturn_train &rt)
   {
      out << " Jumlah iterasi  : " << rt.jml_iterasi << endl;
      out << " Jumlah alpha    : " << rt.jml_alpha << endl;
      out << " Jumlah alpha v1 : " << rt.jml_alpha_v1 << endl;
      out << " Jumlah alpha v2 : " << rt.jml_alpha_v2 << endl;
      out << " rho v1          : " << rt.rho_v1 << endl;
      out << " rho v2          : " << rt.rho_v2 << endl;
      out << " alpha " << setw(10) << " kkt " << setw(10) << " not_nol " << setw(10) << "sv" << endl;
      out << "   1   " << setw(10) << rt.n_kkt_v1 << setw(10) << rt.n_all_sv_v1 << setw(10) << rt.n_sv_v1 << endl;
      out << "   2   " << setw(10) << rt.n_kkt_v2 << setw(10) << rt.n_all_sv_v2 << setw(10) << rt.n_sv_v2 << endl;
      out << "   0   " << setw(10) << rt.n_kkt << setw(10) << rt.n_all_sv << setw(10) << rt.n_sv << endl;

      return out;
   }
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
   Treturn_update_rho _rho_1;

   map<int, vector<string>> _model;
   vector<Tmy_double> _alpha_sv;

public:
   Tmy_svm(Tconfig *v_config);
   ~Tmy_svm();
   Treturn_train train(Tdataframe &df);
   vector<string> test(Tdataframe &df);
   bool examineExample();
   int examineExample(int idx_b);
   bool take_step(int idx_b, int idx_a);

   bool examineExample_v();
   int examineExample_v(int idx_b);
   bool take_step_v(int idx_b, int idx_a, Tmy_double &rho, T_alpha_container &alpha, T_grad_container &grad);
};

#endif
