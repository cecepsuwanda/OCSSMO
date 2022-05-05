#include <string>
#include <thread>
#include "Tmy_kernel.h"
#include "Tmy_list_alpha.h"
#include "Tmy_double.h"
#include "global.h"

using namespace std;

#ifndef Included_Tmy_list_G_H

#define Included_Tmy_list_G_H



class Tmy_list_G
{
private:
   int _jml_data;   
   Tmy_kernel *_kernel;
   Tmy_list_alpha *_alpha;
   vector<Tmy_double> _arr_G;   
   
public:
   Tmy_list_G(int jml_data, Tmy_kernel *kernel, Tmy_list_alpha *alpha);
   ~Tmy_list_G();

   void clear_container();
   void init();
   void update_G(int idx_b, int idx_a, Tmy_double new_alpha_b, Tmy_double new_alpha_a);   
   Tmy_double get_G(int idx);
   Tmy_double get_obj(int idx,Tmy_double rho1,Tmy_double rho2);
   Tmy_double get_dec(int idx,Tmy_double rho1,Tmy_double rho2);   
};

#endif