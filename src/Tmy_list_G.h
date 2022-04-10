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
   int _active_size;
   Tmy_kernel *_kernel;
   Tmy_list_alpha *_alpha;

   vector<Tmy_double> _arr_G;
   vector<Tmy_double> _arr_G_bar;
   vector<int> _active_set;
public:
   Tmy_list_G(int jml_data, Tmy_kernel *kernel, Tmy_list_alpha *alpha);
   ~Tmy_list_G();

   void clear_container();
   void init();
   void update_G(int idx_b, int idx_a, Tmy_double new_alpha_b, Tmy_double new_alpha_a);
   Tmy_double get_G(int idx);
   void reconstruct_gradient();

   void set_active_size(int new_value);
   void reset_active_size();
   void reverse_swap();
   void swap_index(int i, int j);
};

#endif