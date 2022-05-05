#include <iostream>
#include <algorithm>
#include <iterator>
#include <map>
#include "Tmy_double.h"
#include "global.h"

using namespace std;


#ifndef Included_Tmy_list_alpha_H

#define Included_Tmy_list_alpha_H

struct Treturn_alpha_stat
{
   Tmy_double jml_alpha;
   int n_all_sv;
   int n_sv;
   Tmy_double jml_alpha_n_sv;
};



class Tmy_list_alpha
{
private:
   int _jml_data;
   Tmy_double _lb;
   Tmy_double _ub;
   Tmy_double _jml_alpha;
   int _n_all_sv;
   int _n_sv;
   Tmy_double _jml_alpha_n_sv;

   vector<Tmy_double> _alpha;   
   vector<int> _alpha_status;
   
   map<int, Tmy_double> _alpha_sv;

   void update_alpha_status(int idx);
   void update_alpha_sv(int idx);   

public:
   Tmy_list_alpha(int v_jml_data, Tmy_double v_lb, Tmy_double v_ub);
   ~Tmy_list_alpha();

   void clear_container();

   void init(Tmy_double V, Tmy_double eps,int flag);
   void update_alpha(int idx, Tmy_double value);   
   
   bool is_nol(int idx_i,int idx_j);   
   bool is_lower_bound(int idx);
   bool is_upper_bound(int idx);
   bool is_free(int idx);
   
   vector<bool> is_alpha_sv(int idx);
   
   Tmy_double get_alpha(int idx);
   map<int, Tmy_double> get_list_alpha_sv();

   Treturn_alpha_stat get_stat();
   Tmy_double get_ub();
   Tmy_double get_lb();  
};


#endif