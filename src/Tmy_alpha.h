#include <iostream>
#include "Tmy_list_alpha.h"
#include "Tmy_double.h"
#include "global.h"

using namespace std;

#ifndef Included_Tmy_alpha_H

#define Included_Tmy_alpha_H

struct Treturn_is_pass_h
{
  bool is_pass;
  Treturn_is_pass alpha;
  Treturn_is_pass alpha_v1;
  Treturn_is_pass alpha_v2;  
};

class Tmy_alpha
{
private:
	Tconfig *_config;
    Tmy_list_alpha *_my_list_alpha;
    Tmy_list_alpha *_my_list_alpha_v1;
    Tmy_list_alpha *_my_list_alpha_v2;
public:
	Tmy_alpha(Tconfig *v_config);
	~Tmy_alpha();

	void clear_container();

    void init(int jml_data);
    void update_alpha(int idx1,Tmy_double value1,int idx2,Tmy_double value2);

    Tmy_list_alpha* get_alpha();
    Tmy_list_alpha* get_alpha_v1();
    Tmy_list_alpha* get_alpha_v2();

    Treturn_is_pass_h is_pass(int i, int j, Tmy_double delta, Tmy_double delta_v1, Tmy_double delta_v2, int flag);
    
	
};


#endif