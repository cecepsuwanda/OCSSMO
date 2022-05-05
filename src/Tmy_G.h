#include <string>
#include <cstdlib>
#include <ctime>
#include "Tmy_kernel.h"
#include "T_alpha_container.h"
#include "T_grad_container.h"
#include "Tmy_alpha.h"
#include "Tmy_double.h"

using namespace std;

#ifndef Included_Tmy_G_H

#define Included_Tmy_G_H

struct Treturn_update_rho
{
	Tmy_double rho_v1;
	Tmy_double rho_v2;
};

class Tmy_G
{
private:
	int _jml_data;		
	Tmy_double sum_alpha_rho_Q(int i,Tmy_kernel *kernel, T_alpha_container alpha);
public:
	Tmy_G();
	~Tmy_G();

	void init(int jml_data, Tmy_kernel *kernel, T_alpha_container alpha,T_grad_container &grad);
	Treturn_update_rho update_rho(Tmy_kernel *kernel, T_alpha_container alpha,T_grad_container grad);
	
	bool is_kkt(int idx, Treturn_update_rho rho,T_alpha_container alpha,T_grad_container grad);
	int cari_idx_a(int idx_b, Treturn_update_rho rho,T_alpha_container alpha,T_grad_container grad);
	
	int cari_idx_lain(int idx_b, Treturn_update_rho rho,Tmy_kernel *kernel,T_alpha_container alpha,T_grad_container grad,Tmy_alpha *my_alpha);
	
	bool cari_idx(int& idx_b, int& idx_a, Treturn_update_rho rho);

	Tmy_double sum_alpha_diff_Q(T_alpha_container alpha, vector<Tmy_double> diff_Q);
	void update_G(int idx_b, int idx_a, Treturn_is_pass tmp, Tmy_kernel *kernel, T_alpha_container &alpha,T_grad_container &grad);
	 
};

#endif