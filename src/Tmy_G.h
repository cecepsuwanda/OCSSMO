#include <string>
#include <cstdlib>
#include <ctime>
#include "Tmy_list_G.h"
#include "Tmy_kernel.h"
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

struct Treturn_cari_idx
{
	int idx_b;
	int idx_a;
	Tmy_double gmax;
	Tmy_double gmin;
};

class Tmy_G
{
private:
	int _jml_data;	
	Tmy_list_G *_my_list_G;
	Tmy_list_G *_my_list_G_v1;
	Tmy_list_G *_my_list_G_v2;
	Tmy_kernel *_kernel;
	Tmy_alpha *_alphas;

	

	Tmy_double sum_alpha_diff_Q(Tmy_list_alpha* alpha,vector<Tmy_double> diff_Q);
public:
	Tmy_G(int jml_data, Tmy_kernel *kernel, Tmy_alpha *alphas);
	~Tmy_G();

	void clear_container();

	void init();
	Treturn_update_rho update_rho();
	Tmy_list_G* get_list_G();
	Tmy_list_G* get_list_G_v1();
	Tmy_list_G* get_list_G_v2();

	bool is_kkt(int idx, Treturn_update_rho rho);
	Treturn_cari_idx cari_idx(Treturn_update_rho rho);
	int cari_idx_lain(int idx_b,Treturn_update_rho rho);

	void update_G(int idx_b,int idx_a,Treturn_is_pass_h tmp);
	
};

#endif