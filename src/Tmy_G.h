#include <string>
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
	int _active_size;
	Tmy_list_G *_my_list_G;
	Tmy_list_G *_my_list_G_v1;
	Tmy_list_G *_my_list_G_v2;
	Tmy_kernel *_kernel;
	Tmy_alpha *_alphas;

	std::map<int, Tmy_double> _g_exclude;

	bool _unshrink;
	bool be_shrunk(int i, Tmy_double gmax1, Tmy_double gmax2);
	void swap_index(int i, int j);
	Tmy_double sum_alpha_diff_Q(Tmy_list_alpha* alpha,vector<Tmy_double> diff_Q);
public:
	Tmy_G(int jml_data, Tmy_kernel *kernel, Tmy_alpha *alphas);
	~Tmy_G();

	void clear_container();

	void init();
	Treturn_update_rho update_rho();
	Treturn_update_rho update_rho(Treturn_update_rho old_rho);
	Tmy_list_G* get_list_G();
	Tmy_list_G* get_list_G_v1();
	Tmy_list_G* get_list_G_v2();

	bool is_kkt(int idx, Treturn_update_rho rho);
	Treturn_cari_idx cari_idx();
	Treturn_cari_idx cari_idx(Treturn_update_rho rho);
	int cari_idx_lain(int idx_b);
	int cari_idx_lain(int idx_b,Treturn_update_rho rho);

	void do_shrinking();
	void do_shrinking(Treturn_update_rho rho);
	void reconstruct_gradient();
	int get_active_size();
	void reset_active_size();
	void reverse_swap();

	void update_G(int idx_b,int idx_a,Treturn_is_pass_h tmp,Treturn_update_rho rho);

	void insert_idx_exclude(int idx);
	void delete_idx_exclude();
	bool is_include(Tmy_double abs_F,Treturn_update_rho rho);

};

#endif