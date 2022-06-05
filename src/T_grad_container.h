#include <functional>
#include "Tmy_double.h"
#include "T_alpha_container.h"
#include "Tmy_alpha.h"
#include "Tmy_kernel.h"

#ifndef Included_T_grad_container_H
#define Included_T_grad_container_H

struct callback_param
{
	int idx = -1;
	Tmy_double grad = 0.0;
	Tmy_double obj = 0.0;
	Tmy_double dec = 0.0;
};

using namespace std::placeholders;
using callback_type = std::function<bool(callback_param, callback_param, vector<T_alpha_container>)>;
using callback_type1 = std::function<bool(callback_param, callback_param, vector<T_alpha_container>, Tmy_kernel *, Tmy_alpha *)>;

class T_grad_container
{
private:
	vector<Tmy_double> _grad;
	vector<int> _idx;
	vector<bool> _is_kkt;
	bool _cek_kkt;
	bool _filter_delta;

	bool delta_filter(int idx_b, int idx_a, vector<T_alpha_container> alpha, Tmy_double delta);

public:
	T_grad_container();
	~T_grad_container();

	void reserve(size_t n);
	void assign(size_t n, Tmy_double value);

	Tmy_double &operator[](size_t idx);
	Tmy_double obj(size_t idx, Tmy_double rho1, Tmy_double rho2);
	Tmy_double dec(size_t idx, Tmy_double rho1, Tmy_double rho2);

	void mv_idx(int idx, int flag);

	int max(Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, callback_type f);
	int max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, callback_type f);
	int max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f);

	int cari(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f);
	void set_kkt(int idx, bool val);
};

#endif