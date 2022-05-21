#include "T_grad_container.h"

T_grad_container::T_grad_container()
{
}

T_grad_container::~T_grad_container()
{
	_grad.clear();
	_idx.clear();
}

void T_grad_container::reserve(size_t n)
{
	_grad.reserve(n);
	_idx.reserve(n);
}

void T_grad_container::assign(size_t n, Tmy_double value)
{
	_grad.assign(n, value);
	_idx.assign(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		_idx[i] = i;
	}
}

Tmy_double &T_grad_container::operator[](size_t idx)
{
	return _grad.at(idx);
}

Tmy_double T_grad_container::obj(size_t idx, Tmy_double rho1, Tmy_double rho2)
{
	return min((_grad.at(idx) - rho1), (rho2 - _grad.at(idx)));
}

Tmy_double T_grad_container::dec(size_t idx, Tmy_double rho1, Tmy_double rho2)
{
	return ((_grad.at(idx) - rho1) * (rho2 - _grad.at(idx)));
}

void T_grad_container::mv_idx(int idx, int flag)
{
	int i = 0;
	bool ketemu = false;
	while (!ketemu and (i < _idx.size()))
	{
		ketemu = _idx[i] == idx;
		i++;
	}
	if (ketemu)
	{
		_idx.erase(_idx.begin() + (i - 1));
		if (flag == 0)
		{
			_idx.insert(_idx.begin(), idx);
		}
		else
		{
			_idx.push_back(idx);
		}
	}
}

int T_grad_container::max(Tmy_double rho1, Tmy_double rho2, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, callback_type f)
{
	Tmy_double gmax = -HUGE_VAL;
	int idx_max = -1;
	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double F = obj(tmp_idx[i], rho1, rho2);
		Tmy_double abs_F = abs(F);
		bool is_pass = f(-1, tmp_idx[i], F, 0.0, alpha, alpha_v1, alpha_v2);
		if (is_pass)
		{
			if (abs_F >= gmax)
			{
				gmax = abs_F;
				idx_max = tmp_idx[i];
				mv_idx(tmp_idx[i], 0);
			}
		}
	}
	return idx_max;
}

int T_grad_container::max(int idx_b, Tmy_double rho1, Tmy_double rho2, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, Tmy_kernel *kernel, callback_type f)
{
	Tmy_double gmax = -HUGE_VAL;
	int idx_max = -1;

	Tmy_double Fb = obj(idx_b, rho1, rho2);
	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Fa = obj(tmp_idx[i], rho1, rho2);
		Tmy_double diff_F = abs(Fb - Fa);
		bool is_pass = f(idx_b, tmp_idx[i], Fa, Fb, alpha, alpha_v1, alpha_v2);
		if (is_pass)
		{
			if (diff_F >= gmax)
			{
				gmax = diff_F;
				idx_max = tmp_idx[i];
				mv_idx(tmp_idx[i], 0);
			}
		}

	}
	return idx_max;
}

int T_grad_container::max(int idx_b, Tmy_double rho1, Tmy_double rho2, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
	Tmy_double gmax = -HUGE_VAL;
	int idx_max = -1;

	Tmy_double Fb = obj(idx_b, rho1, rho2);
	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Fa = obj(tmp_idx[i], rho1, rho2);
		Tmy_double diff_F = abs(Fb - Fa);
		bool is_pass = f(idx_b, tmp_idx[i], Fa, Fb, alpha, alpha_v1, alpha_v2, kernel, my_alpha);
		if (is_pass)
		{
			if (diff_F >= gmax)
			{
				gmax = diff_F;
				idx_max = tmp_idx[i];
				mv_idx(tmp_idx[i], 0);
			}
		}

	}
	return idx_max;
}

int T_grad_container::cari(int idx_b, Tmy_double rho1, Tmy_double rho2, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
	int idx_a = -1;

	Tmy_double Fb = obj(idx_b, rho1, rho2);
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Fa = obj(_idx[i], rho1, rho2);
		bool is_pass = f(idx_b, _idx[i], Fa, Fb, alpha, alpha_v1, alpha_v2, kernel, my_alpha);
		if (is_pass)
		{
			idx_a = _idx[i];
			break;
		}

	}

	return idx_a;
}
