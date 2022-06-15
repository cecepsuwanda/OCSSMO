#include "T_grad_container.h"

T_grad_container::T_grad_container()
{
}

T_grad_container::~T_grad_container()
{
	_grad.clear();
	_idx.clear();
	_is_kkt.clear();
}

void T_grad_container::reserve(size_t n)
{
	_grad.reserve(n);
	_idx.reserve(n);
	_is_kkt.reserve(n);
}

void T_grad_container::assign(size_t n, Tmy_double value)
{
	_grad.assign(n, value);
	_is_kkt.assign(n, false);
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

Tcek_opt_return T_grad_container::cek_opt(T_alpha_container alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha)
{
	Tcek_opt_return hasil;
	Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL, delta_gmax = -HUGE_VAL;
	int idx_a = -1, idx_b = -1;
	// cout << " ukuran _grad " << _grad.size() << endl;
	for (size_t i = 0; i < _grad.size(); i++)
	{
		// cout << alpha[i] << "," << alpha.ub() << endl;
		if (alpha[i] >= alpha.ub())
		{
			Tmy_double Gb = _grad[i];
			// cout << "Gb " << Gb << endl;
			Tmy_double tmp = -1.0 * Gb;
			if (tmp >= gmax)
			{
				gmax = tmp;
				idx_b = i;
			}
		}
	}

	if (idx_b != -1)
	{
		Tmy_double Gb = _grad[idx_b];
		for (size_t i = 0; i < _grad.size(); i++)
		{
			if (idx_b != i)
			{
				if (alpha[i] <= alpha.lb())
				{
					Tmy_double Ga = _grad[i];
					// cout << "Ga " << Ga << endl;
					if (Ga >= gmax2)
					{
						gmax2 = Ga;
					}

					vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, i);
					Tmy_double delta = (Ga - Gb) * hsl_eta[0];
					bool is_pass = delta_filter(idx_b, i, alpha, delta);
					if (is_pass)
					{
						if (abs(delta) >= delta_gmax)
						{
							idx_a = i;
							delta_gmax = abs(delta);
						}
					}
				}

				if (idx_a == -1)
				{
					if (alpha[i] >= alpha.ub())
					{
						Tmy_double Ga = _grad[i];
						vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, i);
						Tmy_double delta = (Ga - Gb) * hsl_eta[0];
						bool is_pass = delta_filter(idx_b, i, alpha, delta);
						if (is_pass)
						{
							if (abs(delta) >= delta_gmax)
							{
								idx_a = i;
								delta_gmax = abs(delta);
							}
						}
					}
				}

				if (idx_a == -1)
				{
					Tmy_double Ga = _grad[i];
					vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, i);
					Tmy_double delta = (Ga - Gb) * hsl_eta[0];
					bool is_pass = delta_filter(idx_b, i, alpha, delta);
					if (is_pass)
					{
						if (abs(delta) >= delta_gmax)
						{
							idx_a = i;
							delta_gmax = abs(delta);
						}
					}
				}
			}
		}
	}

	hasil.idx_a = idx_a;
	hasil.idx_b = idx_b;
	hasil.diff = gmax + gmax2;
	hasil.is_opt = hasil.diff < 1e-3;

	return hasil;
}

bool T_grad_container::delta_filter(int idx_b, int idx_a, T_alpha_container alpha, Tmy_double delta)
{
	bool is_pass = true;

	if (alpha.is_nol(idx_b))
	{
		is_pass = (delta > 0.0) and (delta <= alpha[idx_a]);
	}
	else
	{
		if (alpha[idx_b] > 0.0)
		{
			if (delta > 0.0)
			{
				is_pass = (delta <= alpha[idx_a]);
			}
			else
			{
				if (delta < 0.0)
				{
					is_pass = (abs(delta) <= alpha[idx_b]);
				}
			}
		}
	}

	return is_pass;
}

void T_grad_container::set_kkt(int idx, bool val)
{
	_is_kkt[idx] = val;
}

bool T_grad_container::get_kkt(int idx)
{
	return _is_kkt[idx];
}

vector<int> T_grad_container::get_rand_idx()
{
	return _idx;
}