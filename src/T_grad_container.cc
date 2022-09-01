#include "T_grad_container.h"

T_grad_container::T_grad_container()
{
}

T_grad_container::~T_grad_container()
{
	_grad.clear();

	_is_kkt.clear();
}

void T_grad_container::reserve(size_t n)
{
	_grad.reserve(n);
	_is_kkt.reserve(n);
}

void T_grad_container::assign(size_t n, Tmy_double value)
{
	_grad.assign(n, value);
	_is_kkt.assign(n, false);
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

Tmy_double T_grad_container::obj(size_t idx, Tmy_double rho)
{
	return (_grad.at(idx) - rho);
}

Tmy_double T_grad_container::dec(size_t idx, Tmy_double rho)
{
	return (_grad.at(idx) - rho);
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
