#include "T_alpha_container.h"


T_alpha_container::T_alpha_container()
{

}

T_alpha_container::~T_alpha_container()
{
	_alpha.clear();
}

void T_alpha_container::reserve(size_t n)
{
	_alpha.reserve(n);
}

void T_alpha_container::assign(size_t n, Tmy_double value)
{
	_alpha.assign(n, value);
}

void T_alpha_container::boundaries(Tmy_double lb, Tmy_double ub)
{
	_lb = lb;
	_ub = ub;
}

Tmy_double T_alpha_container::sum()
{
	Tmy_double total = 0.0;
	for (auto it : _alpha)
	{
		total = total + ((Tmy_double) it);
	}
	return total;
}

int T_alpha_container::n_all_sv()
{
	int n = 0;
	int i = 0;
	for (auto it : _alpha)
	{
		if (_alpha.at(i) != 0.0)
		{
			n = n + 1;
		}
		i = i + 1;
	}
	return n;
}

int T_alpha_container::n_sv()
{
	int n = 0;
	int i = 0;
	for (auto it : _alpha)
	{
		if (is_sv(i))
		{
			n = n + 1;
		}
		i = i + 1;
	}
	return n;
}

bool T_alpha_container::is_sv(size_t idx)
{
	return  (((_lb < _alpha.at(idx)) and (_alpha.at(idx) < _ub)) and _alpha.at(idx) != 0.0);
}

bool T_alpha_container::is_ub(size_t idx)
{
  return _alpha.at(idx) == _ub;
}

bool T_alpha_container::is_lb(size_t idx)
{
  return _alpha.at(idx) == _lb;
}

Tmy_double T_alpha_container::ub()
{
	return _ub;
}

Tmy_double T_alpha_container::lb()
{
	return _lb;
}

Tmy_double& T_alpha_container::operator[](size_t idx)
{
	return _alpha.at(idx);
}