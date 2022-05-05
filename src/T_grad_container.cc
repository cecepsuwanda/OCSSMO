#include "T_grad_container.h"


T_grad_container::T_grad_container()
{

}

T_grad_container::~T_grad_container()
{
	_grad.clear();
}

void T_grad_container::reserve(size_t n)
{
	_grad.reserve(n);
}

void T_grad_container::assign(size_t n, Tmy_double value)
{
	_grad.assign(n, value);
}

Tmy_double& T_grad_container::operator[](size_t idx)
{
  return _grad.at(idx);
}

Tmy_double T_grad_container::obj(size_t idx, Tmy_double rho1,Tmy_double rho2)
{
	return min((_grad.at(idx)-rho1),(rho2-_grad.at(idx)));
}

Tmy_double T_grad_container::dec(size_t idx, Tmy_double rho1,Tmy_double rho2)
{
	return ((_grad.at(idx)-rho1)*(rho2-_grad.at(idx)));
}
