#include "Tmy_double.h"

#ifndef Included_T_grad_container_H

#define Included_T_grad_container_H

class T_grad_container
{
private:
	vector<Tmy_double> _grad;
public:
	T_grad_container();
	~T_grad_container();

	void reserve(size_t n);
	void assign(size_t n, Tmy_double value);

	Tmy_double& operator[](size_t idx);
	Tmy_double obj(size_t idx, Tmy_double rho1,Tmy_double rho2);
	Tmy_double dec(size_t idx, Tmy_double rho1,Tmy_double rho2);

};

#endif