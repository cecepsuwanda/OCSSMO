#include "Tmy_alpha.h"


Tmy_alpha::Tmy_alpha(Tconfig *v_config) {
	_config = v_config;
}

Tmy_alpha::~Tmy_alpha() {
	clear_container();
	delete _my_list_alpha;
}

void Tmy_alpha::clear_container()
{
	_my_list_alpha->clear_container();
}


void Tmy_alpha::init(int jml_data)
{
	Tmy_double ub_v1 = _config->eps1 / (_config->V1 * jml_data);
	Tmy_double ub_v2 = _config->eps2 / (_config->V2 * jml_data);

	_my_list_alpha_v1 = new Tmy_list_alpha(jml_data, 0, ub_v1);
	_my_list_alpha_v1->init(_config->V1, _config->eps1,0);
	
	_my_list_alpha_v2 = new Tmy_list_alpha(jml_data, 0, ub_v2);
	_my_list_alpha_v2->init(_config->V2, _config->eps2,1);
	
	_my_list_alpha = new Tmy_list_alpha(jml_data, (-1.0 * ub_v2), ub_v1);
	
	for (int i = 0; i < jml_data; ++i)
	{
		Tmy_double alpha_v1= _my_list_alpha_v1->get_alpha(i);
      Tmy_double alpha_v2= _my_list_alpha_v2->get_alpha(i);
		Tmy_double diff = alpha_v1 - alpha_v2;
		_my_list_alpha->update_alpha(i, diff);
	}

}

void Tmy_alpha::update_alpha(int idx1, Tmy_double value1, int idx2, Tmy_double value2)
{
	_my_list_alpha->update_alpha(idx1, value1);
	_my_list_alpha->update_alpha(idx2, value2);
}

Tmy_list_alpha* Tmy_alpha::get_alpha()
{
	return _my_list_alpha;
}

Tmy_list_alpha* Tmy_alpha::get_alpha_v1()
{
	return _my_list_alpha_v1;
}

Tmy_list_alpha* Tmy_alpha::get_alpha_v2()
{
	return _my_list_alpha_v2;
}
