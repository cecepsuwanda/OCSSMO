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
	_my_list_alpha_v1->init(_config->V1, _config->eps1, 0);

	_my_list_alpha_v2 = new Tmy_list_alpha(jml_data, 0, ub_v2);
	_my_list_alpha_v2->init(_config->V2, _config->eps2, 1);

	_my_list_alpha = new Tmy_list_alpha(jml_data, (-1.0 * ub_v2), ub_v1);

	for (int i = 0; i < jml_data; ++i)
	{
		Tmy_double alpha_v1 = _my_list_alpha_v1->get_alpha(i);
		Tmy_double alpha_v2 = _my_list_alpha_v2->get_alpha(i);
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

Treturn_is_pass_h Tmy_alpha::is_pass(int i, int j, Tmy_double delta, Tmy_double delta_v1, Tmy_double delta_v2)
{
	Treturn_is_pass tmp_v1 = _my_list_alpha_v1->is_pass(i, j, delta_v1);
	Treturn_is_pass tmp_v2 = _my_list_alpha_v2->is_pass(i, j, delta_v2);
	Treturn_is_pass tmp = _my_list_alpha->is_pass(i, j, delta);

	Tmy_double diff = 0.0;
	auto cek = [&diff](Tmy_double alpha_v1, Tmy_double alpha_v2, Tmy_double alpha) -> bool {
		diff = (alpha_v1 - alpha_v2);
		Tmy_double tmp = diff - alpha;
		return (tmp > -1e-3) and (tmp < 1e-3);
	};

	bool old_alpha_i_pass = cek(tmp_v1.alpha_i, tmp_v2.alpha_i, tmp.alpha_i) == false;
	bool old_alpha_j_pass = cek(tmp_v1.alpha_j, tmp_v2.alpha_j, tmp.alpha_j) == false;
	bool new_alpha_i_pass = cek(tmp_v1.new_alpha_i, tmp_v2.new_alpha_i, tmp.new_alpha_i) == false;
	bool new_alpha_j_pass = cek(tmp_v1.new_alpha_j, tmp_v2.new_alpha_j, tmp.new_alpha_j) == false;
	bool v1_zero_v2_notzero = ((tmp_v1.new_alpha_i == 0.0) and (tmp_v1.new_alpha_j == 0.0)) and ((tmp_v2.new_alpha_i != 0.0) or (tmp_v2.new_alpha_j != 0.0));
	bool v2_zero_v1_notzero = ((tmp_v2.new_alpha_i == 0.0) and (tmp_v2.new_alpha_j == 0.0)) and ((tmp_v1.new_alpha_i != 0.0) or (tmp_v1.new_alpha_j != 0.0));
	bool v1_zero_v2_zero = ((tmp_v1.new_alpha_i == 0.0) and (tmp_v1.new_alpha_j == 0.0)) and ((tmp_v2.new_alpha_i == 0.0) or (tmp_v2.new_alpha_j == 0.0));
	bool v1_notzero_v2_notzero = ((tmp_v1.new_alpha_i != 0.0) or (tmp_v1.new_alpha_j != 0.0)) and ((tmp_v2.new_alpha_i != 0.0) or (tmp_v2.new_alpha_j != 0.0));

	Treturn_is_pass_h hsl;
	hsl.is_pass = true;

	if (new_alpha_i_pass or new_alpha_j_pass)
	{
		if (v1_zero_v2_notzero)
		{

		} else {
			if (v2_zero_v1_notzero)
			{

			} else {
				if (v1_zero_v2_zero)
				{
					hsl.is_pass = false;
				} else {
					if (v1_notzero_v2_notzero)
					{

					}
				}
			}
		}
	}

	if ((tmp_v1.is_pass == false) and (tmp_v2.is_pass == false))
	{
		hsl.is_pass = false;
	} else {
		
	}

	hsl.alpha = tmp;
	hsl.alpha_v1 = tmp_v1;
	hsl.alpha_v2 = tmp_v2;
	return hsl;
}
