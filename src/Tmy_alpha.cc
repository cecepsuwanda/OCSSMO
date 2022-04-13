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

Treturn_is_pass_h Tmy_alpha::is_pass(int i, int j, Tmy_double delta, Tmy_double delta_v1, Tmy_double delta_v2, int flag)
{
	Treturn_is_pass tmp_v1 = _my_list_alpha_v1->is_pass(i, j, delta_v1);
	Treturn_is_pass tmp_v2 = _my_list_alpha_v2->is_pass(i, j, delta_v2);
	Treturn_is_pass tmp = _my_list_alpha->is_pass(i, j, delta);


	auto cek = [](Tmy_double alpha_v1, Tmy_double alpha_v2, Tmy_double alpha) -> bool {
		return  (((alpha_v1 - alpha_v2) - alpha) == 0.0);
	};

	auto hitung_new_alpha = [](Tmy_double alpha_i, Tmy_double alpha_j, Tmy_double diff) -> vector<Tmy_double> {
		vector<Tmy_double> hasil;
		hasil.assign(2, 0);

		alpha_i = alpha_i + diff;
		alpha_j = alpha_j - diff;

		hasil[0] = (alpha_i > 1.0) ? 1.0 : (double) alpha_i;
		hasil[0] = (alpha_i < 0.0) ? 0.0 : (double) alpha_i;
		hasil[1] = (alpha_j > 1.0) ? 1.0 : (double) alpha_j;
		hasil[1] = (alpha_j < 0.0) ? 0.0 : (double) alpha_j;
		return hasil;

	};

	bool old_alpha_i_pass = cek(tmp_v1.alpha_i, tmp_v2.alpha_i, tmp.alpha_i) == false;
	bool old_alpha_j_pass = cek(tmp_v1.alpha_j, tmp_v2.alpha_j, tmp.alpha_j) == false;
	bool new_alpha_i_pass = cek(tmp_v1.new_alpha_i, tmp_v2.new_alpha_i, tmp.new_alpha_i) == false;
	bool new_alpha_j_pass = cek(tmp_v1.new_alpha_j, tmp_v2.new_alpha_j, tmp.new_alpha_j) == false;

	Treturn_is_pass_h hsl;
	hsl.is_pass = true;

	if (new_alpha_i_pass or new_alpha_j_pass)
	{
		if (flag == 1) {
			// cout << " " << tmp_v1.is_pass << " " << tmp_v2.is_pass << " ";
			// cout << "v1[" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] v2[" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] v[" << tmp.alpha_i << "," << tmp.alpha_j << "] ";
			// cout << "new v1[" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] new v2[" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] new v[" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] ";
		}

		if ( (tmp_v1.alpha_i == 0.0) and (tmp_v1.alpha_j == 0.0) and (tmp_v2.alpha_i == 0.0) and (tmp_v2.alpha_j == 0.0) )
		{
			hsl.is_pass = false;
		} else {

			Tmy_double diff_i = tmp.new_alpha_i - tmp.alpha_i;
			Tmy_double diff_j = tmp.new_alpha_j - tmp.alpha_j;

			Tmy_double diff_v1_i = tmp_v1.new_alpha_i - tmp_v1.alpha_i;
			Tmy_double diff_v1_j = tmp_v1.new_alpha_j - tmp_v1.alpha_j;

			Tmy_double diff_v2_i = tmp_v2.new_alpha_i - tmp_v2.alpha_i;
			Tmy_double diff_v2_j = tmp_v2.new_alpha_j - tmp_v2.alpha_j;

			Tmy_double delta_i =  tmp.new_alpha_i - (tmp_v1.alpha_i - tmp_v2.alpha_i);
			Tmy_double delta_j =  tmp.new_alpha_j - (tmp_v1.alpha_j - tmp_v2.alpha_j);

			vector<Tmy_double> hasil = hitung_new_alpha(tmp_v1.alpha_i, tmp_v1.alpha_j, delta_i);

			bool cek_i = cek(hasil[0], tmp_v2.new_alpha_i, tmp.new_alpha_i);
			bool cek_j = cek(hasil[1], tmp_v2.new_alpha_j, tmp.new_alpha_j);

			if (cek_i and cek_j) {
				tmp_v1.new_alpha_i = hasil[0];
				tmp_v1.new_alpha_j = hasil[1];
			} else {
				hasil = hitung_new_alpha(tmp_v2.alpha_i, tmp_v2.alpha_j, delta_i);

				cek_i = cek(tmp_v1.new_alpha_i, hasil[0], tmp.new_alpha_i);
				cek_j = cek(tmp_v1.new_alpha_j, hasil[1], tmp.new_alpha_j);

				if (cek_i and cek_j) {
					tmp_v2.new_alpha_i = hasil[0];
					tmp_v2.new_alpha_j = hasil[1];
				} else {
					hasil = hitung_new_alpha(tmp_v1.alpha_i, tmp_v1.alpha_j, -1.0 * delta_i);

					bool cek_i = cek(hasil[0], tmp_v2.new_alpha_i, tmp.new_alpha_i);
					bool cek_j = cek(hasil[1], tmp_v2.new_alpha_j, tmp.new_alpha_j);

					if (cek_i and cek_j) {
						tmp_v1.new_alpha_i = hasil[0];
						tmp_v1.new_alpha_j = hasil[1];
					} else {

						hasil = hitung_new_alpha(tmp_v2.alpha_i, tmp_v2.alpha_j, -1.0 * delta_i);

						cek_i = cek(tmp_v1.new_alpha_i, hasil[0], tmp.new_alpha_i);
						cek_j = cek(tmp_v1.new_alpha_j, hasil[1], tmp.new_alpha_j);

						if (cek_i and cek_j) {
							tmp_v2.new_alpha_i = hasil[0];
							tmp_v2.new_alpha_j = hasil[1];
						} else {
							hsl.is_pass = false;
						}

					}

				}

			}

			//if (flag == 1) cout << " " << diff_i << " " << diff_j << " " << delta_i << " " << delta_j << " " << diff_v1_i << " " << diff_v1_j << " " << diff_v2_i << " " << diff_v2_j << " ";
		}
	}


	old_alpha_i_pass = cek(tmp_v1.alpha_i, tmp_v2.alpha_i, tmp.alpha_i) == false;
	old_alpha_j_pass = cek(tmp_v1.alpha_j, tmp_v2.alpha_j, tmp.alpha_j) == false;
	new_alpha_i_pass = cek(tmp_v1.new_alpha_i, tmp_v2.new_alpha_i, tmp.new_alpha_i) == false;
	new_alpha_j_pass = cek(tmp_v1.new_alpha_j, tmp_v2.new_alpha_j, tmp.new_alpha_j) == false;

	if (flag == 1) {
		if (old_alpha_i_pass)
			cout << " old " << tmp_v1.alpha_i << "-" << tmp_v2.alpha_i << "=" << tmp.alpha_i << " ";
		if (old_alpha_j_pass)
			cout << " old " << tmp_v1.alpha_j << "-" << tmp_v2.alpha_j << "=" << tmp.alpha_j << " ";
		if (new_alpha_i_pass)
			cout << " new_i " << tmp_v1.new_alpha_i << "-" << tmp_v2.new_alpha_i << "=" << tmp.new_alpha_i << " ";
		if (new_alpha_j_pass)
			cout << " new_j " << tmp_v1.new_alpha_j << "-" << tmp_v2.new_alpha_j << "=" << tmp.new_alpha_j << " ";
	}

	if ((tmp_v1.is_pass == false) and (tmp_v2.is_pass == false) and (tmp.is_pass == false))
	{
		hsl.is_pass = false;
	} else {
		double diff = tmp.new_alpha_i - tmp.alpha_i;
		if ((diff > -1e-5) and (diff < 1e-5))
		{
			hsl.is_pass = false;
		}
	}

	hsl.alpha = tmp;
	hsl.alpha_v1 = tmp_v1;
	hsl.alpha_v2 = tmp_v2;
	return hsl;
}
