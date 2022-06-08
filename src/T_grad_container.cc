#include "T_grad_container.h"

T_grad_container::T_grad_container()
{
	// _cek_kkt = false;
	// _filter_delta = true;
	// _min_rho = false;
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

/* int T_grad_container::max(Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, callback_type f)
{
	Tmy_double gmax = -HUGE_VAL;
	int idx_max = -1;
	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double dec_F = dec(tmp_idx[i], rho1, rho2);
		Tmy_double obj_F = obj(tmp_idx[i], rho1, rho2);

		callback_param var_b;
		var_b.idx = tmp_idx[i];
		var_b.dec = dec_F;
		var_b.obj = obj_F;
		var_b.grad = _grad[tmp_idx[i]];

		callback_param var_a;
		bool is_pass = true;
		if (_cek_kkt)
		{
			is_pass = !_is_kkt[tmp_idx[i]];
		}
		//    cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
		if (is_pass)
		{
			is_pass = f(var_b, var_a, alpha);
		}

		if (is_pass)
		{
			Tmy_double abs_F = abs(_grad[tmp_idx[i]]);
			if (_min_rho)
			{
				abs_F = abs(obj_F);
			}

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

bool T_grad_container::delta_filter(int idx_b, int idx_a, vector<T_alpha_container> alpha, Tmy_double delta)
{
	bool is_pass = true;

	if (alpha[0].is_nol(idx_b))
	{
		if (alpha[0][idx_a] > 0.0)
		{
			is_pass = (delta > 0.0) and (delta <= alpha[0][idx_a]);
		}
		else
		{
			if (alpha[0][idx_a] < 0.0)
			{
				is_pass = (delta < 0.0) and (abs(delta) <= abs(alpha[0][idx_a]));
			}
		}
	}
	else
	{
		if (alpha[0][idx_b] < 0.0)
		{
			if (alpha[0][idx_a] < 0.0)
			{
				if (delta > 0.0)
				{
					is_pass = delta <= abs(alpha[0][idx_b]);
				}
				else
				{
					if (delta < 0.0)
					{
						is_pass = abs(delta) <= abs(alpha[0][idx_a]);
					}
				}
			}
			else
			{
				if (alpha[0][idx_a] > 0.0)
				{
					if (delta > 0.0)
					{
						is_pass = delta <= abs(alpha[0][idx_b]);
						if (is_pass)
						{
							is_pass = delta <= abs(alpha[0][idx_a]);
						}
					}
					else
					{
						if (delta < 0.0)
						{
							is_pass = abs(delta) <= abs(alpha[0][idx_a]);
						}
					}
				}
			}
		}
		else
		{
			if (alpha[0][idx_b] > 0.0)
			{
				if (alpha[0][idx_a] < 0.0)
				{
					if (delta > 0.0)
					{
						is_pass = delta <= abs(alpha[0][idx_a]);
					}
					else
					{
						if (delta < 0.0)
						{
							is_pass = abs(delta) <= alpha[0][idx_b];
							if (is_pass)
							{
								is_pass = abs(delta) <= abs(alpha[0][idx_a]);
							}
						}
					}
				}
				else
				{
					if (alpha[0][idx_a] > 0.0)
					{
						if (delta > 0.0)
						{
							is_pass = delta <= abs(alpha[0][idx_a]);
						}
						else
						{
							if (delta < 0.0)
							{
								is_pass = abs(delta) <= alpha[0][idx_b];
							}
						}
					}
				}
			}
		}
	}

	return is_pass;
}

int T_grad_container::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, callback_type f)
{
	Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL;
	int idx_max = -1;

	Tmy_double Gb = _grad[idx_b];
	Tmy_double dec_Fb = dec(idx_b, rho1, rho2);
	Tmy_double obj_Fb = obj(idx_b, rho1, rho2);

	callback_param var_b;
	var_b.idx = idx_b;
	var_b.dec = dec_Fb;
	var_b.obj = obj_Fb;
	var_b.grad = Gb;

	Tmy_double Fa = 0.0;
	Tmy_double Fb = 0.0;

	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Ga = _grad[tmp_idx[i]];
		Tmy_double dec_Fa = dec(tmp_idx[i], rho1, rho2);
		Tmy_double obj_Fa = obj(tmp_idx[i], rho1, rho2);

		callback_param var_a;
		var_a.idx = tmp_idx[i];
		var_a.dec = dec_Fa;
		var_a.obj = obj_Fa;
		var_a.grad = Ga;

		Fa = Ga;
		Fb = Gb;
		if (_min_rho)
		{
			Fa = obj_Fa;
			Fb = obj_Fb;
		}

		bool is_pass = true;
		if (_cek_kkt)
		{
			is_pass = !_is_kkt[tmp_idx[i]];
		}
		// cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
		if (is_pass)
		{
			is_pass = f(var_b, var_a, alpha);
		}

		if (is_pass)
		{
			Tmy_double abs_Fa = abs(Fa);
			if (abs_Fa >= gmax)
			{
				gmax = abs_Fa;
			}

			Tmy_double diff_F = Fb - Fa;
			Tmy_double abs_diff_F = abs(diff_F);
			// cout << " idx_a " << tmp_idx[i] << " abs_diff_obj " << abs_diff_obj << " gmax2 " << gmax2 << endl;
			if ((abs_diff_F >= gmax2) or !_min_rho)
			{
				Tmy_double diff = Ga - Gb;
				vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, tmp_idx[i]);
				Tmy_double delta = diff * hsl_eta[0];

				bool is_pass = true;
				if (_filter_delta)
				{
					is_pass = delta_filter(idx_b, tmp_idx[i], alpha, delta);
				}
				if (is_pass)
				{
					// cout << " idx_a " << tmp_idx[i] << " delta " << delta << endl;
					gmax2 = abs_diff_F;
					idx_max = tmp_idx[i];
					mv_idx(tmp_idx[i], 0);
				}
			}
		}
	}

	if (idx_max != -1)
	{
		cout << endl
			 << "selisih [" << abs(Fb)
			 << "," << gmax
			 << "," << (abs(Fb) - gmax)
			 << "]" << endl;
	}
	return idx_max;
}

int T_grad_container::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
	Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL;
	int idx_max = -1;

	Tmy_double Gb = _grad[idx_b];
	Tmy_double abs_Gb = abs(Gb);
	Tmy_double dec_Fb = dec(idx_b, rho1, rho2);
	Tmy_double obj_Fb = obj(idx_b, rho1, rho2);

	callback_param var_b;
	var_b.idx = idx_b;
	var_b.dec = dec_Fb;
	var_b.obj = obj_Fb;
	var_b.grad = Gb;

	Tmy_double Fa = 0.0;
	Tmy_double Fb = 0.0;

	vector<int> tmp_idx = _idx;
	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Ga = _grad[tmp_idx[i]];
		Tmy_double abs_Ga = abs(Ga);
		Tmy_double dec_Fa = dec(tmp_idx[i], rho1, rho2);
		Tmy_double obj_Fa = obj(tmp_idx[i], rho1, rho2);

		callback_param var_a;
		var_a.idx = tmp_idx[i];
		var_a.dec = dec_Fa;
		var_a.obj = obj_Fa;
		var_a.grad = Ga;

		Fa = Ga;
		Fb = Gb;
		if (_min_rho)
		{
			Fa = obj_Fa;
			Fb = obj_Fb;
		}

		bool is_pass = true;
		if (_cek_kkt)
		{
			is_pass = !_is_kkt[tmp_idx[i]];
		}
		//    cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
		if (is_pass)
		{
			is_pass = f(var_b, var_a, alpha, kernel, my_alpha);
		}

		if (is_pass)
		{

			Tmy_double abs_Fa = abs(Fa);
			if (abs_Fa >= gmax)
			{
				gmax = abs_Fa;
			}

			Tmy_double diff_F = Fb - Fa;
			Tmy_double abs_diff_F = abs(diff_F);
			// cout << " idx_a " << tmp_idx[i] << " abs_diff_obj " << abs_diff_obj << " gmax2 " << gmax2 << endl;
			if ((abs_diff_F >= gmax2) or !_min_rho)
			{
				Tmy_double diff = Ga - Gb;
				vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, tmp_idx[i]);
				Tmy_double delta = diff * hsl_eta[0];

				bool is_pass = true;
				if (_filter_delta)
				{
					is_pass = delta_filter(idx_b, tmp_idx[i], alpha, delta);
				}
				if (is_pass)
				{
					// cout << " idx_a " << tmp_idx[i] << " delta " << delta << endl;
					gmax2 = abs_diff_F;
					idx_max = tmp_idx[i];
					mv_idx(tmp_idx[i], 0);
				}
			}
		}
	}

	if (idx_max != -1)
	{
		cout << endl
			 << "selisih [" << abs(Fb)
			 << "," << gmax
			 << "," << (abs(Fb) - gmax)
			 << "]" << endl;
	}
	return idx_max;
}

int T_grad_container::cari(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
	// cout << " Cari 1 " << endl;
	int idx_a = -1;

	Tmy_double Gb = _grad[idx_b];
	Tmy_double abs_Gb = abs(Gb);
	Tmy_double dec_Fb = dec(idx_b, rho1, rho2);
	Tmy_double obj_Fb = obj(idx_b, rho1, rho2);

	callback_param var_b;
	var_b.idx = idx_b;
	var_b.dec = dec_Fb;
	var_b.obj = obj_Fb;
	var_b.grad = Gb;

	for (int i = 0; i < _idx.size(); ++i)
	{
		Tmy_double Ga = _grad[_idx[i]];
		Tmy_double abs_Ga = abs(Ga);
		Tmy_double dec_Fa = dec(_idx[i], rho1, rho2);
		Tmy_double obj_Fa = obj(_idx[i], rho1, rho2);

		callback_param var_a;
		var_a.idx = _idx[i];
		var_a.dec = dec_Fa;
		var_a.obj = obj_Fa;
		var_a.grad = Ga;

		bool is_pass = true;
		if (_cek_kkt)
		{
			is_pass = !_is_kkt[_idx[i]];
		}
		//    cout << " idx a " << _idx[i] << " is_pass " << is_pass << endl;
		if (is_pass)
		{
			is_pass = f(var_b, var_a, alpha, kernel, my_alpha);
		}

		if (is_pass)
		{
			Tmy_double diff = Ga - Gb;
			vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, _idx[i]);
			Tmy_double delta = diff * hsl_eta[0];

			bool is_pass = true;
			if (_filter_delta)
			{
				is_pass = delta_filter(idx_b, _idx[i], alpha, delta);
			}
			if (is_pass)
			{
				idx_a = _idx[i];
				break;
			}
		}
	}
	// if (idx_a != -1)
	// {
	// 	cout << endl
	// 		 << "[" << abs_Gb
	// 		 << "," << abs(_grad[idx_a])
	// 		 << "," << (abs_Gb + abs(_grad[idx_a]))
	// 		 << "]" << endl;
	// }
	return idx_a;
} */

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