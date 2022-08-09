#include "Tmy_alpha.h"

Tmy_alpha::Tmy_alpha(Tconfig *v_config)
{
	_config = v_config;
	_split = false;
}

Tmy_alpha::~Tmy_alpha()
{
}

void Tmy_alpha::init(int jml_data, T_alpha_container &alpha, T_alpha_container &alpha_v1, T_alpha_container &alpha_v2)
{

	Tmy_double ub_v1 = _config->eps1 / (_config->V1 * jml_data);
	Tmy_double ub_v2 = _config->eps2 / (_config->V2 * jml_data);

	alpha.boundaries((-1.0 * ub_v2), ub_v1);
	alpha.reserve(jml_data);
	alpha.assign(jml_data, 0.0);

	alpha_v1.boundaries(0.0, ub_v1);
	alpha_v1.reserve(jml_data);
	alpha_v1.assign(jml_data, 0.0);

	alpha_v2.boundaries(0.0, ub_v2);
	alpha_v2.reserve(jml_data);
	alpha_v2.assign(jml_data, 0.0);

	// alpha_v1[0] = _config->eps1;
	// alpha_v2[jml_data - 1] = _config->eps2;

	Tmy_double tmp = _config->V1 * ((double)jml_data);
	int jml_v1 = (int)tmp;

	tmp = _config->V2 * ((double)jml_data);
	int jml_v2 = (int)tmp;

	if (_split)
	{
		if ((jml_v1 + jml_v2) > jml_data)
		{
			if (jml_v1 > jml_v2)
			{
				jml_v1 = jml_data - jml_v2;
			}
			else
			{
				jml_v2 = jml_data - jml_v1;
			}
		}
	}

	for (int i = 0; i < jml_v1; ++i)
	{
		alpha_v1[i] = ub_v1;
	}

	Tmy_double jml_alpha = alpha_v1.sum();

	if (jml_alpha < _config->eps1)
	{
		if (_split)
		{
			alpha_v1[jml_v1 - 1] = alpha_v1[jml_v1 - 1] + (_config->eps1 - jml_alpha);
		}
		else
		{
			alpha_v1[jml_v1] = _config->eps1 - jml_alpha;
		}
	}

	int jml_sv = alpha_v1.n_sv();

	if (jml_sv == 0)
	{
		alpha_v1[0] = ub_v1 / 4.0;
		alpha_v1[1] = alpha_v1[1] + (ub_v1 - (ub_v1 / 4.0));
	}

	int idx = 0;
	if (_split)
	{
		idx = jml_data - 1;
		int i = 0;
		while (i < jml_v2)
		{
			alpha_v2[idx] = ub_v2;
			idx--;
			i++;
		}
	}
	else
	{
		for (int i = 0; i < jml_v2; ++i)
		{
			alpha_v2[i] = ub_v2;
		}
	}

	jml_alpha = alpha_v2.sum();

	if (jml_alpha < _config->eps2)
	{
		if (_split)
		{
			alpha_v2[idx + 1] = alpha_v2[idx + 1] + (_config->eps2 - jml_alpha);
		}
		else
		{
			alpha_v2[jml_v2] = _config->eps2 - jml_alpha;
		}
	}

	if (_split)
	{
		jml_sv = alpha_v2.n_sv();

		if (jml_sv == 0)
		{
			alpha_v2[jml_data - 1] = ub_v2 / 4.0;
			alpha_v2[jml_data - 2] = alpha_v2[jml_data - 2] + (ub_v2 - (ub_v2 / 4.0));
		}
	}
	else
	{
		int jml_sv = alpha_v2.n_sv();

		if (jml_sv == 0)
		{
			alpha_v2[0] = ub_v2 / 4.0;
			alpha_v2[1] = alpha_v2[1] + (ub_v2 - (ub_v2 / 4.0));
		}
	}

	for (int i = 0; i < jml_data; ++i)
	{
		alpha[i] = alpha_v1[i] - alpha_v2[i];
	}
}

vector<Tmy_double> Tmy_alpha::calculateBoundaries(int i, int j, T_alpha_container alpha)
{
	Tmy_double t = alpha[i] + alpha[j];
	Tmy_double diff = 0.0;
	Tmy_double diff1 = 0.0;
	diff = t - alpha.ub();
	diff1 = t + abs(alpha.lb());
	vector<Tmy_double> hasil = {alpha.lb(), alpha.ub()};
	if (((alpha[i] <= alpha.ub()) and (alpha[i] >= alpha.lb())) and ((alpha[j] <= alpha.ub()) and (alpha[j] >= alpha.lb())))
	{
		hasil = {max(diff, alpha.lb()), min(alpha.ub(), diff1)};
	}
	return hasil;
}

vector<Tmy_double> Tmy_alpha::limit_alpha(Tmy_double alpha_a, Tmy_double alpha_b, Tmy_double Low, Tmy_double High, int flag)
{
	vector<Tmy_double> hasil = {alpha_a, alpha_b};
	if (alpha_a > High)
	{
		if (flag == 1)
		{
			Tmy_double s = alpha_a - High;
			hasil[1] = alpha_b + s;
		}
		hasil[0] = High;
	}
	else
	{
		if (alpha_a < Low)
		{
			if (flag == 1)
			{
				Tmy_double s = alpha_a - Low;
				hasil[1] = alpha_b + s;
			}
			hasil[0] = Low;
		}
	}
	return hasil;
}

vector<Tmy_double> Tmy_alpha::calculateNewAlpha(int i, int j, Tmy_double delta, Tmy_double Low, Tmy_double High, T_alpha_container alpha)
{
	Tmy_double alpha_a_new = alpha[i] + delta;
	vector<Tmy_double> tmp = limit_alpha(alpha_a_new, 0, Low, High, 0);
	alpha_a_new = tmp[0];
	Tmy_double alpha_b_new = alpha[j] + (alpha[i] - alpha_a_new);
	tmp = limit_alpha(alpha_b_new, alpha_a_new, alpha.lb(), alpha.ub(), 1);
	alpha_b_new = tmp[0];
	alpha_a_new = tmp[1];
	return {alpha[i], alpha[j], alpha_a_new, alpha_b_new};
}

bool Tmy_alpha::is_pass(Treturn_is_pass &v1, Treturn_is_pass &v2, Treturn_is_pass v)
{
	bool is_pass = false;
	Treturn_is_pass coba_v1 = v1;
	Treturn_is_pass coba_v2 = v2;

	string msg = "";

	if ((coba_v2.new_alpha_i == 0.0) and (coba_v2.new_alpha_j == 0.0))
	{
		msg = " masuk 1 ";
		coba_v1.swap();
		if ((coba_v1 - coba_v2) != v)
		{
			if (v.new_alpha_i > 0.0)
			{
				coba_v1.set(1, 0.0);
				if (coba_v1.new_alpha_i > v.new_alpha_i)
				{
					coba_v1.set(0, v.new_alpha_i);
				}
			}
		}
	}
	else
	{
		if ((coba_v1.new_alpha_i == 0.0) and (coba_v1.new_alpha_j == 0.0))
		{
			msg = " masuk 2 ";
			coba_v2.swap();
			if ((coba_v1 - coba_v2) != v)
			{
				if (v.new_alpha_i < 0.0)
				{
					coba_v2.set(1, 0.0);
					if (coba_v2.new_alpha_i > abs(v.new_alpha_i))
					{
						coba_v2.set(0, abs(v.new_alpha_i));
					}
				}
			}
		}
		else
		{
			msg = " masuk 3 ";
			coba_v1.set(1, 0.0);
			coba_v2.set(0, 0.0);
			if ((coba_v1 - coba_v2) != v)
			{
				msg = msg + " masuk 31 ";
				coba_v1.set(0, 0.0);
				coba_v2.set(1, 0.0);
			}

			if ((coba_v1 - coba_v2) != v)
			{
				msg = msg + " masuk 32 ";
				if (v.new_alpha_i > 0.0)
				{
					coba_v1.set(1, 0.0);
					coba_v2.set(0, 0.0);
					if (coba_v1.new_alpha_i > v.new_alpha_i)
					{
						coba_v1.set(0, v.new_alpha_i);
					}
				}
				else
				{
					if (v.new_alpha_i < 0.0)
					{
						coba_v1.set(0, 0.0);
						coba_v2.set(1, 0.0);
						if (coba_v2.new_alpha_i > abs(v.new_alpha_i))
						{
							coba_v2.set(0, abs(v.new_alpha_i));
						}
					}
				}
			}

			if ((coba_v1 - coba_v2) != v)
			{
				msg = msg + " masuk 33 ";
				if ((coba_v1.is_pass == true) and (coba_v2.is_pass == true))
				{
				}
				else
				{
					if ((coba_v1.is_pass == true) and (coba_v2.is_pass == false))
					{
						// coba_v1.set(1, 0.0);
						// Tmy_double tmp = abs(v.new_alpha_i + coba_v2.new_alpha_i);
						// if (coba_v1.new_alpha_i > tmp)
						// {
						// 	coba_v1.set(0, tmp);
						// }
					}
					else
					{
						if ((coba_v1.is_pass == false) and (coba_v2.is_pass == true))
						{
							// coba_v2.set(1, 0.0);
							// Tmy_double tmp = abs(v.new_alpha_i - coba_v1.new_alpha_i);
							// if (coba_v2.new_alpha_i > tmp)
							// {
							// 	coba_v2.set(0, tmp);
							// }
						}
					}
				}
			}
		}
	}

	if ((coba_v1.is_pass == true) or (coba_v2.is_pass == true))
	{
		if ((coba_v1 - coba_v2) != v)
		{
			msg = " masuk 4 ";
			coba_v1.reset();
			coba_v2.reset();
		}
	}

	if ((coba_v1 - coba_v2) == v)
	{
		v1 = coba_v1;
		v2 = coba_v2;
		is_pass = true;
		msg = msg + " solve !!!";
		// cout << msg << endl;
		//  cout << "v1 " << coba_v1.is_pass << " old [" << coba_v1.alpha_i << "," << coba_v1.alpha_j << "] new [" << coba_v1.new_alpha_i << "," << coba_v1.new_alpha_j << "] " << endl;
		//  cout << "v2 " << coba_v2.is_pass << " old [" << coba_v2.alpha_i << "," << coba_v2.alpha_j << "] new [" << coba_v2.new_alpha_i << "," << coba_v2.new_alpha_j << "] " << endl;
		//  cout << "v " << v.is_pass << " old [" << v.alpha_i << "," << v.alpha_j << "] new [" << v.new_alpha_i << "," << v.new_alpha_j << "] " << endl;
	}
	else
	{
		// cout << "v1 " << coba_v1.is_pass << " old [" << coba_v1.alpha_i << "," << coba_v1.alpha_j << "] new [" << coba_v1.new_alpha_i << "," << coba_v1.new_alpha_j << "] " << endl;
		// cout << "v2 " << coba_v2.is_pass << " old [" << coba_v2.alpha_i << "," << coba_v2.alpha_j << "] new [" << coba_v2.new_alpha_i << "," << coba_v2.new_alpha_j << "] " << endl;
		// cout << "v " << v.is_pass << " old [" << v.alpha_i << "," << v.alpha_j << "] new [" << v.new_alpha_i << "," << v.new_alpha_j << "] " << endl;
		v1.new_alpha_i = v1.alpha_i;
		v1.new_alpha_j = v1.alpha_j;
		v2.new_alpha_i = v2.alpha_i;
		v2.new_alpha_j = v2.alpha_j;
		msg = msg + " failed !!!";
		// cout << msg << endl;
	}

	return is_pass;
}

Treturn_is_pass Tmy_alpha::is_pass(int i, int j, Tmy_double delta, T_alpha_container alpha)
{
	Treturn_is_pass tmp;
	tmp.is_pass = false;
	tmp.alpha_i = alpha[i];
	tmp.alpha_j = alpha[j];
	tmp.new_alpha_i = alpha[i];
	tmp.new_alpha_j = alpha[j];
	tmp.lb = alpha.lb();
	tmp.ub = alpha.ub();
	tmp.b_new_alpha_i = alpha[i];
	tmp.b_new_alpha_j = alpha[j];

	if (i == j)
	{
		return tmp;
	}
	else
	{
		vector<Tmy_double> hsl = calculateBoundaries(i, j, alpha);
		Tmy_double Low = hsl[0], High = hsl[1];
		// cout <<"Low "<<Low<<" High "<<High<<endl;
		if (Low == High)
		{
			return tmp;
		}
		else
		{
			vector<Tmy_double> hsl = calculateNewAlpha(j, i, delta, Low, High, alpha);
			Tmy_double alpha_a_old = hsl[1], alpha_b_old = hsl[0], alpha_a_new = hsl[3], alpha_b_new = hsl[2];
			double diff = alpha_b_new - alpha_b_old;
			// abs(diff)<10e-5
			if (abs(diff) < 1e-5)
			{
				return tmp;
			}
			else
			{
				tmp.is_pass = true;
				tmp.alpha_i = alpha_a_old;
				tmp.alpha_j = alpha_b_old;
				tmp.new_alpha_i = alpha_a_new;
				tmp.new_alpha_j = alpha_b_new;
				tmp.b_new_alpha_i = alpha_a_new;
				tmp.b_new_alpha_j = alpha_b_new;
				// cout<<"alpha_a_new : "<<alpha_a_new<<" alpha_a_old : "<<alpha_a_old<<endl;
				// cout<<"alpha_b_new : "<<alpha_b_new<<" alpha_b_old : "<<alpha_b_old<<endl;
				return tmp;
			}
		}
	}

	return tmp;
}

Treturn_is_pass Tmy_alpha::is_pass(int i, int j, Tmy_double c1, Tmy_double c2, T_alpha_container alpha)
{
	Treturn_is_pass tmp;
	tmp.is_pass = false;
	tmp.alpha_i = alpha[i];
	tmp.alpha_j = alpha[j];
	tmp.new_alpha_i = alpha[i];
	tmp.new_alpha_j = alpha[j];
	tmp.lb = alpha.lb();
	tmp.ub = alpha.ub();
	tmp.b_new_alpha_i = alpha[i];
	tmp.b_new_alpha_j = alpha[j];

	if (i == j)
	{
		return tmp;
	}
	else
	{
		vector<Tmy_double> hsl = calculateBoundaries(i, j, alpha);
		Tmy_double Low = hsl[0], High = hsl[1];
		// cout <<"Low "<<Low<<" High "<<High<<endl;
		if (Low == High)
		{
			return tmp;
		}
		else
		{
			c2 = c2 * alpha[j];
			Tmy_double Lobj = c1 * Low * Low + c2 * Low;
			Tmy_double Hobj = c1 * High * High + c2 * High;

			Tmy_double delta = 0.0;
			if (Lobj > Hobj)
			{
				delta = Low - alpha[j];
			}
			else
			{
				if (Lobj < Hobj)
				{
					delta = High - alpha[j];
				}
			}

			vector<Tmy_double> hsl = calculateNewAlpha(j, i, delta, Low, High, alpha);
			Tmy_double alpha_a_old = hsl[1], alpha_b_old = hsl[0], alpha_a_new = hsl[3], alpha_b_new = hsl[2];
			double diff = alpha_b_new - alpha_b_old;
			// abs(diff)<10e-5
			if (abs(diff) < 1e-5)
			{
				return tmp;
			}
			else
			{
				tmp.is_pass = true;
				tmp.alpha_i = alpha_a_old;
				tmp.alpha_j = alpha_b_old;
				tmp.new_alpha_i = alpha_a_new;
				tmp.new_alpha_j = alpha_b_new;
				tmp.b_new_alpha_i = alpha_a_new;
				tmp.b_new_alpha_j = alpha_b_new;
				// cout<<"alpha_a_new : "<<alpha_a_new<<" alpha_a_old : "<<alpha_a_old<<endl;
				// cout<<"alpha_b_new : "<<alpha_b_new<<" alpha_b_old : "<<alpha_b_old<<endl;
				return tmp;
			}
		}
	}

	return tmp;
}