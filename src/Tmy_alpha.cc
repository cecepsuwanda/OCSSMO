#include "Tmy_alpha.h"


Tmy_alpha::Tmy_alpha(Tconfig *v_config) {
	_config = v_config;
}

Tmy_alpha::~Tmy_alpha() {
	
}




void Tmy_alpha::init(int jml_data, T_alpha_container& alpha)
{
	Tmy_double ub_v1 = _config->eps1 / (_config->V1 * jml_data);
	Tmy_double ub_v2 = _config->eps2 / (_config->V2 * jml_data);

	alpha.boundaries((-1.0 * ub_v2), ub_v1);
	alpha.reserve(jml_data);
	alpha.assign(jml_data, 0.0);

	Tmy_double tmp = _config->V1 * ((double) jml_data);
	int jml = (int) tmp;

	for (int i = 0; i < jml; ++i)
	{
		alpha[i] = ub_v1;
	}

	Tmy_double jml_alpha = alpha.sum();

	if (jml_alpha < _config->eps1)
	{
		alpha[jml] = _config->eps1 - jml_alpha;
	}



	tmp = _config->V2 * ((double) jml_data);
	jml = (int) tmp;

	Tmy_double total = 0.0;
	for (int i = (jml_data - 1); i > ((jml_data - 1) - jml); --i)
	{
		if (alpha[i] != 0.0) {
			alpha[i] = alpha[i] - ub_v2;
		} else {
			alpha[i] = -1.0 * ub_v2;
		}

		total = total + ub_v2;
	}

	if (total < _config->eps2)
	{
		if (alpha[((jml_data - 1) - jml)] != 0.0)
		{
			alpha[((jml_data - 1) - jml)] = alpha[((jml_data - 1) - jml)] - (_config->eps2 - total);
		} else {
			alpha[((jml_data - 1) - jml)] = -1.0 * (_config->eps2 - total);
		}
	}

}



vector<Tmy_double> Tmy_alpha::calculateBoundaries(int i, int j, T_alpha_container alpha)
{
	Tmy_double t      = alpha[i] + alpha[j];
	Tmy_double diff   = t - alpha.ub();
	Tmy_double diff1  = t + abs(alpha.lb());
	vector<Tmy_double> hasil = {alpha.lb(), alpha.ub()};
	if (((alpha[i] <= alpha.ub()) and (alpha[i] >= alpha.lb())) and ((alpha[j] <= alpha.ub()) and (alpha[j] >= alpha.lb()))) {
		hasil = {max(diff, alpha.lb()), min(alpha.ub(), diff1)};
	}
	return hasil;
}

vector<Tmy_double> Tmy_alpha::limit_alpha(Tmy_double alpha_a, Tmy_double alpha_b, Tmy_double Low, Tmy_double High, int flag)
{
	vector<Tmy_double> hasil = {alpha_a, alpha_b};
	if (alpha_a > High) {
		if (flag == 1)
		{
			Tmy_double s = alpha_a - High;
			hasil[1] = alpha_b + s;
		}
		hasil[0] = High;
	} else {
		if (alpha_a < Low) {
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
	//cout<<" alpha_old "<< _alpha.at(i) <<" alpha_a_new "<< alpha_a_new <<endl;
	vector<Tmy_double> tmp = limit_alpha(alpha_a_new, 0, Low, High, 0);
	alpha_a_new = tmp[0];
	Tmy_double alpha_b_new = alpha[j] + (alpha[i] - alpha_a_new);
	// tmp = limit_alpha(alpha_b_new,alpha_a_new,_lb,_ub,1);
	// alpha_b_new = tmp[0];
	// alpha_a_new = tmp[1];
	return {alpha[i], alpha[j], alpha_a_new, alpha_b_new};
}

Treturn_is_pass Tmy_alpha::is_pass(int i, int j, Tmy_double delta, T_alpha_container alpha)
{
	Treturn_is_pass tmp;
	tmp.is_pass = false;
	tmp.alpha_i = alpha[i];
	tmp.alpha_j = alpha[j];
	tmp.new_alpha_i = alpha[i];
	tmp.new_alpha_j = alpha[j];

	if (i == j)
	{
		return tmp;
	} else {
		vector<Tmy_double> hsl = calculateBoundaries(i, j, alpha);
		Tmy_double Low = hsl[0], High = hsl[1];
		//cout <<"Low "<<Low<<" High "<<High<<endl;
		if (Low == High) {
			return tmp;
		} else {
			vector<Tmy_double> hsl = calculateNewAlpha(i, j, delta, Low, High, alpha);
			Tmy_double alpha_a_old = hsl[0], alpha_b_old = hsl[1], alpha_a_new = hsl[2], alpha_b_new = hsl[3];
			double diff = alpha_a_new - alpha_a_old;
			//abs(diff)<10e-5
			if (abs(diff) < 10e-5)
			{
				return tmp;
			} else {
				tmp.is_pass = true;
				tmp.alpha_i = alpha_a_old;
				tmp.alpha_j = alpha_b_old;
				tmp.new_alpha_i = alpha_a_new;
				tmp.new_alpha_j = alpha_b_new;
				// cout<<"alpha_a_new : "<<alpha_a_new<<" alpha_a_old : "<<alpha_a_old<<endl;
				// cout<<"alpha_b_new : "<<alpha_b_new<<" alpha_b_old : "<<alpha_b_old<<endl;
				return tmp;
			}
		}
	}


	return tmp;
}
