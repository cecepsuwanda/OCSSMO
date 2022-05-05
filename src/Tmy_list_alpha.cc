#include "Tmy_list_alpha.h"



Tmy_list_alpha::Tmy_list_alpha(int v_jml_data, Tmy_double v_lb, Tmy_double v_ub) {
  _jml_data = v_jml_data;
  _lb = v_lb;
  _ub = v_ub;
  _jml_alpha = 0;
  _jml_alpha_n_sv = 0;
  _n_all_sv = 0;
  _n_sv = 0;

  _alpha.reserve(_jml_data);
  _alpha_status.reserve(_jml_data);
  _alpha.assign(_jml_data, 0.0);
  _alpha_status.assign(_jml_data, 0);

}

Tmy_list_alpha::~Tmy_list_alpha() {
  clear_container();
  _alpha_sv.clear();
}

void Tmy_list_alpha::clear_container()
{
  _alpha.clear();
  _alpha_status.clear();   
}


void Tmy_list_alpha::init(Tmy_double V, Tmy_double eps, int flag) {
  Tmy_double tmp = V * ((double)_jml_data);
  int jml = (int) tmp;

  
  if (flag == 0) {

    for (int idx = 0; idx < jml; idx++) {
      update_alpha(idx, _ub);
    }

    if (_jml_alpha < eps)
    {
      tmp = eps - _jml_alpha;
      update_alpha(jml, tmp);
      jml = jml + 1;
    }

    for (int idx = jml; idx < _jml_data; idx++) {
      update_alpha(idx, 0.0);
    }
  } else {
    if (flag == 1)
    {
      int alpha_idx = _jml_data - 1;
      int i = 0;
      while (i < jml)
      {
        update_alpha(alpha_idx, _ub);
        i++;
        alpha_idx--;
      }

      if (_jml_alpha < eps)
      {
        tmp = eps - _jml_alpha;
        update_alpha(alpha_idx, tmp);
        alpha_idx--;
      }

      while (alpha_idx > 0)
      {
        update_alpha(alpha_idx, 0.0);
        alpha_idx--;
      }


    }
  }

}

void Tmy_list_alpha::update_alpha(int idx, Tmy_double value)
{
  if (_alpha.at(idx) != 0.0) {
    _jml_alpha = _jml_alpha - _alpha.at(idx);
  }

  vector<bool> hasil = is_alpha_sv(idx);

  if (hasil[0] == true) {
    if (_n_all_sv != 0) {
      _n_all_sv = _n_all_sv - 1;
    }
  }

  if (hasil[1] == true) {
    if (_n_sv != 0) {
      _jml_alpha_n_sv = _jml_alpha_n_sv - _alpha.at(idx);
      _n_sv = _n_sv - 1;
    }
  }

  _alpha.at(idx) = value;
  _jml_alpha = _jml_alpha + _alpha.at(idx);

  hasil = is_alpha_sv(idx);

  if (hasil[0] == true) {
    _n_all_sv = _n_all_sv + 1;
  }

  if (hasil[1] == true) {
    _jml_alpha_n_sv = _jml_alpha_n_sv + _alpha.at(idx);
    _n_sv = _n_sv + 1;
  }

  update_alpha_sv(idx);
  update_alpha_status(idx);  
}

void Tmy_list_alpha::update_alpha_status(int idx)
{

  if (_alpha.at(idx) >= _ub)
  {
    //cout<<"Stat 1 "<<idx<<" "<<_alpha.at(idx)<<endl;
    _alpha_status[idx] = 1;
  } else {
    if (_alpha.at(idx) <= _lb)
    {
      //cout<<"Stat 0 "<<idx<<" "<<_alpha.at(idx)<<endl;
      _alpha_status[idx] = 0;
    } else {
      //cout<<"Free "<<idx<<endl;
      _alpha_status[idx] = 2;
    }
  }
}

void Tmy_list_alpha::update_alpha_sv(int idx)
{

  map<int, Tmy_double>::iterator it;
  it = _alpha_sv.find(idx);
  if (it != _alpha_sv.end())
    _alpha_sv.erase (it);

  vector<bool> hasil = is_alpha_sv(idx);
  if (hasil[0] == true) {
    _alpha_sv[idx] = _alpha.at(idx);
  }
}



vector<bool> Tmy_list_alpha::is_alpha_sv(int idx)
{
  vector<bool> tmp;
  tmp.push_back((_alpha.at(idx) >= _lb) and (_alpha.at(idx) <= _ub) and (_alpha.at(idx) != 0.0));
  tmp.push_back((_alpha.at(idx) > _lb) and (_alpha.at(idx) < _ub) and (_alpha.at(idx) != 0.0));
  tmp.push_back(_alpha.at(idx) == _ub);
  tmp.push_back(_alpha.at(idx) == _lb);
  return tmp;
}

bool Tmy_list_alpha::is_lower_bound(int idx)
{
  return (_alpha_status.at(idx) == 0);
}

bool Tmy_list_alpha::is_upper_bound(int idx)
{
  return (_alpha_status.at(idx) == 1);
}

bool Tmy_list_alpha::is_free(int idx)
{
  return (_alpha_status.at(idx) == 2);
}
bool Tmy_list_alpha::is_nol(int idx_i,int idx_j)
{
  return ((_alpha.at(idx_i)==0.0) and (_alpha.at(idx_j)==0.0));
}



Tmy_double Tmy_list_alpha::get_alpha(int idx)
{
  return _alpha.at(idx);
}

map<int, Tmy_double> Tmy_list_alpha::get_list_alpha_sv()
{
  return _alpha_sv;
}

Treturn_alpha_stat Tmy_list_alpha::get_stat()
{
  Treturn_alpha_stat tmp_stat;
  tmp_stat.jml_alpha = _jml_alpha;
  tmp_stat.n_all_sv  = _n_all_sv;
  tmp_stat.n_sv = _n_sv;
  tmp_stat.jml_alpha_n_sv = _jml_alpha_n_sv;
  return tmp_stat;
}

Tmy_double Tmy_list_alpha::get_ub()
{
  return _ub;
}

Tmy_double Tmy_list_alpha::get_lb()
{
  return _lb;
}
