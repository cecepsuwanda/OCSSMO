#include "Tmy_G.h"


Tmy_G::Tmy_G(int jml_data, Tmy_kernel *kernel, Tmy_alpha *alphas) {
  _jml_data = jml_data;
  _kernel = kernel;
  _alphas = alphas;

}

Tmy_G::~Tmy_G() {
  clear_container();
  delete _my_list_G;
}

void Tmy_G::clear_container()
{
  _my_list_G->clear_container();
}

void Tmy_G::init()
{
  _my_list_G = new Tmy_list_G(_jml_data, _kernel, _alphas->get_alpha());
  _my_list_G_v1 = new Tmy_list_G(_jml_data, _kernel, _alphas->get_alpha_v1());
  _my_list_G_v2 = new Tmy_list_G(_jml_data, _kernel, _alphas->get_alpha_v2());
}

Tmy_double Tmy_G::update_rho()
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_double tmp_rho = 0.0;

  vector<int> alpha_free = list_alpha->get_list_lb_ub(2);
  int active_size = _my_list_G->get_active_size();

  int jml_free = 0;
  if (alpha_free.size() > 0)
  {
    Tmy_double sum_free = 0.0;
    for (auto& it : alpha_free)
    {
      if (it < active_size) {
        sum_free = sum_free + _my_list_G->get_G(it);
        jml_free = jml_free + 1;
      }
    }
    if (jml_free > 0)
      tmp_rho = sum_free / (double) jml_free;

  }

  if (jml_free == 0)
  {
    vector<int> alpha_not_lb = list_alpha->get_list_lb_ub(0);
    vector<int> alpha_not_ub = list_alpha->get_list_lb_ub(1);
    Tmy_double max_ub = -100, min_lb = 100;
    for (auto& it : alpha_not_lb)
    {
      if (it < active_size) {
        if (list_alpha->is_upper_bound(it) == true)
        {
          Tmy_double tmp = _my_list_G->get_G(it);
          if (tmp > max_ub)
          {
            max_ub = tmp;
          }
        }
      }
    }

    for (auto& it : alpha_not_ub)
    {
      if (it < active_size) {
        if (list_alpha->is_lower_bound(it) == true)
        {
          Tmy_double tmp = _my_list_G->get_G(it);
          if (tmp < min_lb)
          {
            min_lb = tmp;
          }
        }
      }
    }

    tmp_rho = (max_ub + min_lb) / 2.0;
  }
  return tmp_rho;
}

Tmy_list_G* Tmy_G::get_list_G()
{
  return _my_list_G;
}

Tmy_list_G* Tmy_G::get_list_G_v1()
{
  return _my_list_G_v1;
}

Tmy_list_G* Tmy_G::get_list_G_v2()
{
  return _my_list_G_v2;
}