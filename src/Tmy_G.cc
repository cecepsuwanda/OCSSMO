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

Treturn_update_rho Tmy_G::update_rho()
{
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  Treturn_update_rho tmp_rho;

  tmp_rho.rho_v1 = 0.0;
  tmp_rho.rho_v2 = 0.0;

  vector<int> alpha_free_v1 = list_alpha_v1->get_list_lb_ub(2);
  vector<int> alpha_free_v2 = list_alpha_v2->get_list_lb_ub(2);

  int active_size = _my_list_G->get_active_size();

  int jml_free = 0;
  if (alpha_free_v1.size() > 0)
  {
    Tmy_double sum_free = 0.0;
    for (auto& it : alpha_free_v1)
    {
      if (it < active_size) {
        sum_free = sum_free + _my_list_G->get_G(it);
        jml_free = jml_free + 1;
      }
    }
    if (jml_free > 0)
      tmp_rho.rho_v1 = sum_free / (double) jml_free;
  }

  jml_free = 0;
  if (alpha_free_v2.size() > 0)
  {
    Tmy_double sum_free = 0.0;
    for (auto& it : alpha_free_v2)
    {
      if (it < active_size) {
        sum_free = sum_free + _my_list_G->get_G(it);
        jml_free = jml_free + 1;
      }
    }
    if (jml_free > 0)
      tmp_rho.rho_v2 = sum_free / (double) jml_free;
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

bool Tmy_G::is_kkt(int idx, Treturn_update_rho rho)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  Tmy_double F = (_my_list_G->get_G(idx) - rho.rho_v1) * (rho.rho_v2 - _my_list_G->get_G(idx));

  vector<bool> hsl = list_alpha->is_alpha_sv(idx);
  bool stat = false;


  if ((hsl[3] == true) and (F >= 1e-3))
  {
    stat = true;
  } else {
    if ((hsl[1] == true) and ((F > -1e-3) and (F < 1e-3)))
    {
      stat = true;
    } else {
      if ((hsl[2] == true) and (F <= 1e-3))
      {
        stat = true;
      }
    }
  }

  return stat;
}


Treturn_cari_idx Tmy_G::cari_idx()
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(0);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(1);

  int active_size = _my_list_G->get_active_size();

  int idx_b = -1, idx_a = -1;
  Tmy_double gmax = -HUGE_VAL, kkt_max = -HUGE_VAL, gmin = -HUGE_VAL, diff_max = -HUGE_VAL, obj_diff_min = HUGE_VAL;

  for (auto& idx : idx_alpha_not_ub)
  {
    if (idx < _active_size) {

      Tmy_double Fb = _my_list_G->get_G(idx);
      Tmy_double diff = Fb;
      Tmy_double tmp = -1.0 * diff;
      if (tmp >= gmax)
      {
        gmax = tmp;
        idx_b = idx;
        _my_list_G->mv_lb_ub(idx, 1);
      }

    }
  }

  if (idx_b != -1)
  {
    Tmy_double Fb = _my_list_G->get_G(idx_b);
    for (auto& idx : idx_alpha_not_lb)
    {
      if (idx < _active_size) {
        Tmy_double Fa = _my_list_G->get_G(idx);
        Tmy_double diff = Fa;
        if (diff >= gmin)
        {
          gmin = diff;
          _my_list_G->mv_lb_ub(idx, 0);
        }

        Tmy_double grad_diff = gmax + Fa;
        if (grad_diff > 1e-3)
        {
          _my_list_G->mv_lb_ub(idx, 0);
          Tmy_double obj_diff;
          vector<Tmy_double> tmp_hsl = _kernel->hit_eta(idx, idx_b, _active_size);
          Tmy_double quad_coef = tmp_hsl[0];
          obj_diff = -1.0 * (grad_diff * grad_diff) * quad_coef;
          if (obj_diff <= obj_diff_min)
          {
            obj_diff_min = obj_diff;
            idx_a = idx;
            _my_list_G->mv_lb_ub(idx, 0);
          }
        }
      }
    }
  }
  return {idx_b, idx_a, gmax, gmin};
}

int Tmy_G::cari_idx_lain(int idx_b)
{
  vector<int> idx_alpha_not_lb = _alpha->get_list_lb_ub(0);
  vector<int> idx_alpha_not_ub = _alpha->get_list_lb_ub(1);

  int idx_a = -1;

  for (auto& idx : idx_alpha_not_lb)
  {
    if (idx < _active_size) {
      Tmy_double Fb = _arr_G[idx_b];
      Tmy_double Fa = _arr_G[idx];
      Tmy_double gmax = -1.0 * Fb;
      Tmy_double gmin = Fa;
      Tmy_double diff = gmax + gmin;
      if (diff > 1e-3) {
        vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size);

        Tmy_double delta = (Fa - Fb) * hsl_eta[0];

        vector<Tmy_double> alpha;
        Treturn_is_pass tmp = _alpha->is_pass(idx_b, idx, delta);
        if (tmp.is_pass == true)
        {
          idx_a = idx;
          break;
        }
      }
    }
  }

  if (idx_a == -1)
  {
    for (auto& idx : idx_alpha_not_ub)
    {
      if (idx < _active_size) {
        Tmy_double Fb = _arr_G[idx_b];
        Tmy_double Fa = _arr_G[idx];
        Tmy_double gmax = -1.0 * Fb;
        Tmy_double gmin = Fa;
        Tmy_double diff = gmax + gmin;
        if (diff > 1e-3) {
          vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size);

          Tmy_double delta = (Fa - Fb) * hsl_eta[0];

          vector<Tmy_double> alpha;
          Treturn_is_pass tmp = _alpha->is_pass(idx_b, idx, delta);
          if (tmp.is_pass == true)
          {
            idx_a = idx;
            break;
          }
        }
      }
    }
  }


  return idx_a;

}