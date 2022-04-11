#include "Tmy_G.h"


Tmy_G::Tmy_G(int jml_data, Tmy_kernel *kernel, Tmy_alpha *alphas) {
  _jml_data = jml_data;
  _active_size = _jml_data;
  _kernel = kernel;
  _alphas = alphas;
  _unshrink = false;
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



  int jml_free = 0;
  if (alpha_free_v1.size() > 0)
  {
    Tmy_double sum_free = 0.0;
    for (auto& it : alpha_free_v1)
    {
      if (it < _active_size) {
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
      if (it < _active_size) {
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
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(3);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(4);

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
        list_alpha->mv_lb_ub(idx, 4);
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
          list_alpha->mv_lb_ub(idx, 3);
        }

        Tmy_double grad_diff = gmax + Fa;
        if (grad_diff > 1e-3)
        {
          list_alpha->mv_lb_ub(idx, 3);
          Tmy_double obj_diff;
          vector<Tmy_double> tmp_hsl = _kernel->hit_eta(idx, idx_b, _active_size);
          Tmy_double quad_coef = tmp_hsl[0];
          obj_diff = -1.0 * (grad_diff * grad_diff) * quad_coef;
          if (obj_diff <= obj_diff_min)
          {
            obj_diff_min = obj_diff;
            idx_a = idx;
            list_alpha->mv_lb_ub(idx, 3);
          }
        }
      }
    }
  }
  return {idx_b, idx_a, gmax, gmin};
}

int Tmy_G::cari_idx_lain(int idx_b)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(3);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(4);

  int idx_a = -1;

  for (auto& idx : idx_alpha_not_lb)
  {
    if (idx < _active_size) {
      Tmy_double Fb = _my_list_G->get_G(idx_b);
      Tmy_double Fa = _my_list_G->get_G(idx);
      Tmy_double gmax = -1.0 * Fb;
      Tmy_double gmin = Fa;
      Tmy_double diff = gmax + gmin;
      if (diff > 1e-3) {
        vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size);
        vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx, _active_size);

        Tmy_double hsl_sum_v1 = sum_alpha_diff_Q(list_alpha_v1, hsl_diff);
        Tmy_double delta_v1 = hsl_eta[0] * hsl_sum_v1;

        Tmy_double hsl_sum_v2 = sum_alpha_diff_Q(list_alpha_v2, hsl_diff);
        Tmy_double delta_v2 = hsl_eta[0] * hsl_sum_v2;

        Tmy_double hsl_sum = sum_alpha_diff_Q(list_alpha, hsl_diff);
        Tmy_double delta = hsl_eta[0] * hsl_sum;

        Treturn_is_pass_h tmp = _alphas->is_pass(idx_b, idx, delta, delta_v1, delta_v2);
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
        Tmy_double Fb = _my_list_G->get_G(idx_b);
        Tmy_double Fa = _my_list_G->get_G(idx);
        Tmy_double gmax = -1.0 * Fb;
        Tmy_double gmin = Fa;
        Tmy_double diff = gmax + gmin;
        if (diff > 1e-3) {
          vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size);
          vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx, _active_size);

          Tmy_double hsl_sum_v1 = sum_alpha_diff_Q(list_alpha_v1, hsl_diff);
          Tmy_double delta_v1 = hsl_eta[0] * hsl_sum_v1;

          Tmy_double hsl_sum_v2 = sum_alpha_diff_Q(list_alpha_v2, hsl_diff);
          Tmy_double delta_v2 = hsl_eta[0] * hsl_sum_v2;

          Tmy_double hsl_sum = sum_alpha_diff_Q(list_alpha, hsl_diff);
          Tmy_double delta = hsl_eta[0] * hsl_sum;

          Treturn_is_pass_h tmp = _alphas->is_pass(idx_b, idx, delta, delta_v1, delta_v2);
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

Tmy_double Tmy_G::sum_alpha_diff_Q(Tmy_list_alpha* alpha, vector<Tmy_double> diff_Q)
{
  Tmy_double hasil = 0.0;

  for (int i = 0; i < diff_Q.size(); ++i)
  {
    hasil = hasil + (alpha->get_alpha(i) * diff_Q[i]);
  }

  return hasil;

}

void Tmy_G::do_shrinking()
{
  cout << "do_shrinking" << endl;
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  Tmy_double gmax1 = -HUGE_VAL, gmax2 = -HUGE_VAL;

  for (int i = 0; i < _active_size; ++i)
  {
    if (list_alpha->is_not_lb_ub(i) == true)
    {
      if ((-1.0 * _my_list_G->get_G(i)) >= gmax1)
      {
        gmax1 = (-1.0 * _my_list_G->get_G(i));
      }
    }

    if (list_alpha->is_nol(i) == false)
    {
      if (_my_list_G->get_G(i) >= gmax2)
      {
        gmax2 = _my_list_G->get_G(i);
      }
    }
  }

  if (_unshrink == false && ((gmax1 + gmax2) <= (1e-3 * 10)))
  {
    cout << "un shrink" << endl;
    _unshrink = true;
    reconstruct_gradient();
    _active_size = _jml_data;
  }

  for (int i = 0; i < _active_size; i++)
  {
    if (be_shrunk(i, gmax1, gmax2))
    {
      _active_size--;
      while (_active_size > i)
      {
        if (!be_shrunk(_active_size, gmax1, gmax2))
        {
          swap_index(i, _active_size);
          break;
        }
        _active_size--;
      }
    }
  }

  _my_list_G->set_active_size(_active_size);
  _my_list_G_v1->set_active_size(_active_size);
  _my_list_G_v2->set_active_size(_active_size);

  cout << "_active_size " << _active_size << " _jml_data " << _jml_data << endl;
}

bool Tmy_G::be_shrunk(int i, Tmy_double gmax1, Tmy_double gmax2)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  if (list_alpha->is_not_lb_ub(i) == false)
  {
    return ((-1.0 * _my_list_G->get_G(i)) > gmax1);
  }
  else if (list_alpha->is_nol(i) == true)
  {
    return (_my_list_G->get_G(i) > gmax2);
  }
  else
    return (false);
}

void Tmy_G::swap_index(int i, int j)
{
  _kernel->swap_index(i, j);
  _my_list_G->swap_index(i, j);
  _my_list_G_v1->swap_index(i, j);
  _my_list_G_v2->swap_index(i, j);

}

void Tmy_G::reconstruct_gradient()
{
  _my_list_G->reconstruct_gradient();
  _my_list_G_v1->reconstruct_gradient();
  _my_list_G_v2->reconstruct_gradient();
}

int Tmy_G::get_active_size()
{
  return _active_size;
}

void Tmy_G::reset_active_size()
{
  _active_size = _jml_data;
  _my_list_G->reset_active_size();
  _my_list_G_v1->reset_active_size();
  _my_list_G_v2->reset_active_size();
}

void Tmy_G::reverse_swap()
{
  _my_list_G->reverse_swap();
  _my_list_G_v1->reverse_swap();
  _my_list_G_v2->reverse_swap();
}