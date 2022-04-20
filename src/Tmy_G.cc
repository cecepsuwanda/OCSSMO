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
  Tmy_list_alpha *list_alpha    = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  Treturn_update_rho tmp_rho;
  tmp_rho.rho_v1 = 0.0;
  tmp_rho.rho_v2 = 0.0;

  int jml_n_v1 = 0;
  Tmy_double jml_G_v1 = 0.0;
  int jml_n_v2 = 0;
  Tmy_double jml_G_v2 = 0.0;

  for (int i = 0; i < _jml_data; ++i)
  {
    vector<bool> is_sv_v1 = list_alpha_v1->is_alpha_sv(i);
    vector<bool> is_sv_v2 = list_alpha_v2->is_alpha_sv(i);

    if (is_sv_v1[0])
    {
      if (is_sv_v1[1])
      {
        jml_G_v1 = jml_G_v1 + _my_list_G->get_G(i);
      }
      jml_n_v1++;
    }

    if (is_sv_v2[0])
    {
      if (is_sv_v2[1])
      {
        jml_G_v2 = jml_G_v2 + _my_list_G->get_G(i);
      }
      jml_n_v2++;
    }
  }

  if (jml_n_v1 > 0)
  {
    tmp_rho.rho_v1 = jml_G_v1 / (1.0 * jml_n_v1);
  }

  if (jml_n_v2 > 0)
  {
    tmp_rho.rho_v2 = jml_G_v2 / (1.0 * jml_n_v2);
  }

  return tmp_rho;
}

Treturn_update_rho Tmy_G::update_rho(Treturn_update_rho old_rho)
{
  Tmy_list_alpha *list_alpha    = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  Treturn_update_rho tmp_rho;
  tmp_rho.rho_v1 = 0.0;
  tmp_rho.rho_v2 = 0.0;

  int jml_n_v1 = 0;
  Tmy_double jml_G_v1 = 0.0;
  int jml_n_v2 = 0;
  Tmy_double jml_G_v2 = 0.0;

  Tmy_double abs_F = 0.0;
  // if (_idx_exclude.size() > 0)
  // {
  //   Tmy_double F = _my_list_G->get_F(_idx_exclude[0], old_rho.rho_v1, old_rho.rho_v2);
  //   abs_F = abs(F);
  // }

  for (int i = 0; i < _jml_data; ++i)
  {
    bool is_pass = true;

    // if (_idx_exclude.size() > 0)
    // {
    //   Tmy_double F = _my_list_G->get_F(i, old_rho.rho_v1, old_rho.rho_v2);
    //   Tmy_double abs_F1 = abs(F);
    //   is_pass = abs_F1 < abs_F;
    // }

    if (is_pass) {
      vector<bool> is_sv_v1 = list_alpha_v1->is_alpha_sv(i);
      vector<bool> is_sv_v2 = list_alpha_v2->is_alpha_sv(i);

      if (is_sv_v1[0])
      {
        if (is_sv_v1[1])
        {
          jml_G_v1 = jml_G_v1 + _my_list_G->get_G(i);
        }
        jml_n_v1++;
      }

      if (is_sv_v2[0])
      {
        if (is_sv_v2[1])
        {
          jml_G_v2 = jml_G_v2 + _my_list_G->get_G(i);
        }
        jml_n_v2++;
      }
    }
  }

  if (jml_n_v1 > 0)
  {
    tmp_rho.rho_v1 = jml_G_v1 / (1.0 * jml_n_v1);
  }

  if (jml_n_v2 > 0)
  {
    tmp_rho.rho_v2 = jml_G_v2 / (1.0 * jml_n_v2);
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


  if ((hsl[3] == true) and (F > 0.0))
  {
    stat = true;
  } else {
    if ((hsl[1] == true) and (abs(F) < 1e-3))
    {
      stat = true;
    } else {
      if ((hsl[2] == true) and (F < 0.0))
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

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(0);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(1);

  int idx_b = -1, idx_a = -1;
  Tmy_double gmax = -HUGE_VAL, kkt_max = -HUGE_VAL, gmin = -HUGE_VAL, diff_max = -HUGE_VAL, obj_diff_min = HUGE_VAL;

  // for (int i = 0; i < _active_size; ++i)
  // {
  //   if (list_alpha->is_neg(i) == false)
  //   {
  //     if (list_alpha->is_upper_bound(i) == false)
  //     {
  //       Tmy_double Fb = _my_list_G->get_G(i);
  //       Tmy_double tmp = -1.0 * Fb;
  //       if (tmp >= gmax)
  //       {
  //         gmax = tmp;
  //         idx_b = i;
  //         list_alpha->mv_lb_ub(i, 1);
  //       }
  //     }
  //   } else {
  //     if (list_alpha->is_lower_bound(i) == false)
  //     {
  //       Tmy_double Fb = _my_list_G->get_G(i);
  //       Tmy_double tmp = Fb;
  //       if (tmp >= gmax)
  //       {
  //         gmax = tmp;
  //         idx_b = i;
  //         list_alpha->mv_lb_ub(i, 0);
  //       }
  //     }

  //   }
  // }

  // for (int i = 0; i < _active_size; ++i)
  // {
  //   if (list_alpha->is_neg(i) == false)
  //   {
  //     if (list_alpha->is_lower_bound(i) == false)
  //     {
  //       Tmy_double Fa = _my_list_G->get_G(i);
  //       Tmy_double tmp = Fa;
  //       if (tmp >= gmin)
  //       {
  //         gmin = tmp;
  //         list_alpha->mv_lb_ub(i, 0);
  //       }

  //       Tmy_double grad_diff = gmax + tmp;
  //       if (grad_diff > 0.0)
  //       {
  //         list_alpha->mv_lb_ub(i, 0);
  //         Tmy_double obj_diff;
  //         vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, i, _active_size, 0);
  //         obj_diff = -1.0 * (grad_diff * grad_diff) * hsl_eta[0];
  //         if (obj_diff <= obj_diff_min)
  //         {
  //           obj_diff_min = obj_diff;
  //           idx_a = i;
  //           list_alpha->mv_lb_ub(i, 0);
  //         }
  //       }
  //     }

  //   } else {
  //     if (list_alpha->is_upper_bound(i) == false)
  //     {
  //       Tmy_double Fa = _my_list_G->get_G(i);
  //       Tmy_double tmp = -1.0 * Fa;
  //       if (tmp >= gmin)
  //       {
  //         gmin = tmp;
  //         list_alpha->mv_lb_ub(i, 1);
  //       }

  //       Tmy_double grad_diff = gmax + tmp;
  //       if (grad_diff > 0.0)
  //       {
  //         list_alpha->mv_lb_ub(i, 1);
  //         Tmy_double obj_diff;
  //         vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, i, _active_size, 1);
  //         obj_diff = -1.0 * (grad_diff * grad_diff) * hsl_eta[0];
  //         if (obj_diff <= obj_diff_min)
  //         {
  //           obj_diff_min = obj_diff;
  //           idx_a = i;
  //           list_alpha->mv_lb_ub(i, 1);
  //         }
  //       }
  //     }
  //   }

  // }

  return {idx_b, idx_a, gmax, gmin};
}

Treturn_cari_idx Tmy_G::cari_idx(Treturn_update_rho rho)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(3);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(4);

  int idx_b = -1, idx_a = -1;
  Tmy_double gmax = -HUGE_VAL, gmin = -HUGE_VAL, gmax1 = -HUGE_VAL, gmax2 = -HUGE_VAL;

  for (auto& idx : idx_alpha_not_ub)
  {
    Tmy_double Fb = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_Fb = abs(Fb);
    Tmy_double diff = -1.0 * abs_Fb;
    if (diff >= gmax)
    {
      gmax = diff;
      idx_b = idx;
      list_alpha->mv_lb_ub(idx, 4);
    }
  }

  if (idx_b != -1)
  {
    Tmy_double Fb = _my_list_G->get_F(idx_b, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_Fb = abs(Fb);
    for (auto& idx : idx_alpha_not_lb)
    {
      Tmy_double Fa = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_Fa = abs(Fa);
      Tmy_double diff = abs_Fa;
      if (diff >= gmin)
      {
        gmin = diff;
        list_alpha->mv_lb_ub(idx, 3);
      }

      Tmy_double diff_F = gmax + diff;
      Tmy_double diff_abs_F = abs(Fb - Fa);
      if (abs(diff_F) > 1e-3)
      {
        list_alpha->mv_lb_ub(idx, 3);
        bool is_pass = true;
        is_pass = !list_alpha->is_nol(idx_b, idx);
        if ((diff_abs_F >= gmax2) and is_pass)
        {
          idx_a = idx;
          gmax2 = diff_abs_F;
          list_alpha->mv_lb_ub(idx, 3);
        }
      }


    }
  }

  return {idx_b, idx_a, gmax, gmin};
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  auto cek = [&, this](int idx_b, int idx_a) -> bool {

    vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx_a, _jml_data, 0);
    vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx_a, _jml_data);

    Tmy_double hsl_sum_v1 = this->sum_alpha_diff_Q(list_alpha_v1, hsl_diff);
    Tmy_double delta_v1 = hsl_eta[0] * hsl_sum_v1;

    Tmy_double hsl_sum_v2 = this->sum_alpha_diff_Q(list_alpha_v2, hsl_diff);
    Tmy_double delta_v2 = hsl_eta[0] * hsl_sum_v2;

    Tmy_double hsl_sum = this->sum_alpha_diff_Q(list_alpha, hsl_diff);
    Tmy_double delta = hsl_eta[0] * hsl_sum;

    Treturn_is_pass_h tmp = _alphas->is_pass(idx_b, idx_a, delta, delta_v1, delta_v2, 0);

    return  (tmp.is_pass);
  };


  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(3);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(4);

  Tmy_double Fb = _my_list_G->get_F(idx_b, rho.rho_v1, rho.rho_v2);
  Tmy_double abs_Fb = abs(Fb);

  Tmy_double gmax2 = -HUGE_VAL;
  int idx_a = -1;
  for (auto& idx : idx_alpha_not_lb)
  {
    Tmy_double Fa = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_Fa = abs(Fa);

    Tmy_double diff_F = (-1.0 * abs_Fb) + abs_Fa ;
    Tmy_double diff_abs_F = abs(Fb - Fa);

    if (abs(diff_F) > 1e-3)
    {
      bool is_pass = true;
      is_pass = !list_alpha->is_nol(idx_b, idx);
      if ((diff_abs_F >= gmax2) and is_pass)
      {
        if (cek(idx_b, idx)) {
          idx_a = idx;
          gmax2  = diff_abs_F;
          list_alpha->mv_lb_ub(idx, 3);
        }
      }
    }
  }

  if (idx_a == -1)
  {
    for (auto& idx : idx_alpha_not_lb)
    {
      Tmy_double Fa = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_Fa = abs(Fa);

      Tmy_double diff_F = (-1.0 * abs_Fb) + abs_Fa ;
      Tmy_double diff_abs_F = abs(Fb - Fa);

      if (abs(diff_F) > 1e-3)
      {
        bool is_pass = true;
        is_pass = !list_alpha->is_nol(idx_b, idx);
        if (cek(idx_b, idx) and is_pass)
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
      Tmy_double Fa = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_Fa = abs(Fa);

      Tmy_double diff_F = (-1.0 * abs_Fb) + abs_Fa ;
      Tmy_double diff_abs_F = abs(Fb - Fa);

      if (abs(diff_F) > 1e-3)
      {
        bool is_pass = true;
        is_pass = !list_alpha->is_nol(idx_b, idx);
        if ((diff_abs_F >= gmax2) and is_pass)
        {
          if (cek(idx_b, idx)) {
            idx_a = idx;
            gmax2  = diff_abs_F;
            list_alpha->mv_lb_ub(idx, 4);
          }
        }
      }
    }
  }

  if (idx_a == -1)
  {
    for (auto& idx : idx_alpha_not_ub)
    {
      Tmy_double Fa = _my_list_G->get_F(idx, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_Fa = abs(Fa);

      Tmy_double diff_F = (-1.0 * abs_Fb) + abs_Fa ;
      Tmy_double diff_abs_F = abs(Fb - Fa);

      if (abs(diff_F) > 1e-3)
      {
        bool is_pass = true;
        is_pass = !list_alpha->is_nol(idx_b, idx);
        if (cek(idx_b, idx) and is_pass)
        {
          idx_a = idx;
          break;
        }
      }
    }
  }

  return idx_a;

}

int Tmy_G::cari_idx_lain(int idx_b)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  vector<int> idx_alpha_not_lb = list_alpha->get_list_lb_ub(0);
  vector<int> idx_alpha_not_ub = list_alpha->get_list_lb_ub(1);

  int idx_a = -1;

  // for (auto& idx : idx_alpha_not_lb)
  // {
  //   if (idx < _active_size) {
  //     Tmy_double Fb = _my_list_G->get_G(idx_b);
  //     Tmy_double Fa = _my_list_G->get_G(idx);

  //     Tmy_double gmax = Fb;
  //     if (list_alpha->is_neg(idx) == false)
  //     {
  //       gmax = -1.0 * Fb;
  //     }

  //     Tmy_double gmin = Fa;
  //     if (list_alpha->is_neg(idx) == true)
  //     {
  //       gmin = -1.0 * Fa;
  //     }

  //     Tmy_double diff = gmax + gmin;
  //     if (diff > 0.0) {
  //       vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size, 0);
  //       vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx, _active_size);

  //       Tmy_double hsl_sum = sum_alpha_diff_Q(list_alpha, hsl_diff);
  //       Tmy_double delta = hsl_eta[0] * hsl_sum;

  //       Treturn_is_pass tmp = list_alpha->is_pass(idx_b, idx, delta);
  //       if (tmp.is_pass == true)
  //       {
  //         idx_a = idx;
  //         break;
  //       }
  //     }
  //   }

  // }

  if (idx_a == -1)
  {
    // for (auto& idx : idx_alpha_not_ub)
    // {
    //   if (idx < _active_size) {
    //     Tmy_double Fb = _my_list_G->get_G(idx_b);
    //     Tmy_double Fa = _my_list_G->get_G(idx);

    //     Tmy_double gmax = Fb;
    //     if (list_alpha->is_neg(idx) == false)
    //     {
    //       gmax = -1.0 * Fb;
    //     }

    //     Tmy_double gmin = Fa;
    //     if (list_alpha->is_neg(idx) == true)
    //     {
    //       gmin = -1.0 * Fa;
    //     }

    //     Tmy_double diff = gmax + gmin;
    //     if (diff > 0.0) {
    //       vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx, _active_size, 0);
    //       vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx, _active_size);

    //       Tmy_double hsl_sum = sum_alpha_diff_Q(list_alpha, hsl_diff);
    //       Tmy_double delta = hsl_eta[0] * hsl_sum;

    //       Treturn_is_pass tmp = list_alpha->is_pass(idx_b, idx, delta);
    //       if (tmp.is_pass == true)
    //       {
    //         idx_a = idx;
    //         break;
    //       }
    //     }
    //   }

    // }
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
    //list_alpha->is_upper_bound(i)==false
    if (list_alpha->is_not_lb_ub(i) == true)
    {
      Tmy_double diff = _my_list_G->get_G(i);
      if (list_alpha->is_neg(i) == false)
      {
        diff = -1.0 * diff;
      }

      if (diff >= gmax1)
      {
        gmax1 = diff;
      }

    }

    //list_alpha->is_lower_bound(i)==false
    if (list_alpha->is_nol(i) == false)
    {
      Tmy_double diff = _my_list_G->get_G(i);
      if (list_alpha->is_neg(i) == true)
      {
        diff = -1.0 * diff;
      }

      if (diff >= gmax2)
      {
        gmax2 = diff;
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

void Tmy_G::do_shrinking(Treturn_update_rho rho)
{
  cout << "do_shrinking" << endl;
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  Tmy_double gmax1 = -HUGE_VAL, gmax2 = -HUGE_VAL;

  for (int i = 0; i < _active_size; ++i)
  {
    //list_alpha->is_upper_bound(i)==false
    if (list_alpha->is_not_lb_ub(i) == true)
    {
      Tmy_double Fb = _my_list_G->get_G(i);
      Tmy_double diff = min(Fb - rho.rho_v1, rho.rho_v2 - Fb);

      if ((-1.0 * diff) >= gmax1)
      {
        gmax1 = (-1.0 * _my_list_G->get_G(i));
      }
    }

    //list_alpha->is_lower_bound(i)==false
    if (list_alpha->is_nol(i) == false)
    {
      Tmy_double Fa = _my_list_G->get_G(i);
      Tmy_double diff = min(Fa - rho.rho_v1, rho.rho_v2 - Fa);
      if (diff >= gmax2)
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

  //cout << "_active_size " << _active_size << " _jml_data " << _jml_data << endl;
}

bool Tmy_G::be_shrunk(int i, Tmy_double gmax1, Tmy_double gmax2)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  //list_alpha->is_upper_bound(i)
  if (list_alpha->is_not_lb_ub(i) == false)
  {
    Tmy_double diff = _my_list_G->get_G(i);
    if (list_alpha->is_neg(i) == false)
    {
      diff = -1.0 * diff;
    }

    return (diff > gmax1);

  }//list_alpha->is_lower_bound(i)
  else if (list_alpha->is_nol(i))
  {

    Tmy_double diff = _my_list_G->get_G(i);
    if (list_alpha->is_neg(i) == true)
    {
      diff = -1.0 * diff;
    }

    return (diff > gmax2);

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

void Tmy_G::update_G(int idx_b, int idx_a, Treturn_is_pass_h tmp, Treturn_update_rho rho)
{
  // int idx_exclude = -1;
  // if (_idx_exclude.size() > 0)
  // {
  //   idx_exclude = _idx_exclude[0];
  // }

  _my_list_G->update_G(idx_b, idx_a, tmp.alpha.new_alpha_i, tmp.alpha.new_alpha_j);
  _my_list_G_v1->update_G(idx_b, idx_a, tmp.alpha_v1.new_alpha_i, tmp.alpha_v1.new_alpha_j);
  _my_list_G_v2->update_G(idx_b, idx_a, tmp.alpha_v2.new_alpha_i, tmp.alpha_v2.new_alpha_j);
}

void Tmy_G::insert_idx_exclude(int idx)
{
  //_idx_exclude.insert(_idx_exclude.begin(), idx);
}

void Tmy_G::delete_idx_exclude()
{
  //_idx_exclude.clear();
}

bool Tmy_G::is_include(Tmy_double abs_F, Treturn_update_rho rho)
{
  bool is_pass = true;

  // if (_idx_exclude.size() > 0)
  // {
  //   Tmy_double tmp_F =  _my_list_G->get_F(_idx_exclude[0], rho.rho_v1, rho.rho_v2);
  //   is_pass = abs_F < abs(tmp_F);
  // }

  return is_pass;
}