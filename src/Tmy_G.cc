#include "Tmy_G.h"

Tmy_G::Tmy_G()
{
  _cek_kkt = false;
  _filter_delta = false;
  _min_rho = false;
}

Tmy_G::~Tmy_G()
{
}

void Tmy_G::init(int jml_data, Tmy_kernel *kernel, T_alpha_container alpha, T_grad_container &grad)
{
  _jml_data = jml_data;

  grad.reserve(jml_data);
  grad.assign(jml_data, 0.0);

  for (int i = 0; i < _jml_data; ++i)
  {
    if (alpha[i] != 0.0)
    {
      vector<Tmy_double> data = kernel->get_Q(i);
      for (int j = 0; j < _jml_data; ++j)
      {
        grad[j] = grad[j] + (alpha[i] * data[j]);
      }
    }
  }
}

Treturn_update_rho Tmy_G::update_rho(Tmy_kernel *kernel, vector<T_alpha_container> alpha, T_grad_container grad)
{

  Treturn_update_rho tmp_rho;
  tmp_rho.rho_v1 = 0.0;
  tmp_rho.rho_v2 = 0.0;

  int jml_n_v1 = 0;
  Tmy_double jml_G_v1 = 0.0;
  // Tmy_double jml_G_v1_1 = 0.0;
  int jml_n_v2 = 0;
  Tmy_double jml_G_v2 = 0.0;
  // Tmy_double jml_G_v2_1 = 0.0;

  for (int i = 0; i < _jml_data; ++i)
  {
    if (alpha[2].is_nol(i))
    {
      if (alpha[1].is_sv(i))
      {
        // jml_G_v1 = jml_G_v1 + sum_alpha_rho_Q(i, kernel, alpha);
        jml_G_v1 = jml_G_v1 + grad[i];
        jml_n_v1++;
      }
    }

    if (alpha[1].is_nol(i))
    {
      if (alpha[2].is_sv(i))
      {
        // jml_G_v2 = jml_G_v2 + sum_alpha_rho_Q(i, kernel, alpha);
        jml_G_v2 = jml_G_v2 + grad[i];
        jml_n_v2++;
      }
    }
  }

  // cout << " rho [" << jml_n_v1 << "," << jml_n_v2 << "," << jml_G_v1 << "," << jml_G_v2 << "]" << endl;
  // cout << " rho [" << jml_G_v1_1 << "," << jml_G_v2_1 << "]" << endl;

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

bool Tmy_G::is_kkt(int idx, Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container grad)
{
  Tmy_double F = grad.dec(idx, rho.rho_v1, rho.rho_v2);
  bool stat = false;

  // if (alpha[1].is_nol(idx) and alpha[2].is_nol(idx))
  // {
  //   stat = (grad[idx] > rho.rho_v1) and (grad[idx] < rho.rho_v2);
  // }

  // if (alpha[1].is_sv(idx) and alpha[2].is_nol(idx))
  // {
  //   stat = (grad[idx] == rho.rho_v1) and (grad[idx] < rho.rho_v2);
  // }

  // if (alpha[1].is_nol(idx) and alpha[2].is_sv(idx))
  // {
  //   stat = (grad[idx] == rho.rho_v2);
  // }

  // if (alpha[1].is_ub(idx) and alpha[2].is_nol(idx))
  // {
  //   stat = (grad[idx] < rho.rho_v1) and (grad[idx] < rho.rho_v2);
  // }

  // if (alpha[1].is_nol(idx) and alpha[2].is_ub(idx))
  // {
  //   stat = (grad[idx] > rho.rho_v1) and (grad[idx] > rho.rho_v2);
  // }

  // if ((alpha[1].is_nol(idx) and alpha[2].is_nol(idx)) and (F > 0.0))
  // {
  //   stat = true;
  // }
  // else
  // {
  //   if (((alpha[1].is_sv(idx) and alpha[2].is_nol(idx)) or (alpha[1].is_nol(idx) and alpha[2].is_sv(idx))) and (F == 0.0))
  //   {
  //     stat = true;
  //   }
  //   else
  //   {
  //     if (((alpha[1].is_ub(idx) and alpha[2].is_nol(idx)) or (alpha[1].is_nol(idx) and alpha[2].is_ub(idx))) and (F < 0.0))
  //     {
  //       stat = true;
  //     }
  //   }
  // }

  if ((alpha[0].is_nol(idx)) and (F > 0.0))
  {
    stat = true;
  }
  else
  {
    if ((alpha[0].is_sv(idx)) and (F == 0.0))
    {
      stat = true;
    }
    else
    {
      if (((alpha[0].is_ub(idx)) or (alpha[0].is_lb(idx))) and (F < 0.0))
      {
        stat = true;
      }
    }
  }

  return stat;
}

void Tmy_G::set_kkt(Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container &grad)
{
  for (int i = 0; i < _jml_data; ++i)
  {
    // if (alpha[0].is_sv(i))
    // {
    grad.set_kkt(i, is_kkt(i, rho, alpha, grad));
    // }
    // else
    // {
    //   grad.set_kkt(i, false);
    // }
  }
}

int Tmy_G::cari_idx_a(int idx_b, Treturn_update_rho rho, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  {
    bool is_pass = true;
    // is_pass = alpha[1].is_sv(var_a.idx) or alpha[2].is_sv(var_a.idx);
    return is_pass;
  };

  int idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, cek_filter);

  return idx_a;
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho, Tmy_kernel *kernel, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_alpha *my_alpha)
{

  int idx_a = -1;

  auto cek1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (is_pass)
    {

      is_pass = alpha[1].is_sv(var_a.idx) or alpha[2].is_sv(var_a.idx);

      if (is_pass)
      {

        vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);

        Tmy_double delta_v1 = hsl_eta[0] * (grad[1][var_a.idx] - grad[1][var_b.idx]);
        Treturn_is_pass tmp_v1 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v1, alpha[1]);

        Tmy_double delta_v2 = hsl_eta[0] * (grad[2][var_a.idx] - grad[2][var_b.idx]);
        Treturn_is_pass tmp_v2 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v2, alpha[2]);

        Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
        // cout << "delta " << delta << endl;
        Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

        is_pass = (tmp_v1.is_pass or tmp_v2.is_pass);

        // if (is_pass)
        // {
        //   if ((tmp_v1 - tmp_v2) != tmp)
        //   {
        //     is_pass = my_alpha->is_pass(tmp_v1, tmp_v2, tmp);
        //   }
        // }
      }
    }

    return is_pass;
  };

  // idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha,grad, kernel, my_alpha, cek1);

  // if (idx_a == -1)
  // {
  //   idx_a = cari(idx_b, rho.rho_v1, rho.rho_v2, alpha,grad, kernel, my_alpha, cek1);
  // }

  auto cek = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (is_pass)
    {
      // is_pass = !(alpha[1][var_a.idx] <= alpha[1].lb()) or !(alpha[2][var_a.idx] <= alpha[2].lb());
      if (is_pass)
      {
        if (alpha[0].is_nol(var_b.idx))
        {
          is_pass = !alpha[0].is_nol(var_a.idx);
        }
      }
    }
    if (is_pass)
    {
      vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);

      Tmy_double delta_v1 = hsl_eta[0] * (grad[1][var_a.idx] - grad[1][var_b.idx]);
      Treturn_is_pass tmp_v1 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v1, alpha[1]);

      Tmy_double delta_v2 = hsl_eta[0] * (grad[2][var_a.idx] - grad[2][var_b.idx]);
      Treturn_is_pass tmp_v2 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v2, alpha[2]);

      Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
      // cout << "delta " << delta << endl;
      Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

      is_pass = (tmp_v1.is_pass or tmp_v2.is_pass);

      // if (is_pass)
      // {
      //   if ((tmp_v1 - tmp_v2) != tmp)
      //   {
      //     is_pass = my_alpha->is_pass(tmp_v1, tmp_v2, tmp);
      //   }
      // }
    }

    return is_pass;
  };

  if (idx_a == -1)
  {
    idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek);
  }

  if (idx_a == -1)
  {
    idx_a = cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek);
  }

  return idx_a;
}

bool Tmy_G::cari_idx(int &idx_b, int &idx_a, Treturn_update_rho rho, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  {
    bool is_pass = true;
    // is_pass = !(alpha[1][var_b.idx] >= alpha[1].ub()) or !(alpha[2][var_b.idx] >= alpha[2].ub());
    return is_pass;
  };

  auto cek_filter1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  {
    bool is_pass = true;
    // is_pass = !(alpha[1][var_a.idx] <= alpha[1].lb()) or !(alpha[2][var_a.idx] <= alpha[2].lb());
    if (is_pass)
    {
      if (alpha[0].is_nol(var_b.idx))
      {
        is_pass = !alpha[0].is_nol(var_a.idx);
      }
    }

    return is_pass;
  };

  idx_b = max(rho.rho_v1, rho.rho_v2, alpha, grad, cek_filter);
  if (idx_b != -1)
  {
    idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, cek_filter1);
  }

  return ((idx_b != -1) and (idx_a != -1));
}

Tmy_double Tmy_G::sum_alpha_diff_Q(T_alpha_container alpha, vector<Tmy_double> diff_Q)
{
  Tmy_double hasil = 0.0;

  for (int i = 0; i < diff_Q.size(); ++i)
  {
    hasil = hasil + (alpha[i] * diff_Q[i]);
  }

  return hasil;
}

Tmy_double Tmy_G::sum_alpha_rho_Q(int i, Tmy_kernel *kernel, T_alpha_container alpha)
{
  Tmy_double hasil = 0.0;
  vector<Tmy_double> hsl_Q = kernel->get_rho_Q(i);
  for (int i = 0; i < _jml_data; ++i)
  {
    hasil = hasil + (alpha[i] * hsl_Q[i]);
  }
  return hasil;
}

void Tmy_G::update_G(int idx_b, int idx_a, Treturn_is_pass tmp, Tmy_kernel *kernel, T_alpha_container &alpha, T_grad_container &grad)
{
  Tmy_double alpha_a = alpha[idx_a];
  Tmy_double alpha_b = alpha[idx_b];

  Tmy_double delta_1 = tmp.new_alpha_i - alpha_b;
  Tmy_double delta_2 = tmp.new_alpha_j - alpha_a;

  vector<Tmy_double> data_a = kernel->get_Q(idx_a);
  vector<Tmy_double> data_b = kernel->get_Q(idx_b);

  for (int i = 0; i < _jml_data; ++i)
  {
    grad[i] = grad[i] + ((data_b[i] * delta_1) + (data_a[i] * delta_2));
  }
}

int Tmy_G::max(Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, callback_type f)
{
  Tmy_double gmax = -HUGE_VAL;
  int idx_max = -1;

  vector<int> tmp_idx = grad[0].get_rand_idx();
  for (int i = 0; i < tmp_idx.size(); ++i)
  {
    Tmy_double dec_F = grad[0].dec(tmp_idx[i], rho1, rho2);
    Tmy_double obj_F = grad[0].obj(tmp_idx[i], rho1, rho2);

    callback_param var_b;
    var_b.idx = tmp_idx[i];
    var_b.dec = dec_F;
    var_b.obj = obj_F;
    var_b.grad = grad[0][tmp_idx[i]];

    callback_param var_a;
    bool is_pass = true;
    if (_cek_kkt)
    {
      is_pass = !grad[0].get_kkt(tmp_idx[i]);
    }
    //    cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      // is_pass = f(var_b, var_a, alpha);
    }

    if (is_pass)
    {
      Tmy_double abs_F = abs(grad[0][tmp_idx[i]]);
      if (_min_rho)
      {
        abs_F = abs(obj_F);
      }

      if (abs_F >= gmax)
      {
        gmax = abs_F;
        idx_max = tmp_idx[i];
        grad[0].mv_idx(tmp_idx[i], 0);
      }
    }
  }

  return idx_max;
}

int Tmy_G::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, callback_type f)
{
  Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL;
  int idx_max = -1;

  Tmy_double Gb = grad[0][idx_b];
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;

  Tmy_double Fa = 0.0;
  Tmy_double Fb = 0.0;

  vector<int> tmp_idx = grad[0].get_rand_idx();
  for (int i = 0; i < tmp_idx.size(); ++i)
  {
    Tmy_double Ga = grad[0][tmp_idx[i]];
    Tmy_double dec_Fa = grad[0].dec(tmp_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(tmp_idx[i], rho1, rho2);

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
      is_pass = !grad[0].get_kkt(tmp_idx[i]);
    }
    // cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad);
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

        vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, tmp_idx[i]);
        Tmy_double diff = Ga - Gb;
        Tmy_double delta = diff * hsl_eta[0];

        diff = grad[1][i] - grad[1][idx_b];
        Tmy_double delta_v1 = diff * hsl_eta[0];

        diff = grad[2][i] - grad[2][idx_b];
        Tmy_double delta_v2 = diff * hsl_eta[0];

        bool is_pass = true;
        if (_filter_delta)
        {

          is_pass = delta_filter(idx_b, tmp_idx[i], alpha, delta);
        }
        if (is_pass)
        {
          // cout << " idx_a " << var_a.idx << " delta " << delta << endl;
          gmax2 = abs_diff_F;
          idx_max = tmp_idx[i];
          grad[0].mv_idx(tmp_idx[i], 0);
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

int Tmy_G::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
  Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL;
  int idx_max = -1;

  Tmy_double Gb = grad[0][idx_b];
  Tmy_double abs_Gb = abs(Gb);
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;

  Tmy_double Fa = 0.0;
  Tmy_double Fb = 0.0;

  vector<int> tmp_idx = grad[0].get_rand_idx();
  for (int i = 0; i < tmp_idx.size(); ++i)
  {
    Tmy_double Ga = grad[0][tmp_idx[i]];
    Tmy_double abs_Ga = abs(Ga);
    Tmy_double dec_Fa = grad[0].dec(tmp_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(tmp_idx[i], rho1, rho2);

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
      is_pass = !grad[0].get_kkt(tmp_idx[i]);
    }
    //    cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad, kernel, my_alpha);
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
        vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, tmp_idx[i]);
        Tmy_double diff = Ga - Gb;
        Tmy_double delta = diff * hsl_eta[0];

        diff = grad[1][i] - grad[1][idx_b];
        Tmy_double delta_v1 = diff * hsl_eta[0];

        diff = grad[2][i] - grad[2][idx_b];
        Tmy_double delta_v2 = diff * hsl_eta[0];

        bool is_pass = true;
        if (_filter_delta)
        {
          is_pass = delta_filter(idx_b, tmp_idx[i], alpha, delta);
        }
        if (is_pass)
        {
          // cout << " idx_a " << var_a.idx << " delta " << delta << endl;
          gmax2 = abs_diff_F;
          idx_max = tmp_idx[i];
          grad[0].mv_idx(tmp_idx[i], 0);
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

int Tmy_G::cari(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
  // cout << " Cari 1 " << endl;
  int idx_a = -1;

  Tmy_double Gb = grad[0][idx_b];
  Tmy_double abs_Gb = abs(Gb);
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;

  vector<int> tmp_idx = grad[0].get_rand_idx();
  for (int i = 0; i < tmp_idx.size(); ++i)
  {
    Tmy_double Ga = grad[0][tmp_idx[i]];
    Tmy_double abs_Ga = abs(Ga);
    Tmy_double dec_Fa = grad[0].dec(tmp_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(tmp_idx[i], rho1, rho2);

    callback_param var_a;
    var_a.idx = tmp_idx[i];
    var_a.dec = dec_Fa;
    var_a.obj = obj_Fa;
    var_a.grad = Ga;

    bool is_pass = true;
    if (_cek_kkt)
    {
      is_pass = !grad[0].get_kkt(tmp_idx[i]);
    }
    //    cout << " idx a " << _idx[i] << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad, kernel, my_alpha);
    }

    if (is_pass)
    {
      vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, tmp_idx[i]);
      Tmy_double diff = Ga - Gb;
      Tmy_double delta = diff * hsl_eta[0];

      diff = grad[1][i] - grad[1][idx_b];
      Tmy_double delta_v1 = diff * hsl_eta[0];

      diff = grad[2][i] - grad[2][idx_b];
      Tmy_double delta_v2 = diff * hsl_eta[0];

      bool is_pass = true;
      if (_filter_delta)
      {
        is_pass = delta_filter(idx_b, tmp_idx[i], alpha, delta);
      }
      if (is_pass)
      {
        cout << " idx_a " << var_a.idx << " delta " << delta << endl;
        idx_a = tmp_idx[i];
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
}

bool Tmy_G::delta_filter(int idx_b, int idx_a, vector<T_alpha_container> alpha, Tmy_double delta)
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
