#include "Tmy_G.h"

Tmy_G::Tmy_G()
{
  _cek_kkt = false;
  _filter_delta = false;
  _min_rho = true;
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

  int jml_n_v2 = 0;
  Tmy_double jml_G_v2 = 0.0;

  for (int i = 0; i < _jml_data; ++i)
  {
    if (alpha[2].is_nol(i))
    {
      if (alpha[1].is_sv(i))
      {
        jml_G_v1 = jml_G_v1 + grad[i];
        jml_n_v1++;
      }
      else
      {
        if (alpha[1].is_sv(i, alpha[2].lb(), alpha[2].ub()))
        {
          jml_G_v1 = jml_G_v1 + grad[i];
          jml_n_v1++;
        }
      }
    }
    if (alpha[1].is_nol(i))
    {
      if (alpha[2].is_sv(i))
      {
        jml_G_v2 = jml_G_v2 + grad[i];
        jml_n_v2++;
      }
      else
      {
        if (alpha[2].is_sv(i, alpha[1].lb(), alpha[1].ub()))
        {
          jml_G_v2 = jml_G_v2 + grad[i];
          jml_n_v2++;
        }
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

Tmy_double Tmy_G::update_rho(Tmy_kernel *kernel, T_alpha_container alpha, T_grad_container grad)
{
  Tmy_double jml_G = 0.0;
  Tmy_double G_max = -HUGE_VAL;
  Tmy_double G_min = HUGE_VAL;

  int jml_sv = 0;
  for (int i = 0; i < _jml_data; ++i)
  {
    if (alpha.is_sv(i))
    {
      jml_G = jml_G + grad[i];
      jml_sv = jml_sv + 1;
    }
    else
    {
      if (alpha.is_lb(i))
      {
        if (grad[i] < G_min)
        {
          G_min = grad[i];
        }
      }
      else
      {
        if (alpha.is_ub(i))
        {
          if (grad[i] > G_max)
          {
            G_max = grad[i];
          }
        }
      }
    }
  }

  Tmy_double tmp_rho = 0.0;
  if (jml_sv > 0)
  {
    tmp_rho = jml_G / (1.0 * jml_sv);
  }
  else
  {
    tmp_rho = (G_min + G_max) / 2.0;
  }

  return tmp_rho;
}

bool Tmy_G::is_kkt(int idx, Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container grad)
{
  Tmy_double F = grad.dec(idx, rho.rho_v1, rho.rho_v2);

  // bool stat1 = (grad[idx] > rho.rho_v1) and (grad[idx] < rho.rho_v2) and alpha[1].is_nol(idx) and alpha[2].is_nol(idx);
  // bool stat2 = (grad[idx] == rho.rho_v1) and (grad[idx] < rho.rho_v2) and alpha[1].is_sv(idx) and alpha[2].is_nol(idx);
  // bool stat3 = (grad[idx] == rho.rho_v2) and alpha[1].is_nol(idx) and alpha[2].is_sv(idx);
  // bool stat4 = (grad[idx] < rho.rho_v1) and (grad[idx] < rho.rho_v2) and alpha[1].is_ub(idx) and alpha[2].is_nol(idx);
  // bool stat5 = (grad[idx] > rho.rho_v1) and (grad[idx] > rho.rho_v2) and alpha[1].is_nol(idx) and alpha[2].is_ub(idx);

  // return (stat1 or stat2 or stat3 or stat4 or stat5);

  bool stat1 = alpha[0].is_nol(idx) and (F > 1e-3);
  bool stat2 = alpha[0].is_sv(idx) and (alpha[1].is_nol(idx) or alpha[2].is_nol(idx)) and (abs(F) <= 1e-3);
  bool stat3 = (alpha[0].is_ub(idx) or alpha[0].is_lb(idx)) and (alpha[1].is_nol(idx) or alpha[2].is_nol(idx)) and (F < -1e-3);
  return (stat1 or stat2 or stat3);
}

bool Tmy_G::is_kkt(int idx, Tmy_double rho, T_alpha_container alpha, T_grad_container grad)
{
  Tmy_double dec = grad.dec(idx, rho);
  bool stat1 = alpha.is_nol(idx) and (dec > 1e-3);
  bool stat2 = alpha.is_sv(idx) and (abs(dec) <= 1e-3);
  bool stat3 = alpha.is_ub(idx) and (dec < -1e-3);
  return (stat1 or stat2 or stat3);
}

void Tmy_G::set_kkt(Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container &grad)
{
  for (int i = 0; i < _jml_data; ++i)
  {
    bool tmp = is_kkt(i, rho, alpha, grad);
    grad.set_kkt(i, tmp);
    if (!tmp)
    {
      grad.mv_idx(i, 0);
    }
    else
    {
      grad.mv_idx(i, 1);
    }
  }
}

void Tmy_G::set_kkt(Tmy_double rho, T_alpha_container alpha, T_grad_container &grad)
{
  for (int i = 0; i < _jml_data; ++i)
  {
    bool tmp = is_kkt(i, rho, alpha, grad);
    grad.set_kkt(i, tmp);

    if (!tmp)
    {
      grad.mv_idx(i, 0);
    }
    else
    {
      grad.mv_idx(i, 1);
    }
  }
}

int Tmy_G::cari_idx_a(int idx_b, Treturn_update_rho rho, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  {
    bool is_pass = true;

    // bool kondisi1 = (((var_b.grad - var_b.rho_v1) < -1e-3) and (alpha[1][var_b.idx] < alpha[1].ub())) or (((var_b.grad - var_b.rho_v1) > 1e-3) and (alpha[1][var_b.idx] > alpha[1].lb()));
    // bool kondisi2 = (((var_b.rho_v2 - var_b.grad) < -1e-3) and (alpha[2][var_b.idx] < alpha[2].ub())) or (((var_b.rho_v2 - var_b.grad) > 1e-3) and (alpha[2][var_b.idx] > alpha[2].lb()));

    // is_pass = kondisi1 and kondisi2;

    // bool stat1 = (grad[0][var_b.idx] > var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat2 = (grad[0][var_b.idx] == var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_sv(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat3 = (grad[0][var_b.idx] == var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_sv(var_b.idx);
    // bool stat4 = (grad[0][var_b.idx] < var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_ub(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat5 = (grad[0][var_b.idx] > var_b.rho_v1) and (grad[0][var_b.idx] > var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_ub(var_b.idx);

    // is_pass = !(stat1 or stat2 or stat3 or stat4 or stat5);
    Tmy_double F = var_b.dec;
    bool stat1 = alpha[0].is_nol(var_b.idx) and (F > 1e-3);
    bool stat2 = alpha[0].is_sv(var_b.idx) and (alpha[1].is_nol(var_b.idx) or alpha[2].is_nol(var_b.idx)) and (abs(F) <= 1e-3);
    bool stat3 = (alpha[0].is_ub(var_b.idx) or alpha[0].is_lb(var_b.idx)) and (alpha[1].is_nol(var_b.idx) or alpha[2].is_nol(var_b.idx)) and (F < -1e-3);
    is_pass = !(stat1 or stat2 or stat3);

    if (is_pass)
    {
      is_pass = !alpha[0].is_nol(var_b.idx); //! alpha[1].is_nol(var_b.idx) or !alpha[2].is_nol(var_b.idx);
    }

    if (is_pass)
    {
      if (alpha[0].is_nol(var_b.idx))
      {
        is_pass = !alpha[0].is_nol(var_a.idx);
      }
    }

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
      is_pass = alpha[0].is_sv(var_a.idx) and (alpha[1].is_nol(var_a.idx) or alpha[2].is_nol(var_a.idx));

      if (is_pass)
      {

        vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);

        Tmy_double delta_v1 = hsl_eta[0] * (grad[1][var_a.idx] - grad[1][var_b.idx]);
        Treturn_is_pass tmp_v1 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v1, alpha[1]);
        tmp_v1.reset();

        Tmy_double delta_v2 = hsl_eta[0] * (grad[2][var_a.idx] - grad[2][var_b.idx]);
        Treturn_is_pass tmp_v2 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v2, alpha[2]);
        tmp_v2.reset();

        Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
        //  cout << "delta " << delta << endl;

        Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

        is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass);

        if (is_pass)
        {
          if ((tmp_v1 - tmp_v2) != tmp)
          {
            is_pass = my_alpha->is_pass(tmp_v1, tmp_v2, tmp);
          }
        }
      }
    }

    return is_pass;
  };

  // cout << " max " << endl;
  idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek1);
  if (idx_a == -1)
  {
    // cout << " cari 1 " << endl;
    idx_a = cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek1);
  }

  auto cek = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (alpha[0].is_nol(var_b.idx))
    {
      is_pass = !alpha[0].is_nol(var_a.idx);
    }

    if (is_pass)
    {
      vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);

      Tmy_double delta_v1 = hsl_eta[0] * (grad[1][var_a.idx] - grad[1][var_b.idx]);
      Treturn_is_pass tmp_v1 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v1, alpha[1]);
      tmp_v1.reset();

      Tmy_double delta_v2 = hsl_eta[0] * (grad[2][var_a.idx] - grad[2][var_b.idx]);
      Treturn_is_pass tmp_v2 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v2, alpha[2]);
      tmp_v2.reset();

      Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
      // cout << "delta " << delta << endl;

      Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

      is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass);

      if (is_pass)
      {
        if ((tmp_v1 - tmp_v2) != tmp)
        {
          is_pass = my_alpha->is_pass(tmp_v1, tmp_v2, tmp);
        }
      }
    }

    return is_pass;
  };

  if (idx_a == -1)
  {
    // cout << " cari 2 " << endl;
    idx_a = cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek);
  }

  return idx_a;
}

bool Tmy_G::cari_idx(int &idx_b, int &idx_a, Treturn_update_rho rho, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  {
    bool is_pass = true;

    // bool kondisi1 = (((var_b.grad - var_b.rho_v1) < -1e-3) and (alpha[1][var_b.idx] < alpha[1].ub())) or (((var_b.grad - var_b.rho_v1) > 1e-3) and (alpha[1][var_b.idx] > alpha[1].lb()));
    // bool kondisi2 = (((var_b.rho_v2 - var_b.grad) < -1e-3) and (alpha[2][var_b.idx] < alpha[2].ub())) or (((var_b.rho_v2 - var_b.grad) > 1e-3) and (alpha[2][var_b.idx] > alpha[2].lb()));
    // is_pass = kondisi1 and kondisi2;

    // bool stat1 = (grad[0][var_b.idx] > var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat2 = (grad[0][var_b.idx] == var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_sv(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat3 = (grad[0][var_b.idx] == var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_sv(var_b.idx);
    // bool stat4 = (grad[0][var_b.idx] < var_b.rho_v1) and (grad[0][var_b.idx] < var_b.rho_v2) and alpha[1].is_ub(var_b.idx) and alpha[2].is_nol(var_b.idx);
    // bool stat5 = (grad[0][var_b.idx] > var_b.rho_v1) and (grad[0][var_b.idx] > var_b.rho_v2) and alpha[1].is_nol(var_b.idx) and alpha[2].is_ub(var_b.idx);

    // is_pass = !(stat1 or stat2 or stat3 or stat4 or stat5);
    Tmy_double F = var_b.dec;
    bool stat1 = alpha[0].is_nol(var_b.idx) and (F > 1e-3);
    bool stat2 = alpha[0].is_sv(var_b.idx) and (alpha[1].is_nol(var_b.idx) or alpha[2].is_nol(var_b.idx)) and (abs(F) <= 1e-3);
    bool stat3 = (alpha[0].is_ub(var_b.idx) or alpha[0].is_lb(var_b.idx)) and (alpha[1].is_nol(var_b.idx) or alpha[2].is_nol(var_b.idx)) and (F < -1e-3);
    is_pass = !(stat1 or stat2 or stat3);

    return is_pass;
  };

  // auto cek_filter1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad) -> bool
  // {
  //   bool is_pass = true;

  //   if (alpha[0].is_nol(var_b.idx))
  //   {
  //     is_pass = !alpha[0].is_nol(var_a.idx);
  //   }

  //   return is_pass;
  // };

  auto cek_filter1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (alpha[0].is_nol(var_b.idx))
    {
      is_pass = !alpha[0].is_nol(var_a.idx);
    }

    if (is_pass)
    {
      vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);

      Tmy_double delta_v1 = hsl_eta[0] * (grad[1][var_a.idx] - grad[1][var_b.idx]);
      Treturn_is_pass tmp_v1 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v1, alpha[1]);
      tmp_v1.reset();

      Tmy_double delta_v2 = hsl_eta[0] * (grad[2][var_a.idx] - grad[2][var_b.idx]);
      Treturn_is_pass tmp_v2 = my_alpha->is_pass(var_b.idx, var_a.idx, delta_v2, alpha[2]);
      tmp_v2.reset();

      Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
      // cout << "delta " << delta << endl;

      Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

      is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass);

      if (is_pass)
      {
        if ((tmp_v1 - tmp_v2) != tmp)
        {
          is_pass = my_alpha->is_pass(tmp_v1, tmp_v2, tmp);
        }
      }
    }

    return is_pass;
  };

  idx_b = max(rho.rho_v1, rho.rho_v2, alpha, grad, cek_filter);
  if (idx_b != -1)
  {
    idx_a = max(idx_b, rho.rho_v1, rho.rho_v2, alpha, grad, kernel, my_alpha, cek_filter1);
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

  Tmy_double Fb = 0.0;

  vector<int> rand_idx = grad[0].get_rand_idx();
  for (int i = 0; i < _jml_data; ++i)
  {
    Tmy_double G = grad[0][rand_idx[i]];
    Tmy_double dec_F = grad[0].dec(rand_idx[i], rho1, rho2);
    Tmy_double obj_F = grad[0].obj(rand_idx[i], rho1, rho2);

    callback_param var_b;
    var_b.idx = rand_idx[i];
    var_b.dec = dec_F;
    var_b.obj = obj_F;
    var_b.grad = G;
    var_b.rho_v1 = rho1;
    var_b.rho_v2 = rho2;

    Fb = G;
    if (_min_rho)
    {
      Fb = obj_F;
    }

    callback_param var_a;
    bool is_pass = true;
    if (_cek_kkt)
    {
      is_pass = !grad[0].get_kkt(rand_idx[i]);
    }
    //    cout << " idx a " << tmp_idx[i] << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad);
    }

    if (is_pass)
    {
      Tmy_double abs_F = abs(obj_F);
      if (abs_F >= gmax)
      {
        gmax = abs_F;
        idx_max = rand_idx[i];
        grad[0].mv_idx(rand_idx[i], 0);
      }
    }
  }

  return idx_max;
}

int Tmy_G::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, callback_type f)
{
  Tmy_double gmax = -HUGE_VAL;
  int idx_max = -1;

  Tmy_double Gb_v1 = grad[1][idx_b];
  Tmy_double Gb_v2 = grad[2][idx_b];
  Tmy_double Gb = grad[0][idx_b];
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;
  var_b.rho_v1 = rho1;
  var_b.rho_v2 = rho2;

  Tmy_double Fa = 0.0;
  Tmy_double Fb = 0.0;

  vector<int> rand_idx = grad[0].get_rand_idx();
  for (int i = 0; i < _jml_data; ++i)
  {
    Tmy_double Ga_v1 = grad[1][rand_idx[i]];
    Tmy_double Ga_v2 = grad[2][rand_idx[i]];
    Tmy_double Ga = grad[0][rand_idx[i]];
    Tmy_double dec_Fa = grad[0].dec(rand_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(rand_idx[i], rho1, rho2);

    callback_param var_a;
    var_a.idx = rand_idx[i];
    var_a.dec = dec_Fa;
    var_a.obj = obj_Fa;
    var_a.grad = Ga;
    var_a.rho_v1 = rho1;
    var_a.rho_v2 = rho2;

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
      is_pass = !grad[0].get_kkt(var_a.idx);
    }
    // cout << " idx a " << var_a.idx << endl; // << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad);
    }

    if (is_pass)
    {

      Tmy_double diff_F = Fb - Fa;
      Tmy_double abs_diff_F = abs(diff_F);
      // cout << " idx_a " << i << " abs_diff_obj " << abs_diff_obj << " gmax2 " << gmax2 << endl;
      if ((abs_diff_F >= gmax))
      {

        vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, var_a.idx);
        Tmy_double diff = Ga - Gb;
        Tmy_double delta = diff * hsl_eta[0];

        diff = Ga_v1 - Gb_v1;
        Tmy_double delta_v1 = diff * hsl_eta[0];

        diff = Ga_v2 - Gb_v2;
        Tmy_double delta_v2 = diff * hsl_eta[0];

        bool is_pass = true;
        if (_filter_delta)
        {
          vector<Tmy_double> tmp_delta;
          tmp_delta.push_back(delta);
          tmp_delta.push_back(delta_v1);
          tmp_delta.push_back(delta_v2);
          is_pass = delta_filter(idx_b, var_a.idx, alpha, tmp_delta);
        }
        if (is_pass)
        {
          // cout << " idx_a " << var_a.idx << " delta " << delta << " delta v1 " << delta_v1 << " delta v2 " << delta_v2 << endl;
          gmax = abs_diff_F;
          idx_max = var_a.idx;
          grad[0].mv_idx(var_a.idx, 0);
        }
      }
    }
  }
  return idx_max;
}

int Tmy_G::max(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
  Tmy_double gmax = -HUGE_VAL;
  int idx_max = -1;

  Tmy_double Gb_v1 = grad[1][idx_b];
  Tmy_double Gb_v2 = grad[2][idx_b];
  Tmy_double Gb = grad[0][idx_b];
  Tmy_double abs_Gb = abs(Gb);
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;
  var_b.rho_v1 = rho1;
  var_b.rho_v2 = rho2;

  Tmy_double Fa = 0.0;
  Tmy_double Fb = 0.0;

  vector<int> rand_idx = grad[0].get_rand_idx();
  for (int i = 0; i < _jml_data; ++i)
  {
    Tmy_double Ga_v1 = grad[1][rand_idx[i]];
    Tmy_double Ga_v2 = grad[2][rand_idx[i]];
    Tmy_double Ga = grad[0][rand_idx[i]];
    Tmy_double abs_Ga = abs(Ga);
    Tmy_double dec_Fa = grad[0].dec(rand_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(rand_idx[i], rho1, rho2);

    callback_param var_a;
    var_a.idx = rand_idx[i];
    var_a.dec = dec_Fa;
    var_a.obj = obj_Fa;
    var_a.grad = Ga;
    var_a.rho_v1 = rho1;
    var_a.rho_v2 = rho2;

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
      is_pass = !grad[0].get_kkt(rand_idx[i]);
    }
    //    cout << " idx a " << i << " is_pass " << is_pass << endl;
    if (is_pass)
    {
      is_pass = f(var_b, var_a, alpha, grad, kernel, my_alpha);
    }

    if (is_pass)
    {

      Tmy_double diff_F = Fb - Fa;
      Tmy_double abs_diff_F = abs(diff_F);
      // cout << " idx_a " << i << " abs_diff_obj " << abs_diff_obj << " gmax2 " << gmax2 << endl;
      if ((abs_diff_F >= gmax))
      {
        vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, rand_idx[i]);
        Tmy_double diff = Ga - Gb;
        Tmy_double delta = diff * hsl_eta[0];

        diff = Ga_v1 - Gb_v1;
        Tmy_double delta_v1 = diff * hsl_eta[0];

        diff = Ga_v2 - Gb_v2;
        Tmy_double delta_v2 = diff * hsl_eta[0];

        bool is_pass = true;
        if (_filter_delta)
        {
          vector<Tmy_double> tmp_delta;
          tmp_delta.push_back(delta);
          tmp_delta.push_back(delta_v1);
          tmp_delta.push_back(delta_v2);
          is_pass = delta_filter(idx_b, rand_idx[i], alpha, tmp_delta);
        }
        if (is_pass)
        {
          // cout << " idx_a " << var_a.idx << " delta " << delta << " delta v1 " << delta_v1 << " delta v2 " << delta_v2 << endl;
          gmax = abs_diff_F;
          idx_max = rand_idx[i];
          grad[0].mv_idx(rand_idx[i], 0);
        }
      }
    }
  }
  return idx_max;
}

int Tmy_G::cari(int idx_b, Tmy_double rho1, Tmy_double rho2, vector<T_alpha_container> alpha, vector<T_grad_container> grad, Tmy_kernel *kernel, Tmy_alpha *my_alpha, callback_type1 f)
{
  // cout << " Cari 1 " << endl;
  int idx_a = -1;

  Tmy_double Gb_v1 = grad[1][idx_b];
  Tmy_double Gb_v2 = grad[2][idx_b];
  Tmy_double Gb = grad[0][idx_b];
  Tmy_double abs_Gb = abs(Gb);
  Tmy_double dec_Fb = grad[0].dec(idx_b, rho1, rho2);
  Tmy_double obj_Fb = grad[0].obj(idx_b, rho1, rho2);

  callback_param var_b;
  var_b.idx = idx_b;
  var_b.dec = dec_Fb;
  var_b.obj = obj_Fb;
  var_b.grad = Gb;
  var_b.rho_v1 = rho1;
  var_b.rho_v2 = rho2;

  Tmy_double Fa = 0.0;
  Tmy_double Fb = 0.0;

  vector<int> tmp_idx = grad[0].get_rand_idx();
  for (int i = 0; i < tmp_idx.size(); ++i)
  {
    Tmy_double Ga_v1 = grad[1][tmp_idx[i]];
    Tmy_double Ga_v2 = grad[2][tmp_idx[i]];
    Tmy_double Ga = grad[0][tmp_idx[i]];
    Tmy_double abs_Ga = abs(Ga);
    Tmy_double dec_Fa = grad[0].dec(tmp_idx[i], rho1, rho2);
    Tmy_double obj_Fa = grad[0].obj(tmp_idx[i], rho1, rho2);

    callback_param var_a;
    var_a.idx = tmp_idx[i];
    var_a.dec = dec_Fa;
    var_a.obj = obj_Fa;
    var_a.grad = Ga;
    var_a.rho_v1 = rho1;
    var_a.rho_v2 = rho2;

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

      diff = Ga_v1 - Gb_v1;
      Tmy_double delta_v1 = diff * hsl_eta[0];

      diff = Ga_v2 - Gb_v2;
      Tmy_double delta_v2 = diff * hsl_eta[0];

      bool is_pass = true;
      if (_filter_delta)
      {
        vector<Tmy_double> tmp_delta;
        tmp_delta.push_back(delta);
        tmp_delta.push_back(delta_v1);
        tmp_delta.push_back(delta_v2);
        is_pass = delta_filter(idx_b, tmp_idx[i], alpha, tmp_delta);
      }
      if (is_pass)
      {
        // cout << " idx_a " << var_a.idx << " delta " << delta << " delta v1 " << delta_v1 << " delta v2 " << delta_v2 << endl;
        idx_a = tmp_idx[i];
        break;
      }
    }
  }

  return idx_a;
}

bool Tmy_G::delta_filter(int idx_b, int idx_a, vector<T_alpha_container> alpha, vector<Tmy_double> delta)
{
  bool is_pass = true;

  return is_pass;
}
