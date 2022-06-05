#include "Tmy_G.h"

Tmy_G::Tmy_G()
{
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

  if ((alpha[1].is_nol(idx) and alpha[2].is_nol(idx)) and (F > 1e-3))
  {
    stat = true;
  }
  else
  {
    if (((alpha[1].is_sv(idx) and alpha[2].is_nol(idx)) or (alpha[1].is_nol(idx) and alpha[2].is_sv(idx))) and (abs(F) < 1e-3))
    {
      stat = true;
    }
    else
    {
      if (((alpha[1].is_ub(idx) and alpha[2].is_nol(idx)) or (alpha[1].is_nol(idx) and alpha[2].is_ub(idx))) and (F < -1e-3))
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
    grad.set_kkt(i, is_kkt(i, rho, alpha, grad));
  }
}

int Tmy_G::cari_idx_a(int idx_b, Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha) -> bool
  {
    bool is_pass = true;
    is_pass = alpha[0].is_sv(var_a.idx);
    return is_pass;
  };

  int idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, cek_filter);

  return idx_a;
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho, Tmy_kernel *kernel, vector<T_alpha_container> alpha, T_grad_container grad, Tmy_alpha *my_alpha)
{

  int idx_a = -1;

  auto cek1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (is_pass)
    {

      is_pass = alpha[0].is_sv(var_a.idx);

      if (is_pass)
      {
        Treturn_is_pass tmp_v1;
        tmp_v1.is_pass = true;
        tmp_v1.alpha_i = alpha[1][var_b.idx];
        tmp_v1.alpha_j = alpha[1][var_a.idx];
        tmp_v1.new_alpha_i = alpha[1][var_b.idx];
        tmp_v1.new_alpha_j = alpha[1][var_a.idx];
        tmp_v1.lb = alpha[1].lb();
        tmp_v1.ub = alpha[1].ub();

        Treturn_is_pass tmp_v2;
        tmp_v2.is_pass = true;
        tmp_v2.alpha_i = alpha[2][var_b.idx];
        tmp_v2.alpha_j = alpha[2][var_a.idx];
        tmp_v2.new_alpha_i = alpha[2][var_b.idx];
        tmp_v2.new_alpha_j = alpha[2][var_a.idx];
        tmp_v2.lb = alpha[2].lb();
        tmp_v2.ub = alpha[2].ub();

        vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);
        Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
        Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

        is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass)

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

  idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, my_alpha, cek1);

  if (idx_a == -1)
  {
    idx_a = grad.cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, my_alpha, cek1);
  }

  auto cek = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha, Tmy_kernel *kernel, Tmy_alpha *my_alpha) -> bool
  {
    bool is_pass = true;

    if (is_pass)
    {
      is_pass = !(alpha[0].is_nol(var_a.idx));
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

      Treturn_is_pass tmp_v1;
      tmp_v1.is_pass = true;
      tmp_v1.alpha_i = alpha[1][var_b.idx];
      tmp_v1.alpha_j = alpha[1][var_a.idx];
      tmp_v1.new_alpha_i = alpha[1][var_b.idx];
      tmp_v1.new_alpha_j = alpha[1][var_a.idx];
      tmp_v1.lb = alpha[1].lb();
      tmp_v1.ub = alpha[1].ub();

      Treturn_is_pass tmp_v2;
      tmp_v2.is_pass = true;
      tmp_v2.alpha_i = alpha[2][var_b.idx];
      tmp_v2.alpha_j = alpha[2][var_a.idx];
      tmp_v2.new_alpha_i = alpha[2][var_b.idx];
      tmp_v2.new_alpha_j = alpha[2][var_a.idx];
      tmp_v2.lb = alpha[2].lb();
      tmp_v2.ub = alpha[2].ub();

      vector<Tmy_double> hsl_eta = kernel->hit_eta(var_b.idx, var_a.idx);
      Tmy_double delta = hsl_eta[0] * (var_a.grad - var_b.grad);
      // cout << "delta " << delta << endl;
      Treturn_is_pass tmp = my_alpha->is_pass(var_b.idx, var_a.idx, delta, alpha[0]);

      is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass)

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
    idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, my_alpha, cek);
  }

  if (idx_a == -1)
  {
    idx_a = grad.cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, my_alpha, cek);
  }

  return idx_a;
}

bool Tmy_G::cari_idx(int &idx_b, int &idx_a, Treturn_update_rho rho, vector<T_alpha_container> alpha, T_grad_container grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha) -> bool
  {
    bool is_pass = true;
    is_pass = ((var_b.dec < -1e-3) and (!alpha[0].is_ub(var_b.idx) or !alpha[0].is_lb(var_b.idx))) or ((var_b.dec > 1e-3) and (!alpha[0].is_nol(var_b.idx)));
    // if (is_pass)
    // {
    //   is_pass = !alpha[0].is_ub(var_b.idx) or !alpha[0].is_lb(var_b.idx);
    // }
    return is_pass;
  };

  auto cek_filter1 = [](callback_param var_b, callback_param var_a, vector<T_alpha_container> alpha) -> bool
  {
    bool is_pass = true;
    // is_pass = !alpha[0].is_nol(var_a.idx);
    return is_pass;
  };

  idx_b = grad.max(rho.rho_v1, rho.rho_v2, alpha, cek_filter);
  if (idx_b != -1)
  {
    idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, kernel, cek_filter1);
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
