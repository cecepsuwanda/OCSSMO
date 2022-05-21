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

Treturn_update_rho Tmy_G::update_rho(Tmy_kernel *kernel, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2)
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
    if (alpha_v1.is_sv(i))
    {
      jml_G_v1 = jml_G_v1 + sum_alpha_rho_Q(i, kernel, alpha);
      jml_n_v1++;
    }

    if (alpha_v2.is_sv(i))
    {
      jml_G_v2 = jml_G_v2 + sum_alpha_rho_Q(i, kernel, alpha);
      jml_n_v2++;
    }
  }

  // cout << "[" << jml_n_v1 << "," << jml_n_v2 << "," << jml_G_v1 << "," << jml_G_v2 << "]" << endl;

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

bool Tmy_G::is_kkt(int idx, Treturn_update_rho rho, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, T_grad_container grad)
{
  Tmy_double F = grad.dec(idx, rho.rho_v1, rho.rho_v2);
  bool stat = false;

  if (((abs(alpha_v1[idx]) < 1e-3) and (abs(alpha_v2[idx]) < 1e-3)) and (F >= 1e-3))
  {
    stat = true;
  }
  else
  {
    if (( (alpha_v1.is_sv(idx) and (abs(alpha_v2[idx]) < 1e-3) ) or ( (abs(alpha_v1[idx]) < 1e-3) and alpha_v2.is_sv(idx)) ) and (abs(F) < 1e-3))
    {
      stat = true;
    }
    else
    {
      if (( (alpha_v1.is_ub(idx) and (abs(alpha_v2[idx]) < 1e-3))  or ((abs(alpha_v1[idx]) < 1e-3) and alpha_v2.is_ub(idx)) ) and (F <= -1e-3))
      {
        stat = true;
      }
    }
  }

  return stat;
}

int Tmy_G::cari_idx_a(int idx_b, Treturn_update_rho rho, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, T_grad_container grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](int idx_b, int idx_a, Tmy_double Fa, Tmy_double Fb, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2) -> bool
  {
    bool is_pass = true;

    if (alpha.is_nol(idx_b))
    {
      is_pass = !alpha.is_nol(idx_a);
    }

    if (is_pass)
    {
      is_pass = alpha.is_sv(idx_a);
    }

    // if (is_pass)
    // {
    //   is_pass = abs(Fb - Fa) > 1e-3;
    // }

    return is_pass;
  };

  int idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, cek_filter);

  return idx_a;
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho, Tmy_kernel *kernel, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, T_grad_container grad, Tmy_alpha *my_alpha)
{

  int idx_a = -1;

  auto cek1 = [](int idx_b, int idx_a, Tmy_double Fa, Tmy_double Fb, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, Tmy_kernel * kernel, Tmy_alpha * my_alpha) -> bool
  {
    bool is_pass = false;

    vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, idx_a);
    vector<Tmy_double> hsl_diff = kernel->get_diff_Q(idx_b, idx_a);

    Tmy_double hsl_sum = 0.0;
    Tmy_double delta = 0.0;

    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha_v1[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp_v1 = my_alpha->is_pass(idx_b, idx_a, delta, alpha_v1);

    hsl_sum = 0.0;
    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha_v2[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp_v2 = my_alpha->is_pass(idx_b, idx_a, delta, alpha_v2);


    hsl_sum = 0.0;
    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp = my_alpha->is_pass(idx_b, idx_a, delta, alpha);

    is_pass = (tmp_v1.is_pass == true) or (tmp_v2.is_pass == true);

    if (is_pass) {
      is_pass = alpha_v1.is_sv(idx_a) or alpha_v2.is_sv(idx_a);
    }

    // if (is_pass)
    // {
    //   if (alpha.is_nol(idx_b))
    //   {
    //     tmp.is_pass = !alpha.is_nol(idx_a);
    //   }
    // }

    // if (tmp.is_pass)
    // {
    //   tmp.is_pass = (abs(Fb - Fa) > 1e-3);
    // }

    return is_pass;
  };

  // idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, my_alpha, cek1);

  // if (idx_a == -1)
  // {
  //   idx_a = grad.cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, my_alpha, cek1);
  // }

  auto cek = [](int idx_b, int idx_a, Tmy_double Fa, Tmy_double Fb, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, Tmy_kernel * kernel, Tmy_alpha * my_alpha) -> bool
  {
    bool is_pass = false;

    vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, idx_a);
    vector<Tmy_double> hsl_diff = kernel->get_diff_Q(idx_b, idx_a);

    Tmy_double hsl_sum = 0.0;
    Tmy_double delta = 0.0;

    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha_v1[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp_v1 = my_alpha->is_pass(idx_b, idx_a, delta, alpha_v1);

    hsl_sum = 0.0;
    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha_v2[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp_v2 = my_alpha->is_pass(idx_b, idx_a, delta, alpha_v2);


    hsl_sum = 0.0;
    for (int i = 0; i < hsl_diff.size(); ++i)
    {
      hsl_sum = hsl_sum + (alpha[i] * hsl_diff[i]);
    }
    delta = hsl_eta[0] * hsl_sum;
    Treturn_is_pass tmp = my_alpha->is_pass(idx_b, idx_a, delta, alpha);

    is_pass = (tmp_v1.is_pass == true) or (tmp_v2.is_pass == true);

    if (is_pass) {
      is_pass = !(alpha_v1[idx_a] <= alpha_v1.lb()) or !(alpha_v2[idx_a] >= alpha_v2.ub());
    }

    // if (is_pass)
    // {
    //   if (alpha.is_nol(idx_b))
    //   {
    //     tmp.is_pass = !alpha.is_nol(idx_a);
    //   }
    // }

    // if (tmp.is_pass)
    // {
    //   tmp.is_pass = (abs(Fb - Fa) > 1e-3);
    // }

    return is_pass;
  };

  // if (idx_a == -1)
  // {
  //   idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, my_alpha, cek);
  // }

  if (idx_a == -1)
  {
    idx_a = grad.cari(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, my_alpha, cek);
  }

  return idx_a;
}

bool Tmy_G::cari_idx(int& idx_b, int& idx_a, Treturn_update_rho rho, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2, T_grad_container grad, Tmy_kernel *kernel)
{
  auto cek_filter = [](int idx_b, int idx_a, Tmy_double Fa, Tmy_double Fb, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2) -> bool
  {
    bool is_pass = true;
    is_pass = !(alpha_v1[idx_a] >= alpha_v1.ub()) or !(alpha_v2[idx_a] <= alpha_v2.lb());
    return is_pass;
  };

  auto cek_filter1 = [](int idx_b, int idx_a, Tmy_double Fa, Tmy_double Fb, T_alpha_container alpha, T_alpha_container alpha_v1, T_alpha_container alpha_v2) -> bool
  {
    bool is_pass = true;
    is_pass = !(alpha_v1[idx_a] <= alpha_v1.lb()) and !(alpha_v2[idx_a] >= alpha_v2.ub());
    //if (is_pass) {
    //  is_pass = alpha_v1.is_sv(idx_a) or alpha_v2.is_sv(idx_a);
    //}
    return is_pass;
  };

  idx_b = grad.max(rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, cek_filter);
  if (idx_b != -1)
  {
    idx_a = grad.max(idx_b, rho.rho_v1, rho.rho_v2, alpha, alpha_v1, alpha_v2, kernel, cek_filter1);
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
