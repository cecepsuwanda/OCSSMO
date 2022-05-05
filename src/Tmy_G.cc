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

Treturn_update_rho Tmy_G::update_rho(Tmy_kernel *kernel, T_alpha_container alpha, T_grad_container grad)
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
    if (alpha.is_sv(i))
    {
      if (alpha[i] > 0.0) {
        jml_G_v1 = jml_G_v1 + sum_alpha_rho_Q(i, kernel, alpha);
        jml_n_v1++;
      }

      if (alpha[i] < 0.0) {
        jml_G_v2 = jml_G_v2 + sum_alpha_rho_Q(i, kernel, alpha);
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



bool Tmy_G::is_kkt(int idx, Treturn_update_rho rho,T_alpha_container alpha,T_grad_container grad)
{
  Tmy_double F = (grad[idx] - rho.rho_v1) * (rho.rho_v2 - grad[idx]);  
  bool stat = false;

  if ((alpha[idx]==0.0) and (F > 0.0))
  {
    stat = true;
  }
  else
  {
    if ((alpha.is_sv(idx)) and (abs(F) < 1e-3))
    {
      stat = true;
    }
    else
    {
      if (( alpha.is_lb(idx) or alpha.is_ub(idx)) and (F < 0.0))
      {
        stat = true;
      }
    }
  }

  return stat;
}

int Tmy_G::cari_idx_a(int idx_b, Treturn_update_rho rho, T_alpha_container alpha, T_grad_container grad)
{
  int idx_a = -1;
  Tmy_double gmax = -HUGE_VAL;

  Tmy_double Fb = grad.obj(idx_b, rho.rho_v1, rho.rho_v2);

  for (int j = 0; j < _jml_data; ++j)
  {
    Tmy_double Fa = grad.obj(j, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_diff = abs(Fb - Fa);

    bool is_pass = true;
    is_pass = alpha.is_sv(j);
    if ((abs_diff >= gmax) and is_pass)
    {
      idx_a = j;
      gmax = abs_diff;
      //cout << " " << min_Fb << " " << Fa << " ";
    }
  }

  return idx_a;
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho, Tmy_kernel *kernel, T_alpha_container alpha, T_grad_container grad, Tmy_alpha *my_alpha)
{

  auto cek = [&, this](int idx_b, int idx_a) -> bool
  {
    vector<Tmy_double> hsl_eta = kernel->hit_eta(idx_b, idx_a);
    vector<Tmy_double> hsl_diff = kernel->get_diff_Q(idx_b, idx_a);

    Tmy_double hsl_sum = this->sum_alpha_diff_Q(alpha, hsl_diff);
    Tmy_double delta = hsl_eta[0] * hsl_sum;

    Treturn_is_pass tmp = my_alpha->is_pass(idx_b, idx_a, delta, alpha);

    return (tmp.is_pass);
  };


  //srand(time(0));
  Tmy_double gmax2 = -HUGE_VAL;
  int idx_a = -1;

  Tmy_double Fb = grad.obj(idx_b, rho.rho_v1, rho.rho_v2);
  for (int j = 0; j < _jml_data; ++j)
  {
    Tmy_double Fa = grad.obj(j, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_diff = abs(Fb - Fa);
    bool is_pass = true;

    if (is_pass)
    {
      is_pass = alpha.is_sv(j);
    }
    if (is_pass)
    {
      is_pass = cek(idx_b, j);
    }
    if ((abs_diff >= gmax2) and is_pass)
    {
      idx_a = j;
      gmax2 = abs_diff;
    }

  }

  //int startIndex = (rand() % _jml_data);
  if (idx_a == -1)
  {
    for (int i = idx_b; i < _jml_data; ++i)
    {
      bool is_pass = true;

      if (is_pass)
      {
        is_pass = alpha.is_sv(i);
      }
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }

    }
  }

  if (idx_a == -1)
  {
    for (int i = 0; i < idx_b; ++i)
    {
      bool is_pass = true;
      if (is_pass)
      {
        is_pass = alpha.is_sv(i);
      }
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }

    }
  }

  if (idx_a == -1)
  {
    gmax2 = -HUGE_VAL;
    for (int j = 0; j < _jml_data; ++j)
    {
      Tmy_double Fa = grad.obj(j, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_diff = abs(Fb - Fa);
      bool is_pass = true;
      if (is_pass)
      {
        is_pass = cek(idx_b, j);
      }
      if ((abs_diff >= gmax2) and is_pass)
      {
        idx_a = j;
        gmax2 = abs_diff;
      }

    }

  }

  //startIndex = (rand() % _jml_data);
  if (idx_a == -1)
  {
    for (int i = idx_b; i < _jml_data; ++i)
    {
      bool is_pass = true;
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }

    }
  }

  if (idx_a == -1)
  {

    for (int i = 0; i < idx_b; ++i)
    {
      bool is_pass = true;
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }

    }
  }

  return idx_a;
}

bool Tmy_G::cari_idx(int& idx_b, int& idx_a, Treturn_update_rho rho)
{
  // Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_double gmax = -HUGE_VAL, gmin = HUGE_VAL, diff_f_max = -HUGE_VAL;

  idx_a = -1;
  idx_b = -1;
  // for (int i = 0; i < _jml_data; ++i)
  // {
  //   Tmy_double Gb = _my_list_G->get_G(i);
  //   Tmy_double Fb = _my_list_G->get_obj(i, rho.rho_v1, rho.rho_v2);
  //   Tmy_double abs_Fb = abs(Fb);
  //   //(abs((Gb - rho.rho_v1)*(rho.rho_v2 - Gb))>1e-3)
  //   bool is_pass = !is_kkt(i, rho);
  //   if (abs_Fb >= gmax)
  //   {
  //     if (is_pass) {
  //       idx_b = i;
  //     }
  //     gmax = abs_Fb;
  //   }

  // }

  // if (idx_b != -1)
  // {
  //   Tmy_double Fb = _my_list_G->get_obj(idx_b, rho.rho_v1, rho.rho_v2);
  //   for (int j = 0; j < _jml_data; ++j)
  //   {
  //     Tmy_double Fa = _my_list_G->get_obj(j, rho.rho_v1, rho.rho_v2);
  //     Tmy_double abs_Fa = abs(Fa);
  //     if (abs_Fa <= gmin)
  //     {
  //       gmin = abs_Fa;
  //     }

  //     Tmy_double diff_f = abs(Fb - Fa);
  //     bool is_pass = !list_alpha->is_nol(idx_b, j);
  //     if (is_pass)
  //     {
  //       is_pass = list_alpha->is_free(j);
  //     }
  //     if ((diff_f >= diff_f_max) and is_pass)
  //     {
  //       idx_a = j;
  //       diff_f_max = diff_f;
  //     }
  //   }
  // }

  Tmy_double diff = gmax - gmin;
  //cout << " " << gmax << "-" << gmin << "=" << diff << " ";
  return (abs(diff) > 1e-3) and ( (idx_b != -1) and (idx_a != -1) );
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

  alpha[idx_a] = tmp.new_alpha_j;
  alpha[idx_b] = tmp.new_alpha_i;
}



