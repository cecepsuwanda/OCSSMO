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

Treturn_cari_idx Tmy_G::cari_idx(Treturn_update_rho rho)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();

  int idx_b = -1, idx_a = -1;
  Tmy_double gmax = -HUGE_VAL, gmax2 = -HUGE_VAL, diff_max = -HUGE_VAL;


  for (int i = 0; i < _jml_data; ++i)
  {
    Tmy_double Gb = _my_list_G->get_G(i);
    Tmy_double Fb = _my_list_G->get_F(i, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_Fb = abs(Fb);
    if (abs(Gb) > 1e-3)
    {
      if (abs_Fb >= gmax)
      {
        idx_b = i;
        gmax = abs_Fb;
      }
    }
  }

  if (idx_b != -1)
  {
    Tmy_double Fb = _my_list_G->get_F(idx_b, rho.rho_v1, rho.rho_v2);
    for (int j = 0; j < _jml_data; ++j)
    {
      Tmy_double Fa = _my_list_G->get_F(j, rho.rho_v1, rho.rho_v2);
      Tmy_double abs_diff = abs(Fb - Fa);
      bool is_pass = !list_alpha->is_nol(idx_b, j);
      if ((abs_diff >= gmax2) and is_pass)
      {
        idx_a = j;
        gmax2 = abs_diff;
      }
    }
  }

  return {idx_b, idx_a, gmax, gmax2};
}

int Tmy_G::cari_idx_lain(int idx_b, Treturn_update_rho rho)
{
  Tmy_list_alpha *list_alpha = _alphas->get_alpha();
  Tmy_list_alpha *list_alpha_v1 = _alphas->get_alpha_v1();
  Tmy_list_alpha *list_alpha_v2 = _alphas->get_alpha_v2();

  auto cek = [&, this](int idx_b, int idx_a) -> bool {

    vector<Tmy_double> hsl_eta = _kernel->hit_eta(idx_b, idx_a);
    vector<Tmy_double> hsl_diff = _kernel->get_diff_Q(idx_b, idx_a);

    Tmy_double hsl_sum_v1 = this->sum_alpha_diff_Q(list_alpha_v1, hsl_diff);
    Tmy_double delta_v1 = hsl_eta[0] * hsl_sum_v1;

    Tmy_double hsl_sum_v2 = this->sum_alpha_diff_Q(list_alpha_v2, hsl_diff);
    Tmy_double delta_v2 = hsl_eta[0] * hsl_sum_v2;

    Tmy_double hsl_sum = this->sum_alpha_diff_Q(list_alpha, hsl_diff);
    Tmy_double delta = hsl_eta[0] * hsl_sum;

    Treturn_is_pass_h tmp = _alphas->is_pass(idx_b, idx_a, delta, delta_v1, delta_v2, 0);

    return  (tmp.is_pass);
  };

  srand(time(0));
  Tmy_double gmax2 = -HUGE_VAL;
  int idx_a = -1;

  Tmy_double Fb = _my_list_G->get_F(idx_b, rho.rho_v1, rho.rho_v2);
  for (int j = 0; j < _jml_data; ++j)
  {
    Tmy_double Fa = _my_list_G->get_F(j, rho.rho_v1, rho.rho_v2);
    Tmy_double abs_diff = abs(Fb - Fa);
    bool is_pass = !list_alpha->is_nol(idx_b, j);
    if (is_pass) {
      is_pass = list_alpha->is_free(j);
    }
    if ((abs_diff >= gmax2) and is_pass)
    {
      idx_a = j;
      gmax2 = abs_diff;
    }
  }


  int startIndex = (rand() % _jml_data);
  for (int i = startIndex; i < _jml_data; ++i)
  {
    if (list_alpha->is_free(i))
    {
      bool is_pass = !list_alpha->is_nol(idx_b, i);
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }
    }
  }

  if (idx_a == -1)
  {
    for (int i = 0; i < startIndex; ++i)
    {
      if (list_alpha->is_free(i))
      {
        bool is_pass = !list_alpha->is_nol(idx_b, i);
        if (cek(idx_b, i) and is_pass)
        {
          idx_a = i;
          break;
        }
      }
    }
  }

  startIndex = (rand() % _jml_data);
  if (idx_a == -1)
  {
    for (int i = startIndex; i < _jml_data; ++i)
    {
      bool is_pass = !list_alpha->is_nol(idx_b, i);
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
      }
    }
  }

  if (idx_a == -1)
  {

    for (int i = 0; i < startIndex; ++i)
    {
      bool is_pass = !list_alpha->is_nol(idx_b, i);
      if (cek(idx_b, i) and is_pass)
      {
        idx_a = i;
        break;
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

void Tmy_G::update_G(int idx_b, int idx_a, Treturn_is_pass_h tmp)
{
  _my_list_G->update_G(idx_b, idx_a, tmp.alpha.new_alpha_i, tmp.alpha.new_alpha_j);
  _my_list_G_v1->update_G(idx_b, idx_a, tmp.alpha_v1.new_alpha_i, tmp.alpha_v1.new_alpha_j);
  _my_list_G_v2->update_G(idx_b, idx_a, tmp.alpha_v2.new_alpha_i, tmp.alpha_v2.new_alpha_j);
}


