#include "Tmy_svm.h"

Tmy_svm::Tmy_svm(Tconfig *v_config)
{
   _config = v_config;
   _my_alpha = new Tmy_alpha(_config);
}

Tmy_svm::~Tmy_svm()
{
   delete _my_alpha;
   delete _my_kernel;
   _model.clear();
   _alpha_sv.clear();
}

bool Tmy_svm::take_step(int idx_b, int idx_a)
{
   if (idx_b == idx_a)
   {
      return false;
   }
   else
   {
      cout << endl
           << idx_b << "," << idx_a << endl;

      vector<Tmy_double> hsl_eta = _my_kernel->hit_eta(idx_b, idx_a);
      // vector<Tmy_double> hsl_diff = _my_kernel->get_diff_Q(idx_b, idx_a);

      // Tmy_double hsl_sum = _my_G.sum_alpha_diff_Q(_alpha_v1, hsl_diff);
      Tmy_double hsl_sum = (_grad_v1[idx_a] - _grad_v1[idx_b]);
      Tmy_double delta_v1 = hsl_eta[0] * hsl_sum;
      Treturn_is_pass tmp_v1 = _my_alpha->is_pass(idx_b, idx_a, delta_v1, _alpha_v1);
      tmp_v1.reset();

      // hsl_sum = _my_G.sum_alpha_diff_Q(_alpha_v2, hsl_diff);
      hsl_sum = (_grad_v2[idx_a] - _grad_v2[idx_b]);
      Tmy_double delta_v2 = hsl_eta[0] * hsl_sum;
      Treturn_is_pass tmp_v2 = _my_alpha->is_pass(idx_b, idx_a, delta_v2, _alpha_v2);
      tmp_v2.reset();

      // hsl_sum = _my_G.sum_alpha_diff_Q(_alpha, hsl_diff);
      hsl_sum = (_grad[idx_a] - _grad[idx_b]);
      Tmy_double delta = hsl_eta[0] * hsl_sum;

      cout << "delta " << delta << " delta v1 " << delta_v1 << " delta v2 " << delta_v2 << endl;
      Treturn_is_pass tmp = _my_alpha->is_pass(idx_b, idx_a, delta, _alpha);

      cout << tmp_v1.is_pass << " old [" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] new [" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] " << endl;
      cout << tmp_v2.is_pass << " old [" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] new [" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] " << endl;
      cout << tmp.is_pass << " old [" << tmp.alpha_i << "," << tmp.alpha_j << "] new [" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] " << endl;

      // bool is_pass = true;
      bool is_pass = tmp.is_pass; // and (tmp_v1.is_pass or tmp_v2.is_pass);

      if (is_pass == false)
      {
         return false;
      }
      else
      {
         if ((tmp_v1 - tmp_v2) != tmp)
         {
            cout << "tak sama !!!" << endl;
            if (_my_alpha->is_pass(tmp_v1, tmp_v2, tmp) == true)
            {
               cout << "solve 1 !!!" << endl;
            }
            else
            {
               cout << "failed 1 !!!" << endl;
               is_pass = false;
            }
         }

         cout << tmp_v1.is_pass << " old [" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] new [" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] " << endl;
         cout << tmp_v2.is_pass << " old [" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] new [" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] " << endl;
         cout << tmp.is_pass << " old [" << tmp.alpha_i << "," << tmp.alpha_j << "] new [" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] " << endl;

         if (is_pass)
         {

            // if (tmp_v1.is_pass)
            // {
            _my_G.update_G(idx_b, idx_a, tmp_v1, _my_kernel, _alpha_v1, _grad_v1);

            _alpha_v1[idx_a] = tmp_v1.new_alpha_j;
            _alpha_v1[idx_b] = tmp_v1.new_alpha_i;
            _grad_v1.mv_idx(idx_a, 1);
            _grad_v1.mv_idx(idx_b, 1);
            // }
            // else
            // {
            //    tmp_v1.new_alpha_i = tmp_v1.alpha_i;
            //    tmp_v1.new_alpha_j = tmp_v1.alpha_j;
            // }

            // if (tmp_v2.is_pass)
            // {
            _my_G.update_G(idx_b, idx_a, tmp_v2, _my_kernel, _alpha_v2, _grad_v2);

            _alpha_v2[idx_a] = tmp_v2.new_alpha_j;
            _alpha_v2[idx_b] = tmp_v2.new_alpha_i;
            _grad_v2.mv_idx(idx_a, 1);
            _grad_v2.mv_idx(idx_b, 1);
            // }
            // else
            // {
            //    tmp_v2.new_alpha_i = tmp_v2.alpha_i;
            //    tmp_v2.new_alpha_j = tmp_v2.alpha_j;
            // }

            tmp.new_alpha_j = tmp_v1.new_alpha_j - tmp_v2.new_alpha_j;
            tmp.new_alpha_i = tmp_v1.new_alpha_i - tmp_v2.new_alpha_i;

            _my_G.update_G(idx_b, idx_a, tmp, _my_kernel, _alpha, _grad);

            _alpha[idx_a] = tmp.new_alpha_j;
            _alpha[idx_b] = tmp.new_alpha_i;
            _grad.mv_idx(idx_a, 1);
            _grad.mv_idx(idx_b, 1);

            vector<T_alpha_container> tmp_alpha;
            tmp_alpha.push_back(_alpha);
            tmp_alpha.push_back(_alpha_v1);
            tmp_alpha.push_back(_alpha_v2);

            // cout << " old rho [" << _rho.rho_v1 << "," << _rho.rho_v2 << "]" << endl;
            _rho = _my_G.update_rho(_my_kernel, tmp_alpha, _grad);
            // cout << " new rho [" << _rho.rho_v1 << "," << _rho.rho_v2 << "]" << endl;

            cout << tmp_v1.is_pass << " old [" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] new [" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] " << endl;
            cout << tmp_v2.is_pass << " old [" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] new [" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] " << endl;
            cout << tmp.is_pass << " old [" << tmp.alpha_i << "," << tmp.alpha_j << "] new [" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] " << endl;

            _rho_1.rho_v1 = _my_G.update_rho(_my_kernel, _alpha_v1, _grad_v1);
            _rho_1.rho_v2 = _my_G.update_rho(_my_kernel, _alpha_v2, _grad_v2);

            _my_G.set_kkt(_rho, tmp_alpha, _grad);
            _my_G.set_kkt(_rho_1.rho_v1, _alpha_v1, _grad_v1);
            _my_G.set_kkt(_rho_1.rho_v2, _alpha_v2, _grad_v2);
            return true;
         }
         else
         {
            return false;
         }
      }
   }
}

bool Tmy_svm::take_step_v(int idx_b, int idx_a, Tmy_double &rho, T_alpha_container &alpha, T_grad_container &grad)
{
   bool stat = false;
   if (idx_b == idx_a)
   {
   }
   else
   {
      // cout << " take step idx_b " << idx_b << " idx_a " << idx_a << endl;
      Tmy_double Fa = grad.dec(idx_a, rho);
      Tmy_double Fb = grad.dec(idx_b, rho);
      vector<Tmy_double> hsl_eta = _my_kernel->hit_eta(idx_b, idx_a);
      Tmy_double delta = hsl_eta[0] * (Fa - Fb);

      // cout << " delta " << delta << endl;

      Treturn_is_pass hsl = _my_alpha->is_pass(idx_b, idx_a, delta, alpha);

      // cout << hsl.is_pass << " old [" << hsl.alpha_i << "," << hsl.alpha_j << "] new [" << hsl.new_alpha_i << "," << hsl.new_alpha_j << "] " << endl;

      if (!hsl.is_pass)
      {
         // cout << " not pass !!!" << endl;
      }
      else
      {
         _my_G.update_G(idx_b, idx_a, hsl, _my_kernel, alpha, grad);

         alpha[idx_a] = hsl.new_alpha_j;
         alpha[idx_b] = hsl.new_alpha_i;
         grad.mv_idx(idx_a, 1);
         grad.mv_idx(idx_b, 1);

         rho = _my_G.update_rho(_my_kernel, alpha, grad);
         _my_G.set_kkt(rho, alpha, grad);
         stat = true;
         // cout << " pass !!!" << endl;
      }
   }
   return stat;
}

int Tmy_svm::examineExample(int idx_b)
{
   // cout << " examine Example 2 " << endl;
   int hasil = 0;
   int idx_a = -1;

   Tmy_double Gb = _grad[idx_b];
   Tmy_double obj = _grad.obj(idx_b, _rho.rho_v1, _rho.rho_v2);
   Tmy_double dec = _grad.dec(idx_b, _rho.rho_v1, _rho.rho_v2);

   vector<T_alpha_container> tmp_alpha;
   tmp_alpha.push_back(_alpha);
   tmp_alpha.push_back(_alpha_v1);
   tmp_alpha.push_back(_alpha_v2);

   vector<T_grad_container> tmp_grad;
   tmp_grad.push_back(_grad);
   tmp_grad.push_back(_grad_v1);
   tmp_grad.push_back(_grad_v2);

   // is_pass = !_my_G.is_kkt(idx_b, _rho, tmp_alpha, _grad);

   // cout << " idx_b " << idx_b << endl;
   idx_a = _my_G.cari_idx_a(idx_b, _rho, tmp_alpha, tmp_grad, _my_kernel);
   if (idx_a != -1)
   {
      // cout << " idx_a " << idx_a << " " << endl;
      bool is_pass = take_step(idx_b, idx_a);
      if (!is_pass)
      {

         vector<T_alpha_container> tmp_alpha;
         tmp_alpha.push_back(_alpha);
         tmp_alpha.push_back(_alpha_v1);
         tmp_alpha.push_back(_alpha_v2);

         vector<T_grad_container> tmp_grad;
         tmp_grad.push_back(_grad);
         tmp_grad.push_back(_grad_v1);
         tmp_grad.push_back(_grad_v2);
         cout << " cari idx_a lain 3 " << endl;
         idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, tmp_alpha, tmp_grad, _my_alpha);
         if (idx_a != -1)
         {
            cout << " idx_b " << idx_b << " idx_a " << idx_a << " ";
            is_pass = take_step(idx_b, idx_a);
            if (is_pass)
            {
               cout << "         sukses cari idx_a lain 3 " << endl;
               hasil = 1;
            }
            else
            {
            }
         }
      }
      else
      {
         cout << "           sukses 3 " << endl;
         hasil = 1;
      }
   }

   return hasil;
}

bool Tmy_svm::examineExample()
{
   // cout << " examine Example 1 " << endl;
   bool hasil = false;
   int idx_a = -1;
   int idx_b = -1;

   vector<T_alpha_container> tmp_alpha;
   tmp_alpha.push_back(_alpha);
   tmp_alpha.push_back(_alpha_v1);
   tmp_alpha.push_back(_alpha_v2);

   vector<T_grad_container> tmp_grad;
   tmp_grad.push_back(_grad);
   tmp_grad.push_back(_grad_v1);
   tmp_grad.push_back(_grad_v2);

   bool is_pass = _my_G.cari_idx(idx_b, idx_a, _rho, tmp_alpha, tmp_grad, _my_kernel);
   if (idx_a != -1)
   {
      cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
      is_pass = take_step(idx_b, idx_a);
      if (!is_pass)
      {
         vector<T_alpha_container> tmp_alpha;
         tmp_alpha.push_back(_alpha);
         tmp_alpha.push_back(_alpha_v1);
         tmp_alpha.push_back(_alpha_v2);

         vector<T_grad_container> tmp_grad;
         tmp_grad.push_back(_grad);
         tmp_grad.push_back(_grad_v1);
         tmp_grad.push_back(_grad_v2);

         cout << " cari idx_a lain 1 " << endl;
         idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, tmp_alpha, tmp_grad, _my_alpha);
         if (idx_a != -1)
         {
            cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
            is_pass = take_step(idx_b, idx_a);
            if (is_pass)
            {
               hasil = true;
            }
            else
            {
               cout << "Out 1" << endl;
            }
         }
         else
         {
            cout << "Out 2" << endl;
         }
      }
      else
      {
         hasil = true;
      }
   }
   else
   {
      if (idx_b != -1)
      {
         cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
         cout << " cari idx_a lain 2 " << endl;
         idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, tmp_alpha, tmp_grad, _my_alpha);
         if (idx_a != -1)
         {
            cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
            is_pass = take_step(idx_b, idx_a);
            if (is_pass)
            {
               hasil = true;
            }
            else
            {
               cout << "Out 3" << endl;
            }
         }
         else
         {
            cout << "Out 4" << endl;
         }
      }
      else
      {
         cout << "Out 5" << endl;
      }
   }

   // cout << " idx_b " << idx_b << " idx_a " << idx_a << " ";
   return hasil;
}

int Tmy_svm::examineExample_v(int idx_b)
{
   int tmp_v1 = 0;
   int idx_a_v1 = _my_G.cari_idx_a(idx_b, _rho_1.rho_v1, _alpha_v1, _grad_v1, _my_kernel);
   if (idx_a_v1 != -1)
   {
      bool is_pass_v1 = take_step_v(idx_b, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1);
      bool is_pass = take_step(idx_b, idx_a_v1);
      if (!is_pass_v1)
      {
         idx_a_v1 = _my_G.cari_idx_lain(idx_b, _rho_1.rho_v1, _my_kernel, _alpha_v1, _grad_v1, _my_alpha);
         if (idx_a_v1 != -1)
         {
            is_pass_v1 = take_step_v(idx_b, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1);
            is_pass = take_step(idx_b, idx_a_v1);
            if (is_pass_v1)
            {
               tmp_v1 = 1;
            }
         }
      }
      else
      {
         tmp_v1 = 1;
      }
   }

   int tmp_v2 = 0;
   int idx_a_v2 = _my_G.cari_idx_a(idx_b, _rho_1.rho_v2, _alpha_v2, _grad_v2, _my_kernel);
   if (idx_a_v2 != -1)
   {
      bool is_pass_v2 = take_step_v(idx_b, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2);
      bool is_pass = take_step(idx_b, idx_a_v2);
      if (!is_pass_v2)
      {
         idx_a_v2 = _my_G.cari_idx_lain(idx_b, _rho_1.rho_v2, _my_kernel, _alpha_v2, _grad_v2, _my_alpha);
         if (idx_a_v2 != -1)
         {
            is_pass_v2 = take_step_v(idx_b, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2);
            is_pass = take_step(idx_b, idx_a_v2);
            if (is_pass_v2)
            {
               tmp_v2 = 1;
            }
         }
      }
      else
      {
         tmp_v2 = 1;
      }
   }

   int tmp = 0;
   if ((tmp_v1 == 1) or (tmp_v2 == 1))
   {
      tmp = 1;
   }

   return tmp;
}

bool Tmy_svm::examineExample_v()
{
   int idx_b_v1 = -1;
   int idx_a_v1 = -1;
   bool examineAll_v1 = false;
   bool is_pass_v1 = _my_G.cari_idx(idx_b_v1, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1, _my_kernel);
   if (idx_a_v1 != -1)
   {
      // cout << " idx_b " << idx_b << " idx_a " << idx_a << endl;
      is_pass_v1 = take_step_v(idx_b_v1, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1);
      bool is_pass = take_step(idx_b_v1, idx_a_v1);
      if (!is_pass_v1)
      {
         idx_a_v1 = _my_G.cari_idx_lain(idx_b_v1, _rho_1.rho_v1, _my_kernel, _alpha_v1, _grad_v1, _my_alpha);
         if (idx_a_v1 != -1)
         {
            is_pass_v1 = take_step_v(idx_b_v1, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1);
            is_pass = take_step(idx_b_v1, idx_a_v1);
            if (is_pass_v1)
            {
               examineAll_v1 = true;
            }
            else
            {
               cout << "v1 Out 1" << endl;
            }
         }
         else
         {
            cout << "v1 Out 2" << endl;
         }
      }
      else
      {
         examineAll_v1 = true;
      }
   }
   else
   {
      if (idx_b_v1 != -1)
      {
         idx_a_v1 = _my_G.cari_idx_lain(idx_b_v1, _rho_1.rho_v1, _my_kernel, _alpha_v1, _grad_v1, _my_alpha);
         if (idx_a_v1 != -1)
         {
            is_pass_v1 = take_step_v(idx_b_v1, idx_a_v1, _rho_1.rho_v1, _alpha_v1, _grad_v1);
            bool is_pass = take_step(idx_b_v1, idx_a_v1);
            if (is_pass_v1)
            {
               examineAll_v1 = true;
            }
            else
            {
               cout << "v1 Out 3" << endl;
            }
         }
         else
         {
            cout << "v1 Out 4" << endl;
         }
      }
      else
      {
         cout << "v1 Out 5" << endl;
      }
   }

   int idx_b_v2 = -1;
   int idx_a_v2 = -1;
   bool examineAll_v2 = false;
   bool is_pass_v2 = _my_G.cari_idx(idx_b_v2, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2, _my_kernel);
   if (idx_a_v2 != -1)
   {
      // cout << " idx_b " << idx_b << " idx_a " << idx_a << endl;
      is_pass_v2 = take_step_v(idx_b_v2, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2);
      bool is_pass = take_step(idx_b_v2, idx_a_v2);
      if (!is_pass_v2)
      {
         idx_a_v2 = _my_G.cari_idx_lain(idx_b_v2, _rho_1.rho_v2, _my_kernel, _alpha_v2, _grad_v2, _my_alpha);
         if (idx_a_v2 != -1)
         {
            is_pass_v2 = take_step_v(idx_b_v2, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2);
            is_pass = take_step(idx_b_v2, idx_a_v2);
            if (is_pass_v2)
            {
               examineAll_v2 = true;
            }
            else
            {
               cout << "v2 Out 1" << endl;
            }
         }
         else
         {
            cout << "v2 Out 2" << endl;
         }
      }
      else
      {
         examineAll_v2 = true;
      }
   }
   else
   {
      if (idx_b_v2 != -1)
      {
         idx_a_v2 = _my_G.cari_idx_lain(idx_b_v2, _rho_1.rho_v2, _my_kernel, _alpha_v2, _grad_v2, _my_alpha);
         if (idx_a_v2 != -1)
         {
            is_pass_v2 = take_step_v(idx_b_v2, idx_a_v2, _rho_1.rho_v2, _alpha_v2, _grad_v2);
            bool is_pass = take_step(idx_b_v2, idx_a_v2);
            if (is_pass_v2)
            {
               examineAll_v2 = true;
            }
            else
            {
               cout << "v2 Out 3" << endl;
            }
         }
         else
         {
            cout << "v2 Out 4" << endl;
         }
      }
      else
      {
         cout << "v2 Out 5" << endl;
      }
   }
   return (examineAll_v1 or examineAll_v2);
}

Treturn_train Tmy_svm::train(Tdataframe &df)
{
   cout << " masuk " << endl;
   int jml_data = df.getjmlrow_svm();
   _my_kernel = new Tmy_kernel(df, _config->gamma);
   _my_alpha->init(jml_data, _alpha, _alpha_v1, _alpha_v2);
   _my_G.init(jml_data, _my_kernel, _alpha_v1, _grad_v1);
   _my_G.init(jml_data, _my_kernel, _alpha_v2, _grad_v2);
   _my_G.init(jml_data, _my_kernel, _alpha, _grad);
   vector<T_alpha_container> tmp_alpha;
   tmp_alpha.push_back(_alpha);
   tmp_alpha.push_back(_alpha_v1);
   tmp_alpha.push_back(_alpha_v2);
   _rho = _my_G.update_rho(_my_kernel, tmp_alpha, _grad);
   _rho_1.rho_v1 = _my_G.update_rho(_my_kernel, _alpha_v1, _grad_v1);
   _rho_1.rho_v2 = _my_G.update_rho(_my_kernel, _alpha_v2, _grad_v2);

   _my_G.set_kkt(_rho, tmp_alpha, _grad);
   _my_G.set_kkt(_rho_1.rho_v1, _alpha_v1, _grad_v1);
   _my_G.set_kkt(_rho_1.rho_v2, _alpha_v2, _grad_v2);

   int max_iter = jml_data * 10;
   int counter = min(jml_data, 1000) + 1;
   bool stop_iter = false;
   int jml_out = 0;

   cout << " masuk " << endl;
   for (int i = 0; i < jml_data; ++i)
   {
      cout << i << setw(15) << _alpha_v1[i] << setw(15) << _alpha_v2[i] << setw(15) << _alpha[i] << setw(15) << _grad[i] << setw(15) << _rho.rho_v1 << setw(15) << _rho.rho_v2 << setw(15) << _grad.get_kkt(i) << endl;

      // cout << _grad_v1[i] << setw(15) << _grad_v2[i] << setw(15) << _grad[i] << setw(15) << _rho.rho_v1 << setw(15) << _rho.rho_v2 << endl;
   }

   int iter = 0;
   bool examineAll = true;
   while ((iter < max_iter) and !stop_iter)
   {

      if (examineAll)
      {
         int i = 0;
         bool ulangi = true;
         while ((i < counter) and ulangi)
         {
            // cout << " Iter : " << (iter + 1) << endl;
            ulangi = examineExample();
            i = i + 1;
            iter = iter + 1;
         }
         if (!ulangi)
         {
            jml_out = jml_out + 1;
            if (jml_out == 10)
            {
               cout << " jml out " << jml_out << endl;
               stop_iter = true;
            }
            else
            {
               examineAll = false;
            }
         }
      }
      else
      {
         vector<int> rand_idx = _grad.get_rand_idx();
         int jml_pass = 0;
         for (size_t i = 0; i < jml_data; i++)
         {
            // cout << " Iter : " << (iter + 1) << endl;
            jml_pass = jml_pass + examineExample(rand_idx[i]);
            iter = iter + 1;
         }

         if (jml_pass > 0)
         {
            examineAll = true;
         }
         else
         {
            cout << " jml out " << jml_out << endl;
            stop_iter = true;
         }
      }
   }

   tmp_alpha.clear();
   tmp_alpha.push_back(_alpha);
   tmp_alpha.push_back(_alpha_v1);
   tmp_alpha.push_back(_alpha_v2);
   _rho = _my_G.update_rho(_my_kernel, tmp_alpha, _grad);
   _rho_1.rho_v1 = _my_G.update_rho(_my_kernel, _alpha_v1, _grad_v1);
   _rho_1.rho_v2 = _my_G.update_rho(_my_kernel, _alpha_v2, _grad_v2);

   _my_G.set_kkt(_rho, tmp_alpha, _grad);
   _my_G.set_kkt(_rho_1.rho_v1, _alpha_v1, _grad_v1);
   _my_G.set_kkt(_rho_1.rho_v2, _alpha_v2, _grad_v2);

   cout << endl;
   for (int i = 0; i < jml_data; ++i)
   {
      cout << i << setw(15) << _alpha_v1[i] << setw(15) << _alpha_v2[i] << setw(15) << _alpha[i] << setw(15) << _grad[i] << setw(15) << _rho.rho_v1 << setw(15) << _rho.rho_v2 << setw(15) << _grad.get_kkt(i) << endl;
   }

   Treturn_train tmp_train;

   tmp_train.jml_iterasi = iter;
   tmp_train.jml_alpha = _alpha.sum();
   tmp_train.jml_alpha_v1 = _alpha_v1.sum_if([](Tmy_double val) -> bool
                                             { return val != 0.0; });
   tmp_train.jml_alpha_v2 = _alpha_v2.sum_if([](Tmy_double val) -> bool
                                             { return val != 0.0; });
   tmp_train.n_all_sv = _alpha.n_all_sv();
   tmp_train.n_sv = _alpha.n_sv();

   tmp_train.n_all_sv_v1 = _alpha_v1.n_all_sv();
   tmp_train.n_sv_v1 = _alpha_v1.n_sv() + _alpha_v1.n_sv(_alpha_v2.lb(), _alpha_v2.ub());

   tmp_train.n_all_sv_v2 = _alpha_v2.n_all_sv();
   tmp_train.n_sv_v2 = _alpha_v2.n_sv() + _alpha_v2.n_sv(_alpha_v1.lb(), _alpha_v1.ub());

   tmp_train.rho_v1 = _rho.rho_v1;
   tmp_train.rho_v2 = _rho.rho_v2;

   int n_kkt = 0;
   int n_kkt_v1 = 0;
   int n_kkt_v2 = 0;
   int i = 0;
   _alpha_sv.reserve(tmp_train.n_all_sv);
   for (int idx = 0; idx < jml_data; ++idx)
   {
      if (_alpha[idx] != 0.0)
      {
         if (_my_G.is_kkt(idx, _rho, tmp_alpha, _grad) == true)
         {
            n_kkt = n_kkt + 1;
         }
         vector<string> data = df.goto_rec(idx);
         _alpha_sv.push_back(_alpha[idx]);
         _model.insert(pair<int, vector<string>>(i, data));
         i = i + 1;
      }

      if (_alpha_v1[idx] != 0.0)
      {
         if (_my_G.is_kkt(idx, _rho_1.rho_v1, _alpha_v1, _grad_v1) == true)
         {
            n_kkt_v1 = n_kkt_v1 + 1;
         }
      }
      if (_alpha_v2[idx] != 0.0)
      {
         if (_my_G.is_kkt(idx, _rho_1.rho_v2, _alpha_v2, _grad_v2) == true)
         {
            n_kkt_v2 = n_kkt_v2 + 1;
         }
      }
   }
   tmp_train.n_kkt = n_kkt;
   tmp_train.n_kkt_v1 = n_kkt_v1;
   tmp_train.n_kkt_v2 = n_kkt_v2;

   _my_kernel->clear_container();
   // _my_alpha->clear_container();
   // _my_G->clear_container();

   return tmp_train;
}

vector<string> Tmy_svm::test(Tdataframe &df)
{

   int jml_data = df.getjmlrow_svm();
   std::vector<string> hasil;
   hasil.reserve(jml_data);
   hasil.assign(jml_data, "inside");

   vector<string> tmp_data;

   int j = 0;
   df.reset_file();
   while (!df.is_eof())
   {

      tmp_data = df.get_record_svm();

      Tmy_double sum = 0.0;
      for (map<int, vector<string>>::iterator it = _model.begin(); it != _model.end(); ++it)
      {
         Tmy_double tmp = _my_kernel->kernel_function_f(tmp_data, it->second);
         Tmy_double alpha = _alpha_sv[it->first];
         tmp = alpha * tmp;
         sum = sum + tmp;
      }

      sum = (sum - _rho.rho_v1) * (_rho.rho_v2 - sum);
      if (sum >= 0.0)
      {
      }
      else
      {
         if (sum < 0.0)
         {
            hasil[j] = "outside";
         }
      }

      df.next_record();
      j = j + 1;
   }

   return hasil;
}