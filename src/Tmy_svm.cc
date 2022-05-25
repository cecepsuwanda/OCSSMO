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

int Tmy_svm::take_step(int idx_b, int idx_a)
{
   if (idx_b == idx_a)
   {
      return 0;
   }
   else
   {

      auto tambah = [](Treturn_is_pass & v1, Treturn_is_pass & v2, Treturn_is_pass v)->bool
      {
         bool is_pass = false;
         Treturn_is_pass coba_v1;
         Treturn_is_pass coba_v2;

         if (v.new_alpha_i > 0.0)
         {
            if (v1.is_pass and v2.is_pass) {
               coba_v2 = v2;
               coba_v2 = 0.0;
               coba_v1 = v1;
               coba_v1 = v.new_alpha_i;
            } else {
               if (v2.is_pass) {
                  coba_v1 = v1;
                  coba_v2 = v2;
                  Tmy_double delta_i = v1.new_alpha_i - v.new_alpha_i;
                  coba_v2 = delta_i;
               } else {
                  if (v1.is_pass) {
                     coba_v1 = v1;
                     coba_v2 = v2;
                     Tmy_double delta_i = v2.new_alpha_i + v.new_alpha_i;
                     coba_v1 = delta_i;
                  }
               }
            }
         } else{
            if (v.new_alpha_i < 0.0)
            {
               if (v1.is_pass and v2.is_pass) {
                  coba_v1 = v1;
                  coba_v1 = 0.0;
                  coba_v2 = v2;
                  coba_v2 = -1.0 * v.new_alpha_i;
               } else {
                  if (v1.is_pass) {
                     coba_v1 = v1;
                     coba_v2 = v2;
                     Tmy_double delta_i = v2.new_alpha_i + v.new_alpha_i;
                     coba_v1 = delta_i;
                  } else {
                     if (v2.is_pass) {
                        coba_v1 = v1;
                        coba_v2 = v2;
                        Tmy_double delta_i = v1.new_alpha_i - v.new_alpha_i;
                        coba_v2 = delta_i;
                     }
                  }
               }
            }
         }


         if ( (coba_v1 - coba_v2) == v )
         {
            v1 = coba_v1;
            v2 = coba_v2;
            is_pass = true;
         } else {
            v1.new_alpha_i = v1.alpha_i;
            v1.new_alpha_j = v1.alpha_j;
            v2.new_alpha_i = v2.alpha_i;
            v2.new_alpha_j = v2.alpha_j;
         }

         return is_pass;
      };

      cout << endl << idx_b << "," << idx_a << endl;

      vector<Tmy_double> hsl_eta = _my_kernel->hit_eta(idx_b, idx_a);
      vector<Tmy_double> hsl_diff = _my_kernel->get_diff_Q(idx_b, idx_a);

      Tmy_double hsl_sum = 0.0;
      //Tmy_double delta = 0.0;

      hsl_sum = _my_G.sum_alpha_diff_Q(_alpha_v1, hsl_diff);
      Tmy_double delta_v1 = hsl_eta[0] * hsl_sum;
      cout << "delta v1 " << delta_v1 << endl;
      Treturn_is_pass tmp_v1 = _my_alpha->is_pass(idx_b, idx_a, delta_v1, _alpha_v1);

      hsl_sum = _my_G.sum_alpha_diff_Q(_alpha_v2, hsl_diff);
      Tmy_double delta_v2 = hsl_eta[0] * hsl_sum;
      cout << "delta v2 " << delta_v2 << endl;
      Treturn_is_pass tmp_v2 = _my_alpha->is_pass(idx_b, idx_a, delta_v2, _alpha_v2);

      hsl_sum = _my_G.sum_alpha_diff_Q(_alpha, hsl_diff);
      Tmy_double delta = hsl_eta[0] * hsl_sum;
      cout << "delta " << delta << " delta_v1-delta_v2 " << (delta_v1 - delta_v2) << endl;
      Treturn_is_pass tmp = _my_alpha->is_pass(idx_b, idx_a, delta, _alpha);

      cout << tmp_v1.is_pass << " old [" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] new [" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] " << endl;
      cout << tmp_v2.is_pass << " old [" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] new [" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] " << endl;
      cout << tmp.is_pass << " old [" << tmp.alpha_i << "," << tmp.alpha_j << "] new [" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] " << endl;

      bool is_pass = tmp.is_pass and (tmp_v1.is_pass or tmp_v2.is_pass);


      // if (is_pass == true)
      // {

      //    if ((tmp.new_alpha_i < 0.0) or (tmp.new_alpha_j < 0.0) ) {
      //       is_pass = ((tmp_v2.new_alpha_i > 1e-3) or (tmp_v2.new_alpha_j > 1e-3)) and ((tmp_v1.new_alpha_i > 1e-3) or (tmp_v1.new_alpha_j > 1e-3));

      //    } else {
      //       if (tmp.new_alpha_i > 0.0)
      //       {
      //          is_pass = (tmp_v1.new_alpha_i > 1e-3) or (tmp_v1.new_alpha_j > 1e-3);
      //       }
      //    }
      // }


      if (is_pass == false)
      {
         return 0;
      }
      else
      {
         if ((tmp_v1 - tmp_v2) != tmp)
         {
            cout << "tak sama !!!" << endl;
            if (tambah(tmp_v1, tmp_v2, tmp) == true)
            {
               cout << "solve 1 !!!" << endl;
            } else {
               cout << "failed 1 !!!" << endl;
               is_pass = false;
            }
         }
         if (is_pass) {
            _my_G.update_G(idx_b, idx_a, tmp, _my_kernel, _alpha_v1, _grad_v1);

            _alpha_v1[idx_a] = tmp_v1.new_alpha_j;
            _alpha_v1[idx_b] = tmp_v1.new_alpha_i;

            _my_G.update_G(idx_b, idx_a, tmp, _my_kernel, _alpha_v2, _grad_v2);

            _alpha_v2[idx_a] = tmp_v2.new_alpha_j;
            _alpha_v2[idx_b] = tmp_v2.new_alpha_i;


            tmp.new_alpha_j = tmp_v1.new_alpha_j - tmp_v2.new_alpha_j;
            tmp.new_alpha_i = tmp_v1.new_alpha_i - tmp_v2.new_alpha_i;

            _my_G.update_G(idx_b, idx_a, tmp, _my_kernel, _alpha, _grad);

            _alpha[idx_a] = tmp.new_alpha_j;
            _alpha[idx_b] = tmp.new_alpha_i;
            _grad.mv_idx(idx_a, 1);
            _grad.mv_idx(idx_b, 1);

            cout << tmp_v1.is_pass << " old [" << tmp_v1.alpha_i << "," << tmp_v1.alpha_j << "] new [" << tmp_v1.new_alpha_i << "," << tmp_v1.new_alpha_j << "] " << endl;
            cout << tmp_v2.is_pass << " old [" << tmp_v2.alpha_i << "," << tmp_v2.alpha_j << "] new [" << tmp_v2.new_alpha_i << "," << tmp_v2.new_alpha_j << "] " << endl;
            cout << tmp.is_pass << " old [" << tmp.alpha_i << "," << tmp.alpha_j << "] new [" << tmp.new_alpha_i << "," << tmp.new_alpha_j << "] " << endl;

            _rho = _my_G.update_rho(_my_kernel, _alpha, _alpha_v1, _alpha_v2);
            return 1;
         } else {
            return 0;
         }


      }
   }
}

int Tmy_svm::examineExample(int idx_b)
{
   int hasil = 0;
   int idx_a = -1;

   Tmy_double Gb = _grad[idx_b];
   Tmy_double obj = _grad.obj(idx_b, _rho.rho_v1, _rho.rho_v2);
   Tmy_double dec = _grad.dec(idx_b, _rho.rho_v1, _rho.rho_v2);

   // cout << " [" << Gb << "," << (Gb-_rho.rho_v1) << "," << (_rho.rho_v2-Gb) << "] " << endl;

   bool is_pass = true;

   // if (is_pass)
   // {
   //    is_pass = ((dec < -1e-3) and ((_alpha[idx_b] < _alpha.ub()) or (_alpha[idx_b] > _alpha.lb()))) or ((dec > 1e-3) and ((_alpha[idx_b] < 0.0) or (_alpha[idx_b] > 0.0)) );
   // }

   // if (is_pass)
   // {
   //    is_pass = ((dec < -1e-3) and (((_alpha_v1[idx_b]<_alpha_v1.ub()) and (_alpha_v2[idx_b]>_alpha_v2.lb())) or ((_alpha_v1[idx_b]>_alpha_v1.lb()) and (_alpha_v2[idx_b]<_alpha_v2.ub())))) or ((dec > 1e-3) and ((_alpha_v1[idx_b]>_alpha_v1.lb()) or (_alpha_v2[idx_b]>_alpha_v2.lb())));
   // }


   // if ((_alpha_v1.is_nol(idx_b) and _alpha_v2.is_nol(idx_b)) and (dec >= 1e-3))
   // {
   //    is_pass = false;
   // }
   // else
   // {
      if (((_alpha_v1.is_sv(idx_b) and _alpha_v2.is_nol(idx_b)) or (_alpha_v1.is_nol(idx_b) and _alpha_v2.is_sv(idx_b))) and (abs(dec) < 1e-3))
      {
         is_pass = false;
      }
      else
      {
         if (((_alpha_v1.is_ub(idx_b) and _alpha_v2.is_nol(idx_b)) or (_alpha_v1.is_nol(idx_b) and _alpha_v2.is_ub(idx_b)))  and (dec <= -1e-3))
         {
            is_pass = false;
         }
      }
   // }


   if (is_pass)
   {
      idx_a = _my_G.cari_idx_a(idx_b, _rho, _alpha, _alpha_v1, _alpha_v2, _grad, _my_kernel);
      if (idx_a != -1)
      {
         cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
         hasil = take_step(idx_b, idx_a);
         if (hasil == 0)
         {
            idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, _alpha, _alpha_v1, _alpha_v2, _grad, _my_alpha);
            if (idx_a != -1)
            {
               cout << " idx_b " << idx_b << " idx_a " << idx_a << " ";
               hasil = take_step(idx_b, idx_a);
               if (hasil == 1)
               {
                  cout << "         sukses cari idx_a lain 1 " << endl;
               }
            }
         } else {
            cout << "           sukses 1 " << endl;
         }
      }
      else
      {
         if (idx_b != -1)
         {
            cout << " idx_b " << idx_b << " idx_a " << idx_a << " " << endl;
            idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, _alpha, _alpha_v1, _alpha_v2, _grad, _my_alpha);
            if (idx_a != -1)
            {
               cout << " idx_b " << idx_b << " idx_a " << idx_a << " ";
               hasil = take_step(idx_b, idx_a);
               if (hasil == 1)
               {
                  cout << "          sukses cari idx_a lain 2 " << endl;
               }
            }
         }
      }
   }

   return hasil;
}

int Tmy_svm::examineExample()
{
   int hasil = 0;
   int idx_a = -1;
   int idx_b = -1;

   bool is_pass = _my_G.cari_idx(idx_b, idx_a, _rho, _alpha, _alpha_v1, _alpha_v2, _grad, _my_kernel);
   if (idx_a != -1)
   {
      hasil = take_step(idx_b, idx_a);
      if (hasil == 0)
      {
         //cout << " cari idx_a lain ";
         idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, _alpha, _alpha_v1, _alpha_v2, _grad, _my_alpha);
         if (idx_a != -1)
         {
            hasil = take_step(idx_b, idx_a);
         }
      }
   }
   else
   {
      if (idx_b != -1)
      {
         //cout << " cari idx_a lain ";
         idx_a = _my_G.cari_idx_lain(idx_b, _rho, _my_kernel, _alpha, _alpha_v1, _alpha_v2, _grad, _my_alpha);
         if (idx_a != -1)
         {
            hasil = take_step(idx_b, idx_a);
         }
      }
   }

   //cout << " idx_b " << idx_b << " idx_a " << idx_a << " ";
   return hasil;
}

Treturn_train Tmy_svm::train(Tdataframe & df)
{
   int jml_data = df.getjmlrow_svm();
   _my_kernel = new Tmy_kernel(df, _config->gamma);
   _my_alpha->init(jml_data, _alpha, _alpha_v1, _alpha_v2);
   _my_G.init(jml_data, _my_kernel, _alpha_v1, _grad_v1);
   _my_G.init(jml_data, _my_kernel, _alpha_v2, _grad_v2);
   _my_G.init(jml_data, _my_kernel, _alpha, _grad);
   _rho = _my_G.update_rho(_my_kernel, _alpha, _alpha_v1, _alpha_v2);

   int iter = 0;
   int max_iter = 100;//max(10000000, jml_data > INT_MAX / 100 ? INT_MAX : 100 * jml_data);

   int is_alpha_changed = 1;
   bool examineAll = true;

   Treturn_train tmp_train;

   while ((iter < max_iter) and (is_alpha_changed > 0))
   {
      // if((iter%100)==0){
      //   cetak(".");
      // }
      is_alpha_changed = 0;
      cout << endl << "iterasi ke - " << iter << endl;
      if (examineAll)
      {
         for (int idx_b = 0; idx_b < jml_data; idx_b++)
         {
            cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2;
            is_alpha_changed += examineExample(idx_b);
            cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << endl;
            // if (is_alpha_changed != (idx_b + 1))
            // {
            //    break;
            // }
         }
         // if (is_alpha_changed == 0)
         // {
         //    is_alpha_changed = 1;
         //    examineAll = !examineAll;
         // }
      }
      else
      {
         for (int idx_b = 0; idx_b < jml_data; idx_b++)
         {
            bool is_pass = true;
            is_pass = (!_alpha_v1.is_ub(idx_b) and _alpha_v2.is_nol(idx_b)) and (_alpha_v1.is_nol(idx_b) and !_alpha_v2.is_ub(idx_b));
            if (is_pass)
            {
               cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2;
               is_alpha_changed += examineExample(idx_b);
               cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << endl;
               // if (is_alpha_changed != 0)
               // {
               //    break;
               // }
            }
         }
         // if (is_alpha_changed != 0)
         // {
         //    examineAll = !examineAll;
         // }
      }
      examineAll = !examineAll;
      iter = iter + 1;
   }
   // cetak("\n");

   // _alpha.nol_kan();
   // _alpha_v1.nol_kan();
   // _alpha_v2.nol_kan();
   _rho = _my_G.update_rho(_my_kernel, _alpha, _alpha_v1, _alpha_v2);

   tmp_train.jml_iterasi = iter;
   tmp_train.jml_alpha = _alpha.sum();
   tmp_train.jml_alpha_v1 = _alpha_v1.sum_if([](Tmy_double val) -> bool
   { return val != 0.0; });
   tmp_train.jml_alpha_v2 = _alpha_v2.sum_if([](Tmy_double val) -> bool
   { return val != 0.0; });
   tmp_train.n_all_sv = _alpha.n_all_sv();
   tmp_train.n_sv = _alpha.n_sv();
   tmp_train.rho_v1 = _rho.rho_v1;
   tmp_train.rho_v2 = _rho.rho_v2;

   int n_kkt = 0;
   int i = 0;
   _alpha_sv.reserve(tmp_train.n_all_sv);
   for (int idx = 0; idx < jml_data; ++idx)
   {
      if (_alpha[idx] != 0.0)
      {
         if (_my_G.is_kkt(idx, _rho, _alpha, _alpha_v1, _alpha_v2, _grad) == true)
         {
            n_kkt = n_kkt + 1;
         }
         vector<string> data = df.goto_rec(idx);
         _alpha_sv.push_back(_alpha[idx]);
         _model.insert(pair<int, vector<string>>(i, data));
         i = i + 1;
      }
   }
   tmp_train.n_kkt = n_kkt;

   _my_kernel->clear_container();
   // _my_alpha->clear_container();
   // _my_G->clear_container();

   return tmp_train;
}

vector<string> Tmy_svm::test(Tdataframe & df)
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