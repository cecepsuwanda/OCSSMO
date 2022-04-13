#include "Tmy_svm.h"


Tmy_svm::Tmy_svm(Tconfig *v_config) {
   _config = v_config;
   _my_alpha = new Tmy_alpha(_config);
}

Tmy_svm::~Tmy_svm() {
   delete _my_alpha;
   delete _my_kernel;
   delete _my_G;
   _model.clear();
}

Treturn_cari_alpha Tmy_svm::cari_idx_alpha()
{
   vector<Tmy_double> gmax_gmin;
   Treturn_cari_idx tmp_cari_idx = _my_G->cari_idx();

   Tmy_double diff = tmp_cari_idx.gmax + tmp_cari_idx.gmin;
   cout << tmp_cari_idx.gmax << "+" << tmp_cari_idx.gmin << "=" << diff << endl;
   bool stat = true;
   if ((diff < 1e-3) or ((tmp_cari_idx.idx_b == -1) or (tmp_cari_idx.idx_a == -1)))
   {  //(diff<1e-3) or
      stat = false;
   }

   return {stat, tmp_cari_idx.idx_b, tmp_cari_idx.idx_a};
}

bool Tmy_svm::cari_idx_a_lain(int idx_b, int *idx_alpha)
{
   //cout << "cari_idx_a_lain" << endl;
   *idx_alpha = _my_G->cari_idx_lain(idx_b);

   return (*idx_alpha != -1);
}

Tmy_double Tmy_svm::sum_alpha_diff_Q(Tmy_list_alpha* alpha, vector<Tmy_double> diff_Q)
{
   Tmy_double hasil = 0.0;

   for (int i = 0; i < diff_Q.size(); ++i)
   {
      hasil = hasil + (alpha->get_alpha(i) * diff_Q[i]);
   }

   return hasil;

}

bool Tmy_svm::take_step(int idx_b, int idx_a)
{
   if (idx_b == idx_a) {
      return false;
   } else {
      Tmy_list_G* _my_list_G = _my_G->get_list_G();
      Tmy_list_G* _my_list_G_v1 = _my_G->get_list_G_v1();
      Tmy_list_G* _my_list_G_v2 = _my_G->get_list_G_v2();

      Tmy_list_alpha* _my_list_alpha = _my_alpha->get_alpha();
      Tmy_list_alpha* _my_list_alpha_v1 = _my_alpha->get_alpha_v1();
      Tmy_list_alpha* _my_list_alpha_v2 = _my_alpha->get_alpha_v2();

      int active_size = _my_G->get_active_size();

      vector<Tmy_double> hsl_eta = _my_kernel->hit_eta(idx_b, idx_a, active_size);
      vector<Tmy_double> hsl_diff = _my_kernel->get_diff_Q(idx_b, idx_a, active_size);

      //cout<<" "<< hsl_eta[0] <<" ";

      Tmy_double hsl_sum_v1 = sum_alpha_diff_Q(_my_list_alpha_v1, hsl_diff);
      Tmy_double delta_v1 = hsl_eta[0] * hsl_sum_v1;
      //cout<<" "<< delta_v1 <<" ";

      Tmy_double hsl_sum_v2 = sum_alpha_diff_Q(_my_list_alpha_v2, hsl_diff);
      Tmy_double delta_v2 = hsl_eta[0] * hsl_sum_v2;
      //cout<<" "<< delta_v2 <<" ";

      Tmy_double hsl_sum = sum_alpha_diff_Q(_my_list_alpha, hsl_diff);
      //cout<<" "<< hsl_sum <<" ";
      Tmy_double delta = hsl_eta[0] * hsl_sum;
      //cout<<" "<< delta <<" ";

      Treturn_is_pass_h tmp = _my_alpha->is_pass(idx_b, idx_a, delta, delta_v1, delta_v2, 1);

      if (tmp.is_pass == false)
      {
         return false;
      } else {
         _my_list_G->update_G(idx_b, idx_a, tmp.alpha.new_alpha_i, tmp.alpha.new_alpha_j);
         _my_list_G_v1->update_G(idx_b, idx_a, tmp.alpha_v1.new_alpha_i, tmp.alpha_v1.new_alpha_j);
         _my_list_G_v2->update_G(idx_b, idx_a, tmp.alpha_v2.new_alpha_i, tmp.alpha_v2.new_alpha_j);
         //_rho = _my_G->update_rho(idx_b,idx_a);
         return true;
      }

   }
}


Treturn_train Tmy_svm::train(Tdataframe &df) {
   int jml_data = df.getjmlrow_svm();
   _my_kernel = new Tmy_kernel(df, _config->gamma);
   _my_alpha->init(df.getjmlrow_svm());
   _my_G = new Tmy_G(df.getjmlrow_svm(), _my_kernel, _my_alpha);
   _my_G->init();
   _rho = _my_G->update_rho();

   int iter = 0;
   int max_iter = 100;//max(10000000, jml_data > INT_MAX / 100 ? INT_MAX : 100 * jml_data);
   int counter = min(jml_data, 1000) + 1;

   bool is_alpha_changed = true;
   //(is_alpha_changed==true) and

   Treturn_train tmp_train;

   while ((iter < max_iter))
   {
      // if((iter%100)==0){
      //   cetak(".");
      // }

      if (--counter == 0)
      {
         counter = min(jml_data, 1000);
         _my_G->do_shrinking();
      }


      //cetak("iterasi ke - %d \n",iter);

      Treturn_cari_alpha hasil_cari = cari_idx_alpha();
      tmp_train.is_optimum = !hasil_cari.is_pass;

      if (hasil_cari.is_pass == false)
      {
         _my_G->reconstruct_gradient();
         _my_G->reset_active_size();
         hasil_cari = cari_idx_alpha();
         tmp_train.is_optimum = !hasil_cari.is_pass;
         if (hasil_cari.is_pass == false)
         {
            break;
         } else {
            counter = 1;
         }
      }


      if (hasil_cari.idx_a != -1.0)
      {
         cout << "iterasi ke - " << iter << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << " idx_b " << hasil_cari.idx_b << " idx_a " << hasil_cari.idx_a;
         is_alpha_changed = take_step(hasil_cari.idx_b, hasil_cari.idx_a);
         cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << endl;

         if (is_alpha_changed == false)
         {
            //     counter=1;
            int idx_a;
            bool pass = cari_idx_a_lain(hasil_cari.idx_b, &idx_a);
            if (pass == true)
            {
               cout << "iterasi ke - " << iter << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << " idx_b " << hasil_cari.idx_b << " idx_a " << idx_a;
               is_alpha_changed = take_step(hasil_cari.idx_b, idx_a);
               cout << " rho v1 " << _rho.rho_v1 << " rho v2 " << _rho.rho_v2 << endl;
               if (is_alpha_changed == false)
               {
                  counter = 1;
               }
            } else {
               cout << "break 1" << endl;
               break;
            }
         }
      } else {
         cout << "break 2" << endl;
         break;
      }

      iter = iter + 1;


   }
   //cetak("\n");

   if (iter >= max_iter)
   {
      if (_my_G->get_active_size() > jml_data)
      {
         _my_G->reconstruct_gradient();
         _my_G->reset_active_size();
      }
   }
   _rho = _my_G->update_rho();
   _my_G->reverse_swap();

   _my_kernel->clear_container();
   _my_alpha->clear_container();
   _my_G->clear_container();

   Tmy_list_alpha *list_alpha = _my_alpha->get_alpha();
   map<int, Tmy_double> alpha_sv = list_alpha->get_list_alpha_sv();
   Treturn_alpha_stat alpha_stat = list_alpha->get_stat();

   Tmy_list_alpha *list_alpha_v1 = _my_alpha->get_alpha_v1();
   map<int, Tmy_double> alpha_sv_v1 = list_alpha_v1->get_list_alpha_sv();
   Treturn_alpha_stat alpha_stat_v1 = list_alpha_v1->get_stat();

   Tmy_list_alpha *list_alpha_v2 = _my_alpha->get_alpha_v2();
   map<int, Tmy_double> alpha_sv_v2 = list_alpha_v2->get_list_alpha_sv();
   Treturn_alpha_stat alpha_stat_v2 = list_alpha_v2->get_stat();



   tmp_train.jml_iterasi = iter;
   tmp_train.jml_alpha = alpha_stat.jml_alpha;
   tmp_train.jml_alpha_v1 = alpha_stat_v1.jml_alpha;
   tmp_train.jml_alpha_v2 = alpha_stat_v2.jml_alpha;
   tmp_train.n_all_sv = alpha_stat.n_all_sv;
   tmp_train.n_sv = alpha_stat.n_sv;
   tmp_train.jml_alpha_n_sv = alpha_stat.jml_alpha_n_sv;
   tmp_train.rho_v1 = _rho.rho_v1;
   tmp_train.rho_v2 = _rho.rho_v2;


   int n_kkt = 0;

   _model.reserve(alpha_sv.size());
   for (map<int, Tmy_double>::iterator it = alpha_sv.begin(); it != alpha_sv.end(); ++it)
   {
      int idx = it->first;
      if (_my_G->is_kkt(idx, _rho) == true)
      {
         n_kkt = n_kkt + 1;
      }
      vector<string> data = df.goto_rec(idx);
      _model.push_back(data);
   }
   tmp_train.n_kkt = n_kkt;

   return tmp_train;
}

vector<string> Tmy_svm::test(Tdataframe &df)
{
   Tmy_list_alpha *list_alpha = _my_alpha->get_alpha();
   map<int, Tmy_double> alpha_sv = list_alpha->get_list_alpha_sv();
   int jml_data = df.getjmlrow_svm();
   std::vector<string> hasil;
   hasil.reserve(jml_data);
   hasil.assign(jml_data, "inside");

   vector<string> tmp_data;

   int j = 0;
   df.reset_file();
   while (!df.is_eof()) {

      tmp_data = df.get_record_svm();

      Tmy_double sum = 0.0;
      int i = 0;
      Tmy_double jml_alpha;
      for (map<int, Tmy_double>::iterator it = alpha_sv.begin(); it != alpha_sv.end(); ++it)
      {
         Tmy_double tmp = _my_kernel->kernel_function_f(tmp_data, _model[i]);
         Tmy_double alpha = it->second;
         tmp = alpha * tmp;
         sum = sum + tmp;
         i = i + 1;
      }

      sum = (sum - _rho.rho_v1) * (_rho.rho_v2 - sum);
      if (sum >= 0.0) {

      } else {
         if (sum < 0.0) {
            hasil[j] = "outside";
         }
      }

      df.next_record();
      j = j + 1;
   }

   return hasil;
}