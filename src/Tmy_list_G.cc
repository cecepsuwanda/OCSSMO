#include "Tmy_list_G.h"


Tmy_list_G::Tmy_list_G(int jml_data, Tmy_kernel *kernel, Tmy_list_alpha *alpha)
{
    _jml_data = jml_data;
    _kernel = kernel;
    _alpha = alpha;
    _arr_G.reserve(_jml_data);
    _arr_G.assign(_jml_data, 0.0);    
    init();
}

Tmy_list_G::~Tmy_list_G()
{
    clear_container();
}

void Tmy_list_G::clear_container()
{
    _arr_G.clear();       
}

void Tmy_list_G::init()
{
    //cetak("Start init G : \n");
    int i = 0;
    for (int idx=0;idx<_jml_data;idx++)
    {
        // if((i%100)==0){
        //    cetak(".");
        //   }
        Tmy_double alpha = _alpha->get_alpha(idx);
        vector<Tmy_double> data = _kernel->get_Q(idx);
        if(alpha!=0.0){
           for (int j = 0; j < _jml_data; ++j)
           {
            _arr_G.at(j) = _arr_G.at(j) + (alpha * data[j]);
           }
        }

        i = i + 1;

    }
    //cetak("\nEnd init G : \n");

}

void Tmy_list_G::update_G(int idx_b, int idx_a, Tmy_double new_alpha_b, Tmy_double new_alpha_a)
{

    Tmy_double alpha_a = _alpha->get_alpha(idx_a);
    Tmy_double alpha_b = _alpha->get_alpha(idx_b);

    Tmy_double delta_1 = new_alpha_b - alpha_b;
    Tmy_double delta_2 = new_alpha_a - alpha_a;

    vector<Tmy_double> data_a = _kernel->get_Q(idx_a);
    vector<Tmy_double> data_b = _kernel->get_Q(idx_b);

    for (int i = 0; i < _jml_data; ++i)
    {
        _arr_G.at(i) = _arr_G.at(i) + ((data_b[i] * delta_1) + (data_a[i] * delta_2));
    }
    
    _alpha->update_alpha(idx_a, new_alpha_a);
    _alpha->update_alpha(idx_b, new_alpha_b);
    
}

Tmy_double Tmy_list_G::get_G(int idx)
{
    return _arr_G.at(idx);
}

Tmy_double Tmy_list_G::get_obj(int idx, Tmy_double rho1, Tmy_double rho2)
{
    Tmy_double G = _arr_G.at(idx);
    Tmy_double tmp_F = min((G - rho1), (rho2 - G));
    return tmp_F;
}

Tmy_double Tmy_list_G::get_dec(int idx, Tmy_double rho1, Tmy_double rho2)
{
    Tmy_double G = _arr_G.at(idx);
    Tmy_double tmp_F = (G - rho1)*(rho2 - G);
    return tmp_F;
}


