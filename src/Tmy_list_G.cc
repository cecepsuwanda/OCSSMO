#include "Tmy_list_G.h"


Tmy_list_G::Tmy_list_G(int jml_data, Tmy_kernel *kernel, Tmy_list_alpha *alpha)
{
    _jml_data = jml_data;
    _active_size = jml_data;
    _kernel = kernel;
    _alpha = alpha;
    _arr_G.reserve(_jml_data);
    _arr_G.assign(_jml_data, 0.0);
    _arr_G_bar.reserve(_jml_data);
    _arr_G_bar.assign(_jml_data, 0.0);
    _active_set.reserve(_jml_data);
    _active_set.assign(_jml_data, 0);
    init();
}

Tmy_list_G::~Tmy_list_G()
{
    clear_container();
}

void Tmy_list_G::clear_container()
{
    _arr_G.clear();
    _arr_G_bar.clear();
    _active_set.clear();
}

void Tmy_list_G::init()
{
    //cetak("Start init G : \n");

    for (int i = 0; i < _jml_data; ++i)
    {
        _active_set[i] = i;
    }

    vector<int> idx_alpha = _alpha->get_list_lb_ub(3);

    int i = 0;
    for (auto& idx : idx_alpha)
    {
        // if((i%100)==0){
        //    cetak(".");
        //   }
        Tmy_double alpha = _alpha->get_alpha(idx);
        vector<Tmy_double> data = _kernel->get_Q(idx, _jml_data);
        for (int j = 0; j < _jml_data; ++j)
        {
            _arr_G.at(j) = _arr_G.at(j) + (alpha * data[j]);
        }

        if ( _alpha->is_upper_bound(idx) == true)
        {
            for (int j = 0; j < _jml_data; ++j)
            {
                _arr_G_bar[j] = _arr_G_bar[j] + (_alpha->get_ub() * data[j]);
            }
        }

        if ( _alpha->is_lower_bound(idx) == true)
        {
            for (int j = 0; j < _jml_data; ++j)
            {
                _arr_G_bar[j] = _arr_G_bar[j] + (_alpha->get_lb() * data[j]);
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

    vector<Tmy_double> data_a = _kernel->get_Q(idx_a, _active_size);
    vector<Tmy_double> data_b = _kernel->get_Q(idx_b, _active_size);

    for (int i = 0; i < _active_size; ++i)
    {
        _arr_G.at(i) = _arr_G.at(i) + ((data_b[i] * delta_1) + (data_a[i] * delta_2));
    }

    bool is_alpha_a_ub = _alpha->is_upper_bound(idx_a);
    bool is_alpha_b_ub = _alpha->is_upper_bound(idx_b);
    bool is_alpha_a_lb = _alpha->is_lower_bound(idx_a);
    bool is_alpha_b_lb = _alpha->is_lower_bound(idx_b);
    _alpha->update_alpha(idx_a, new_alpha_a);
    _alpha->update_alpha(idx_b, new_alpha_b);

    if (is_alpha_a_ub != _alpha->is_upper_bound(idx_a))
    {
        vector<Tmy_double> data_a = _kernel->get_Q(idx_a, _jml_data);
        if (is_alpha_a_ub == true)
        {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] - (data_a[i] * _alpha->get_ub());
            }
        } else {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] + (data_a[i] * _alpha->get_ub());
            }
        }
    }

    if (is_alpha_b_ub != _alpha->is_upper_bound(idx_b))
    {
        vector<Tmy_double> data_b = _kernel->get_Q(idx_b, _jml_data);
        if (is_alpha_b_ub == true)
        {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] - (data_b[i] * _alpha->get_ub());
            }
        } else {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] + (data_b[i] * _alpha->get_ub());
            }
        }
    }


    if (is_alpha_a_lb != _alpha->is_lower_bound(idx_a))
    {
        vector<Tmy_double> data_a = _kernel->get_Q(idx_a, _jml_data);
        if (is_alpha_a_lb == true)
        {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] - (data_a[i] * _alpha->get_lb());
            }
        } else {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] + (data_a[i] * _alpha->get_lb());
            }
        }
    }

    if (is_alpha_b_lb != _alpha->is_lower_bound(idx_b))
    {
        vector<Tmy_double> data_b = _kernel->get_Q(idx_b, _jml_data);
        if (is_alpha_b_lb == true)
        {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] - (data_b[i] * _alpha->get_lb());
            }
        } else {
            for (int i = 0; i < _jml_data; ++i)
            {
                _arr_G_bar[i] = _arr_G_bar[i] + (data_b[i] * _alpha->get_lb());
            }
        }
    }
}

Tmy_double Tmy_list_G::get_G(int idx)
{
    return _arr_G.at(idx);
}

Tmy_double Tmy_list_G::get_F(int idx,Tmy_double rho1,Tmy_double rho2)
{
   Tmy_double G = _arr_G.at(idx);
   Tmy_double tmp_F = min((G - rho1), (rho2 - G));
   return tmp_F; 
}

void Tmy_list_G::swap_index(int i, int j)
{
    _alpha->swap_index(i, j);
    swap(_arr_G.at(i), _arr_G.at(j));
    swap(_active_set[i], _active_set[j]);
    swap(_arr_G_bar[i], _arr_G_bar[j]);
}

void Tmy_list_G::reverse_swap()
{
    for (int i = 0; i < _active_set.size(); ++i)
    {
        if (_active_set[i] != i)
        {
            _alpha->swap_index(i, _active_set[i]);
            swap(_arr_G.at(i), _arr_G.at(_active_set[i]));
            _active_set[_active_set[i]] = _active_set[i];
            _active_set[i] = i;
        }
    }
}

void Tmy_list_G::reconstruct_gradient()
{
    if (_active_size != _jml_data) {
        cout << "start reconstruct gradient" << endl;

        int i, j;
        int nr_free = 0;

        for (j = _active_size; j < _jml_data; j++)
            _arr_G.at(j) = _arr_G_bar[j];

        for (j = 0; j < _active_size; j++)
            if (_alpha->is_free(j))
                nr_free++;


        if (nr_free * _jml_data > 2 * _active_size * (_jml_data - _active_size))
        {
            for (i = _active_size; i < _jml_data; i++)
            {
                vector<Tmy_double> Q_i = _kernel->get_Q(i, _active_size);
                for (j = 0; j < _active_size; j++)
                    if (_alpha->is_free(j)) {
                        Tmy_double alpha_j = _alpha->get_alpha(j);
                        _arr_G.at(i) = _arr_G.at(i) + (alpha_j * Q_i[j]);
                    }
            }
        }
        else
        {
            for (i = 0; i < _active_size; i++)
                if (_alpha->is_free(i))
                {
                    vector<Tmy_double> Q_i = _kernel->get_Q(i, _jml_data);
                    Tmy_double alpha_i = _alpha->get_alpha(i);
                    for (j = _active_size; j < _jml_data; j++)
                        _arr_G.at(j) = _arr_G.at(j) + (alpha_i * Q_i[j]);
                }
        }

        cout << "end reconstruct gradient" << endl;
    }
}

void Tmy_list_G::set_active_size(int new_value)
{

    _active_size=new_value;
}

void Tmy_list_G::reset_active_size()
{
    _active_size = _jml_data;
}
