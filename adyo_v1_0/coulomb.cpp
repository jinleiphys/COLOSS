#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath> 

using namespace std;

#define SIGN(a) (((a) < 0) ? (-1) : (1))

const double precision = 1E-10, sqrt_precision = 1E-5;

#include "complex_functions.H"
#include "cwfcomp.cpp"
#include "test_rec_rel.cpp"

// extern function for fortran
extern "C" {
    void compute_coulomb_wave_functions(
        const char* is_normalized,
        complex<double>* eta,
        complex<double>* z,
        complex<double>* l,
        int* Nl,
        complex<double>* F,
        complex<double>* DF,
        complex<double>* G,
        complex<double>* DG
    ) {
        bool is_it_normalized = (string(is_normalized) == "true");
        complex<double> eta_val = *eta;
        complex<double> l_val = *l;
        int Nl_val = *Nl;
        complex<double> z_val = *z;

        Coulomb_wave_functions cwf(is_it_normalized, l_val, eta_val);

        complex<double> *F_tab = new complex<double>[Nl_val];
        complex<double> *dF_tab = new complex<double>[Nl_val];
        complex<double> *G_tab = new complex<double>[Nl_val];
        complex<double> *dG_tab = new complex<double>[Nl_val];
        complex<double> *Hp_tab = new complex<double>[Nl_val];
        complex<double> *dHp_tab = new complex<double>[Nl_val];
        complex<double> *Hm_tab = new complex<double>[Nl_val];
        complex<double> *dHm_tab = new complex<double>[Nl_val];

        cwf_l_tables_recurrence_relations(l_val, Nl_val, eta_val, is_it_normalized, z_val, F_tab, dF_tab, G_tab, dG_tab, Hp_tab, dHp_tab, Hm_tab, dHm_tab);

        for (int i = 0; i < Nl_val; i++) {
            F[i] = F_tab[i];
            DF[i] = dF_tab[i];
            G[i] = G_tab[i];
            DG[i] = dG_tab[i];
        }

        delete[] F_tab;
        delete[] dF_tab;
        delete[] G_tab;
        delete[] dG_tab;
        delete[] Hp_tab;
        delete[] dHp_tab;
        delete[] Hm_tab;
        delete[] dHm_tab;
    }
}
