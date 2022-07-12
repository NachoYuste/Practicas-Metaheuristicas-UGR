#ifndef SOLUCIONBMB_H_INCLUDED
#define SOLUCIONBMB_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <math.h>

using Random = effolkronium::random_static;
using namespace std;

class SolucionBMB {
    
    private:
        arma::mat datos;
        arma::mat mejor_sel;
        int contador;
        int m;

        arma::mat BLocal(arma::mat sel);
        bool Int(arma::mat &sel, double * costeOriginal, pair<int,int> intercambio);
        arma::mat generaSAleatoria();
        void mezcla(vector<pair<int,int>> & index);

    public:
        SolucionBMB(int semilla, int m, arma::mat datos);
        arma::irowvec BMB();
};

#endif