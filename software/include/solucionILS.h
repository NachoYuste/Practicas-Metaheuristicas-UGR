#ifndef SOLUCIONILS_H_INCLUDED
#define SOLUCIONILS_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <math.h>

using Random = effolkronium::random_static;
using namespace std;

class SolucionILS {
    
    private:
        arma::mat datos;
        int contador;
        int m;

        arma::mat BLocal(arma::mat sel);
        bool Int(arma::mat &sel, double * costeOriginal, pair<int,int> intercambio);
        arma::mat generaSelAleatoria();
        arma::mat mutacion(arma::mat sel);
        void mezcla(vector<pair<int,int>> & index);


    public:
        SolucionILS(int semilla, int m, arma::mat datos);
        arma::irowvec ILS();
};

#endif