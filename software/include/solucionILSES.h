#ifndef SOLUCIONILSES_H_INCLUDED
#define SOLUCIONILSES_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <math.h>

using Random = effolkronium::random_static;
using namespace std;

class SolucionILSES {
    
    private:
        arma::mat datos;
        int contador;
        int m;
        double T;
        const double Tf = pow(10, -3);

        arma::mat BES(arma::mat sel);
        bool Int(arma::mat &sel, double delta, double & coste_sel, int interA, int interB, arma::mat datos);
        arma::mat generaSelAleatoria();
        arma::mat mutacion(arma::mat sel);

    public:
        SolucionILSES(int semilla, int m, arma::mat datos);
        arma::irowvec ILSES();
};

#endif