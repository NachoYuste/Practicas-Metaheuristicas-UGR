#ifndef SOLUCIONBLM_H_INCLUDED
#define SOLUCIONBLM_H_INCLUDED

#include <vector>
#include <utility>
#include <armadillo>
#include "random.hpp"

using Random = effolkronium::random_static;
using namespace std;


arma::irowvec blm(int semilla, int m, arma::Mat<double> datos);

/*
    Comprueba si el conjunto resultante de intercambiar es mejor
    que el conjunto original
    Si lo es, devuelve true, si no false.
*/
bool Int(arma::Mat<double>&  sel,double * costeOriginal, pair <int,int> intercambio, arma::Mat<double> datos);

/*
    Mezcla aleatoriamente el vector de pares que conforma el vecindario
*/
void mezcla(vector<pair<int,int>> & index);
#endif