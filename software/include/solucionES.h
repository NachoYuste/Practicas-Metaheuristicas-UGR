#ifndef SOLUCIONES_H_INCLUDED
#define SOLUCIONES_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <math.h>

using Random = effolkronium::random_static;
using namespace std;

arma::irowvec ES(int semilla, int m, arma::mat datos);
bool Int(arma::mat &sel, double delta, double * coste_sel, int interA, int interB, arma::mat datos);

#endif