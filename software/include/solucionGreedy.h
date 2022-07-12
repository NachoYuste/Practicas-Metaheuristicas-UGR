#ifndef SOLUCOINGREEDY_H_INCLUDED
#define SOLUCOINGREEDY_H_INCLUDED


#include <armadillo>
#include <vector>
#include <algorithm>
#include "random.hpp"

using Random = effolkronium::random_static;
using namespace std;

arma::irowvec greedy(int semilla, int m, arma::Mat<double> datos);

pair<int,int> funcionHeuristica(arma::Mat<double> sel, arma::Mat<double> rcl, arma::Mat<double> datos);


#endif