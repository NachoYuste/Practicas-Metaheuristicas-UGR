#ifndef SOLUCIONGENETICA_H_INCLUDED
#define SOLUCIONGENETICA_H_INCLUDED

#include "random.hpp"
#include <armadillo>


using Random = effolkronium::random_static;
using namespace std;

enum class tipo_cruce{UNIFORME, POSICION};
enum class tipo_genetico{ESTACIONARIO, UNIFORME};

class Genetico{

    private:
        arma::mat datos;
        arma::mat poblacion;
        arma::rowvec costePoblacion;

        const int numCromosomas = 50;
        int m;

        tipo_cruce cruce;
        tipo_genetico algoritmo;

        arma::rowvec cruceUniforme(arma::rowvec padre1, arma::rowvec padre2);
        arma::mat crucePosicion(arma::rowvec padre1, arma::rowvec padre2);
        double calculaCoste(arma::rowvec cromosoma);
        double calculaDistanciaTotal(arma::rowvec cromosoma);

    public:

        Genetico(int semilla, int m, arma::mat datos, tipo_cruce cruce);

        arma::irowvec generacional();
        arma::irowvec estacionario();
        

};



#endif