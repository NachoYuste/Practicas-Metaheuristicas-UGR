#ifndef SOLUCIONMEMETICA_H_INCLUDED
#define SOLUCIONMEMETICA_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <vector>
#include <utility>

using Random = effolkronium::random_static;
using namespace std;

enum class tipo_ejecucion{AM1, AM2, AM3};

class Memetico{

    private:
        arma::mat datos;
        arma::mat poblacion;
        arma::mat sel;
        arma::rowvec costePoblacion;

        tipo_ejecucion tipoEjecucion;

        const int numCromosomas = 10;
        int m;
        int t;
        double costeOriginal;

        arma::rowvec cruceUniforme(arma::rowvec padre1, arma::rowvec padre2);
        double calculaCoste(arma::rowvec cromosoma);
        double calculaDistanciaTotal(arma::rowvec cromosoma);
        
        //Métodos para la búsqueda local
        arma::irowvec blm(arma::rowvec cromosoma);
        bool Int(pair <int,int> intercambio);
        void mezclaIndex(vector<pair<int,int>> &index);

    public:

        Memetico(int semilla, int m, arma::mat datos, tipo_ejecucion tipo);

        arma::irowvec generacional();
        

};


#endif