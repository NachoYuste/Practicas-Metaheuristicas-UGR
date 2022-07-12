#ifndef SOLUCIONCVOACEC_H_INCLUDED
#define SOLUCIONCVOACEC_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <iostream>
#include <vector>

extern "C" {
    #include "cec17.h"
}

using Random = effolkronium::random_static;
using namespace std;

class CVOA_CEC {
    
    private:
        /*--- Datos ---*/

        arma::mat datos;
        int dim;
        int seed;
        /*--- Parámetros de configuración ---*/

        //Probabilidades
        const int P_DIE = 5;
        const int P_SUPERSPREADER = 10;
        const int P_REINFECTION = 14;
        const int P_ISOLATION = 50;
        const int P_TRAVEL = 10;

        //Ratios de contagio
        const int ORDINARY_RATE = 5;
        const int SUPERSPREADER_RATE = 15;
        const int SOCIAL_DISTANCING_WEEK = 5;

        //Numero de iteraciones
        const int PANDEMIC_DURATION = 20;
        int time;
        int contador;

        //Hibridacion con BLM
        bool memetico;

        /*--- Poblaciones y costes ---*/

        //Poblaciones
        arma::mat infectedPopulation, newInfectedPopulation;
        arma::mat dead, recovered;

        /*--- Funciones del algoritmo ---*/

        //Función de infeccion
        void infect(arma::rowvec infected);
        void newInfection(arma::rowvec infected, bool superspreader, bool traveler);
        arma::mat replicate(arma::rowvec infected, int spread_rate, bool traveler);

        void die();
        int selectBestIndividual(arma::mat infectedPopulation);
        int inPoblation(arma::mat pob, arma::rowvec individual);

        //Funciones para BLM
        arma::rowvec blm(arma::rowvec cromosoma);
        bool Int(arma::mat &  sel, double * costeOriginal, pair <int,int> intercambio, arma::mat datos);
        void mezclaIndex(vector<pair<int,int>> &index);
        

    public:
        //Constructor
        CVOA_CEC(int semilla, int dimension, bool meme);

        //Función principal
        void cvoa_cec();


};

#endif
