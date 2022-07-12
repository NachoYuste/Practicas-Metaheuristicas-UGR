#ifndef SOLUCIONCVOA_H_INCLUDED
#define SOLUCIONCVOA_H_INCLUDED

#include "random.hpp"
#include <armadillo>
#include <iostream>

using Random = effolkronium::random_static;
using namespace std;

enum class tipo_memetico {BLM, ES, NONE};

class CVOA {
    
    private:
        /*--- Datos ---*/

        arma::mat datos;
        int m;

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
        const int SOCIAL_DISTANCING_WEEK = 3;

        //Numero de iteraciones
        const int PANDEMIC_DURATION = 15;
        int time;
        int contador;

        //Hibridacion con BLM
        tipo_memetico memetico;

        /*--- Poblaciones y costes ---*/

        //Poblaciones
        arma::mat infectedPopulation, newInfectedPopulation;
        arma::mat dead, recovered;

        //Costes
        arma::mat costInfectedPopulation, costNewInfectedPopulation;
        arma::rowvec costInfectedIndividuals, costNewInfectedIndividuals;


        /*--- Funciones del algoritmo ---*/

        //Función de infeccion
        void infect(arma::rowvec infected, arma::rowvec costInfected);
        void newInfection(arma::rowvec infected, arma::rowvec costInfected, bool superspreader, bool traveler);
        arma::mat replicate(arma::rowvec infected, arma::rowvec costInfected, arma::mat & costAux, int spread_rate, bool traveler);

        void die();
        int selectBestIndividual(arma::mat costInfectedPopulation);
        int inPoblation(arma::mat pob, arma::rowvec individual);

        double cost(arma::rowvec row);

        //Funciones para BLM
        arma::rowvec blm(arma::rowvec cromosoma);
        bool Int(arma::mat &  sel, double * costeOriginal, pair <int,int> intercambio, arma::mat datos);
        void mezclaIndex(vector<pair<int,int>> &index);

        //Funciones para ES
        arma::rowvec BES(arma::rowvec cromosoma);
        bool IntES(arma::mat &sel, double delta, double & coste_sel, int interA, int interB, arma::mat datos);

        

    public:
        //Constructor
        CVOA(int semilla, int m, arma::mat datos, tipo_memetico memetico);

        //Función principal
        arma::irowvec cvoa();


};

#endif
