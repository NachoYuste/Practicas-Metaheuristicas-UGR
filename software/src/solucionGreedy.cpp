#include "solucionGreedy.h"
#include <iostream>


arma::irowvec greedy(int semilla, int m, arma::Mat<double> datos){
    
    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    //Inicialización de variables y conjuntos de datos
    arma::mat solucion;

    arma::mat CL (datos.n_rows, 2);  //Almacena los puntos y su valor heurístico
    CL.col(0) = arma::regspace(0,datos.n_rows-1);    //Inicializar con todos los puntos
    
    int filaRandom = Random::get(0, (int)(datos.n_rows-1));
    solucion.insert_rows(solucion.n_rows, arma::rowvec({(double)filaRandom, 0}));
    
    //Eliminamos la fila seleccionado (quitamos el punto inicial del espacio de búsqueda)
    CL.shed_row(filaRandom);

    pair<int,int> u_index;

    //Bucle principal
    while(solucion.n_rows < m){
        //Obtenemos el mejor punto a añadir y lo guardamos
        u_index = funcionHeuristica(solucion, CL, datos);
        solucion.insert_rows(solucion.n_rows, arma::rowvec({(double)u_index.first, 0}));

        //Eliminamos el punto utilizado del espacio de búsqueda
        CL.shed_row(u_index.second);
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(solucion.col(0)));

    return ret;
}


pair<int,int> funcionHeuristica(arma::Mat<double> sel, arma::Mat<double> rcl, arma::Mat<double> datos){

    //Cálculo del valor heurístico de RCL
    for (int u=0; u<rcl.n_rows; u++){
        for(int v=0; v<sel.n_rows; v++)      
            rcl(u,1)+=datos((int)rcl(u,0), (int)sel(v,0));
         
        for(int v=0; v<sel.n_rows; v++)
            sel(v,1) = rcl(u,1) + datos((int)rcl(u,0), (int)sel(v,0));


        rcl(u,1) = std::max(rcl(u,1), sel.col(1).max()) - std::min(rcl(u,1), sel.col(1).min());
    }


    //Coeficiente con menor valor heurístico
    int index= rcl.col(1).index_min();

    return make_pair(rcl.col(0)(index), index);
}