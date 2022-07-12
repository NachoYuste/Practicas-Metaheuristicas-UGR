#include "solucionBLM.h"

arma::irowvec blm(int semilla, int m, arma::mat datos){

    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    //Crear sel original y calcular su coste
    arma::mat sel(1,2);
    sel(0) = Random::get(0, (int)(datos.n_rows-1)), 0;

    while(sel.n_rows < m){
        
        int rdn = Random::get(0, (int)(datos.n_rows-1));

        if(arma::all(sel.col(0)!=rdn))
            sel.insert_rows(sel.n_rows, arma::rowvec({(double)rdn, 0}));

    }

    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++)
            if(i!=j)
                sel(i,1) += datos(sel(i,0), sel(j,0));


    double costeOriginal = sel.col(1).max() - sel.col(1).min();
    double * refCosteOriginal = &costeOriginal;

    //Crear índice de espacio de búsqueda
    vector<pair<int,int>> index;

    for(int i=0; i<datos.n_rows;i++)
        for(int j=0; j<datos.n_rows;j++)
            if(i!=j)
                index.push_back(make_pair(i,j));



    bool inSel, outSel;                         //Variables auxiliares para comprobar los índices
    bool seIntercambia = true;                  //Variable de comprobación de intercambio
    int contador = 0;                           //Contador de iteraciones
    vector<pair<int,int>>::iterator it_index;   //Iterador del vecindario 

    
    //Iteramos mientras no recorramos el vecindario al completo sin encontrar una solución mejor
    //o mientras no pasemos de las 100000 iteraciones en total
    
    while(contador < 100000 && seIntercambia){

        mezcla(index);  //"Creamos" un nuevo vecindario.
        seIntercambia = false;

        for(it_index = index.begin(); it_index!=index.end() && !seIntercambia; ++it_index){     //Recorremos el vecindario hasta que encontremos una solución mejor o lo visitimos al completo

            inSel = false;
            outSel = true;

            inSel = arma::any(sel.col(0) == (*it_index).first);
            outSel = arma::all(sel.col(0) != (*it_index).second);
            
            if(inSel && outSel){

                contador++;
                seIntercambia = Int(sel, refCosteOriginal, (*it_index), datos);  //Calculamos si es una solucion aceptable y hacemos el intercambio
                
                if(seIntercambia){
                    index.erase(it_index);      //Borramos el intercambio del espacio de búsqueda
                }
            }
        }        
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(sel.col(0)));    
    return ret;
}

bool Int(arma::mat &  sel, double * costeOriginal, pair <int,int> intercambio, arma::mat datos){
    
    bool seIntercambia = false;
    double distancia_v = 0;
    arma::mat sel_int = sel;


    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<sel.n_rows; i++){
        if(sel_int(i,0) != intercambio.first) distancia_v += datos(intercambio.second, (int)sel_int(i,0));
        sel_int(i,1) -= datos((int)sel_int(i,0), intercambio.first);
        sel_int(i,1) += datos((int)sel_int(i,0), intercambio.second);
    }

    float diferencia_coste = ((distancia_v, sel_int.col(1).max()) - std::min(distancia_v, sel_int.col(1).min())) - *costeOriginal;

    //Condicion de aceptacion de nueva solución
    if(diferencia_coste < 0){
        seIntercambia = true;
        sel = sel_int;  //Actualizamos los costes del sel original
        *costeOriginal += diferencia_coste;  //Actualizamos el coste 
        bool salir = false;

        //Hacemos el intercambio
        for(int i=0; i<sel.n_rows && !salir; i++){
            if(sel(i,0) == intercambio.first){
                sel.row(i) = arma::rowvec({(double)intercambio.second, distancia_v});
                salir = true;
            }
        }    
    }

    return seIntercambia;

} 

void mezcla(vector<pair<int,int>> & index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size());
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }
}