#include "solucionBMB.h"

SolucionBMB::SolucionBMB(int semilla, int m_sel, arma::mat dat){
    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    m = m_sel;
    datos = dat;
}

arma::irowvec SolucionBMB::BMB(){
    
    arma::mat s;
    arma::mat s0 = generaSAleatoria();
    
    mejor_sel = s0;
    double coste_mejor_sel = mejor_sel.col(1).max() - mejor_sel.col(1).min();
    double coste_s;

    for(int i=0; i<10; i++){
        s = BLocal(s0);

        coste_s = s.col(1).max() - s.col(1).min(); 
        if(coste_s < coste_mejor_sel){
            mejor_sel = s;
            coste_mejor_sel = coste_s; 
        }

        s0 = generaSAleatoria();
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(mejor_sel.col(0)));

    return ret;
}


arma::mat SolucionBMB::BLocal(arma::mat sel){

    double costeOriginal = sel.col(1).max() - sel.col(1).min();
    double * refCosteOriginal = &costeOriginal;

    //Crear índice de espacio de búsqueda
    vector<pair<int,int>> index;

    for(int i=0; i<m;i++)
        for(int j=0; j<m;j++)
            if(i!=j)
                index.push_back(make_pair(i,j));



    bool inSel, outSel;                         //Variables auxiliares para comprobar los índices
    bool seIntercambia = true;                  //Variable de comprobación de intercambio
    int contador = 0;                           //Contador de iteraciones
    vector<pair<int,int>>::iterator it_index;   //Iterador del vecindario 
    
    //Iteramos mientras no recorramos el vecindario al completo sin encontrar una solución mejor
    //o mientras no pasemos de las 100000 iteraciones en total
    
    while(contador < 10000 && seIntercambia){

        mezcla(index);  //"Creamos" un nuevo vecindario.
        seIntercambia = false;

        for(it_index = index.begin(); it_index!=index.end() && !seIntercambia; ++it_index){     //Recorremos el vecindario hasta que encontremos una solución mejor o lo visitimos al completo

            inSel = arma::any(sel.col(0) == (*it_index).first);     //Comprobamos si el primer índice está en sel
            outSel = arma::all(sel.col(0) != (*it_index).second);   //Comprobamos si el segundo índice no está en sel
                
            if(inSel && outSel){
                contador++;
                seIntercambia = Int(sel, refCosteOriginal, (*it_index));  //Calculamos si es una solucion aceptable y hacemos el intercambio
                if(seIntercambia){
                    index.erase(it_index);      //Borramos el intercambio del espacio de búsqueda
                }

            }
        }        
    }
    
    return sel;
}

bool SolucionBMB::Int(arma::mat &sel, double * costeOriginal, pair<int,int> intercambio){
    
    bool seIntercambia = false;
    double distancia_v = 0;
    arma::Mat<double> sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<m; i++){ 
        if(sel_int(i,0) != intercambio.first) distancia_v  += datos(intercambio.second, (int)sel_int(i,0));
        sel_int(i,1) -= datos((int)sel_int(i,0), intercambio.first);
        sel_int(i,1) += datos((int)sel_int(i,0), intercambio.second);
    }

    float diferencia_coste = (std::max(distancia_v, sel_int.col(1).max()) - std::min(distancia_v, sel_int.col(1).min())) - *costeOriginal;
   
    //Condicion de aceptacion de nueva solución
    if(diferencia_coste < 0){
        seIntercambia = true;
        sel = sel_int;  //Actualizamos los costes del sel original
        *costeOriginal += diferencia_coste;  //Actualizamos el coste 
        bool salir = false;

        //Hacemos el intercambio
        for(int i=0; i<m && !salir; i++){
            if(sel(i,0) == intercambio.first){
                sel.row(i) = arma::rowvec({(double)intercambio.second, distancia_v});
                salir = true;
            }
        }
    }

    return seIntercambia;
}

arma::mat SolucionBMB::generaSAleatoria(){

    arma::mat sel(1,2);
    sel(0) = Random::get(0, (int)(datos.n_rows-1)), 0;
    while(sel.n_rows < m){
        
        int rdn = Random::get(0, (int)(datos.n_rows-1));

        if(arma::all(sel.col(0)!=rdn))
            sel.insert_rows(sel.n_rows, arma::rowvec({(double)rdn, 0}));

    }

    double costeOriginal = 0;
    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++)
            if(i!=j){
                sel(i,1) += datos(sel(i,0), sel(j,0));
            }    

    costeOriginal = sel.col(1).max() - sel.col(1).min();

    return sel;
}

void SolucionBMB::mezcla(vector<pair<int,int>> & index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size());
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }

}