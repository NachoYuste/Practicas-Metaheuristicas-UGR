#include "solucionILS.h"

SolucionILS::SolucionILS(int semilla, int m_sel, arma::mat dat){
    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    m = m_sel;
    datos = dat;
}

arma::irowvec SolucionILS::ILS(){
    
    //Generación de solución inicial
    arma::mat s0(generaSelAleatoria());
    arma::mat mejor_sel(BLocal(s0));

    //Declaración de variables para el bucle
    arma::mat sprima, sprima2;
    double coste_mejor_sel = mejor_sel.col(1).max() - mejor_sel.col(1).min();
    double coste_sprima2;

    for(int i=0; i<10; i++){
        //Operador de mutación
        sprima = mutacion(mejor_sel);

        //Aplico la búsqueda local al sel mutado
        sprima2 = BLocal(sprima);

        //Coste de la nueva solución
        coste_sprima2 = sprima2.col(1).max() - sprima2.col(1).min();

        //Si la solución actual es mejor que la mejor de todas, la actualizo
        if(coste_sprima2 < coste_mejor_sel){

            mejor_sel = sprima2;
            coste_mejor_sel = coste_sprima2;
        }
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(mejor_sel.col(0)));
    return ret;
}


arma::mat SolucionILS::BLocal(arma::mat sel){

    double costeOriginal = sel.col(1).max() - sel.col(1).min();

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
    double * refCoste = &costeOriginal;

    while(contador < 10000 && seIntercambia){

        mezcla(index);  //"Creamos" un nuevo vecindario.
        seIntercambia = false;

        for(it_index = index.begin(); it_index!=index.end() && !seIntercambia; ++it_index){     //Recorremos el vecindario hasta que encontremos una solución mejor o lo visitimos al completo

            inSel = arma::any(sel.col(0) == (*it_index).first);     //Comprobamos si el primer índice está en sel
            outSel = arma::all(sel.col(0) != (*it_index).second);   //Comprobamos si el segundo índice no está en sel
            
            if(inSel && outSel){
                contador++;
                seIntercambia = Int(sel, refCoste, (*it_index));  //Calculamos si es una solucion aceptable y hacemos el intercambio
                
                if(seIntercambia){
                    index.erase(it_index);      //Borramos el intercambio del espacio de búsqueda
                }

            }
        }        
    }
    
    return sel;
}

bool SolucionILS::Int(arma::mat &sel, double * costeOriginal, pair<int,int> intercambio){
    
    bool seIntercambia = false;
    double distancia_v = 0;
    arma::mat sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<m; i++){ 
        if(sel_int(i,0) != intercambio.first)    distancia_v  += datos(intercambio.second, (int)sel_int(i,0));
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

arma::mat SolucionILS::generaSelAleatoria(){

    arma::mat sel(1,2);
    sel(0) = Random::get(0, (int)(datos.n_rows-1)), 0;

    while(sel.n_rows < m){
        
        int rdn = Random::get(0, (int)(datos.n_rows-1));

        if(arma::all(sel.col(0)!=rdn))
            sel.insert_rows(sel.n_rows, arma::rowvec({(double)rdn, 0}));

    }

    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++)
                sel(i,1) += datos(sel(i,0), sel(j,0));
    
    return sel;
}

arma::mat SolucionILS::mutacion(arma::mat sel){
    int numMutaciones = 0.3 * m;
    int interA, interB;
    bool inSel, outSel, salir;
    double distancia_v;
    arma::mat sel_int;

    for(int i=0; i< numMutaciones; i++){

        //Generamos el intercambio aleatorio
        //El primer número debe estar en sel
        //El segundo número no puede estar en sel
        do{
            interA = Random::get(0, (int)(datos.n_rows-1));
            interB = Random::get(0, (int)(datos.n_rows-1));
            inSel = arma::any(sel.col(0) == interA);
            outSel = arma::all(sel.col(0) != interB);
            
        }while(!inSel || !outSel || interA==interB);


        //Operador de intercambio
        distancia_v = 0;
        sel_int = sel;

        //Cálculo factorizado de coste de intercambio
        for(int i=0; i<m; i++){ 
            if(sel_int(i,0) != interA)    distancia_v  += datos(interB, (int)sel_int(i,0));
            sel_int(i,1) -= datos((int)sel_int(i,0), interA);
            sel_int(i,1) += datos((int)sel_int(i,0), interB);
        }

        //Hacemos el intercambio
        salir = false;

        sel.col(1) = sel_int.col(1);

        for(int i=0; i<m && !salir; i++){
            if(sel(i,0) == interA){
                sel.row(i) = arma::rowvec({(double)interB, distancia_v});
                salir = true;
            }
        }
    }

    return sel;
}

void SolucionILS::mezcla(vector<pair<int,int>> & index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size());
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }

}