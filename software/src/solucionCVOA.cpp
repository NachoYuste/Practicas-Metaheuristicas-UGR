#include "solucionCVOA.h"

CVOA::CVOA(int semilla, int seleccion, arma::mat dat, tipo_memetico meme){
    //Almacenar datos
    datos = dat;
    m = seleccion;
    memetico = meme;

    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);
}

arma::irowvec CVOA::cvoa(){

    int dim = 10;

    arma::rowvec PZ(m);
    arma::rowvec costPZ(m, arma::fill::zeros);

    arma::rowvec bestSolution, currentBestSolution;

    double bestCost, currentBestCost;
    int indexCurrentBestSolution;

    int appended = 0;

    //Patient zero generation
    while(appended < m){

        int rdn = Random::get(0, (int)(datos.n_rows-1));

        if(arma::all(PZ != rdn)){
            PZ(appended) = rdn;
            appended++;
        }
    }

    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++)
            if(i!=j)
                costPZ(i) += datos(PZ(i), PZ(j));
    
    int maxinfected = 0, maxrecovered = 0, maxdead = 0;

    
    infectedPopulation.insert_rows(infectedPopulation.n_rows, PZ);
    costInfectedPopulation.insert_rows(costInfectedPopulation.n_rows, costPZ);

    bestSolution = PZ;
    bestCost = costPZ.max() - costPZ.min();

    time = 0;
    contador = 0;

    while (time < PANDEMIC_DURATION && !infectedPopulation.empty() && contador < 100000){
        die();

        newInfectedPopulation.clear();
        costNewInfectedPopulation.clear();

        for(int i=0; i<infectedPopulation.n_rows; i++)
            infect(infectedPopulation.row(i), costInfectedPopulation.row(i));
            
        if (!newInfectedPopulation.empty()){
            
            if(memetico == tipo_memetico::ES){
                indexCurrentBestSolution = selectBestIndividual(costNewInfectedPopulation);
                currentBestSolution = BES(newInfectedPopulation.row(indexCurrentBestSolution));
                currentBestCost = cost(currentBestSolution);
            }
            else if (memetico == tipo_memetico::BLM){
                indexCurrentBestSolution = selectBestIndividual(costNewInfectedPopulation);
                currentBestSolution = blm(newInfectedPopulation.row(indexCurrentBestSolution));
                currentBestCost = cost(currentBestSolution);
            }
            else if (memetico == tipo_memetico::NONE){
                indexCurrentBestSolution = selectBestIndividual(costNewInfectedPopulation);
                currentBestCost = costNewInfectedPopulation.row(indexCurrentBestSolution).max() - costNewInfectedPopulation.row(indexCurrentBestSolution).min();
                currentBestSolution = newInfectedPopulation.row(indexCurrentBestSolution);
            }

            if(currentBestCost < bestCost){
                bestSolution = currentBestSolution;
                bestCost = currentBestCost;
            }
        }

        //Actualizar lista de individuos recuperados
        recovered.insert_rows(recovered.n_rows, infectedPopulation);

        infectedPopulation = newInfectedPopulation;
        costInfectedPopulation = costNewInfectedPopulation;
        
        time++;
        //cout << "Iter: " << time << endl;
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(bestSolution));  
    return ret;
}

void CVOA::infect(arma::rowvec infected, arma::rowvec costInfected){

    int rdn1, rdn2;

    rdn1 = Random::get(0, 101);
    rdn2 = Random::get(0, 101);

    if(rdn1 < P_TRAVEL){

        if(rdn2 < P_SUPERSPREADER)
            newInfection(infected, costInfected, true, true);
        
        else
            newInfection(infected, costInfected,false, true);

    }

    else{

        if(rdn2 < P_SUPERSPREADER)
            newInfection(infected, costInfected, true, false);
        
        else
            newInfection(infected, costInfected, false, false);
    }


}

void CVOA::newInfection(arma::rowvec infected, arma::rowvec costInfected, bool superspreader, bool traveler){

    arma::mat aux, costAux;

    int spread_rate;

    if(superspreader)   spread_rate = Random::get(0, SUPERSPREADER_RATE);
    else                spread_rate = Random::get(0, ORDINARY_RATE); 

    int rdn3 = Random::get(0, 101);
    int rdn4 = Random::get(0, 101);

    aux = replicate(infected, costInfected, costAux, spread_rate, traveler);

    for(int i=0; i < aux.n_rows; i++){

        if((inPoblation(dead, aux.row(i)) == -1) && (inPoblation(newInfectedPopulation, aux.row(i)) == -1)){

            int ind_aux = inPoblation(recovered, aux.row(i));

            if(ind_aux == -1){
                
                if((rdn4 < P_ISOLATION || time < SOCIAL_DISTANCING_WEEK) && (inPoblation(infectedPopulation, aux.row(i)) == -1 )){
                    newInfectedPopulation.insert_rows(newInfectedPopulation.n_rows, aux.row(i));
                    costNewInfectedPopulation.insert_rows(costNewInfectedPopulation.n_rows, costAux.row(i));
                }
                
                else
                    recovered.insert_rows(recovered.n_rows, aux.row(i));
            }

            else if(rdn3 < P_REINFECTION){
                newInfectedPopulation.insert_rows(newInfectedPopulation.n_rows, aux.row(i));
                costNewInfectedPopulation.insert_rows(costNewInfectedPopulation.n_rows, costAux.row(i));
                recovered.shed_row(ind_aux);
            }
        }
    }


}

void CVOA::die(){
    int rdn5;

    for (int i=0; i<infectedPopulation.n_rows; i++){
        rdn5 = Random::get(0, 101);

        if(rdn5 <= P_DIE){
            dead.insert_rows(dead.n_rows, infectedPopulation.row(i));
            infectedPopulation.shed_row(i);    
        }
    }
}


arma::mat CVOA::replicate(arma::rowvec infected, arma::rowvec costInfected, arma::mat & costRet, int spread_rate, bool traveler){

    arma::mat replicates;

    //Definicion de variables auxiliares para bucle
    arma::rowvec aux, costAux;
    int interA, interB, travel_distance;
    bool outSel;
    double distancia_v = 0;

    if(traveler)    travel_distance = Random::get(1, (int)infected.n_cols);
    else            travel_distance = 1;

    //Contagia a spread_rate individuos
    for(int i=0; i<spread_rate; i++){
            aux = infected;
            costAux = costInfected;

            for(int j=0; j<travel_distance; j++){
                
                //Generamos el intercambio aleatorio
                //El primer número debe estar en sel
                outSel = false;
                //El segundo número no puede estar en sel
                interA = Random::get(0, (int)(aux.n_cols-1));
                do{
                    interB = Random::get(0, (int)(datos.n_rows-1));
                    outSel = arma::all(aux != interB);
                }while(!outSel);


                //Calculo el nuevo coste de manera factorizado
                for(int k=0; k<costAux.n_cols; k++){
                    if(k != interA)    distancia_v += datos(interB, aux(k));
                    costAux(k) -= datos(aux(k), aux(interA));
                    costAux(k) += datos(aux(k), interB);
                }

                //Intercambiar genes
                aux(interA) = interB;
                costAux(interA) = distancia_v;

            }

            //Introducir a individuos replicados
            replicates.insert_rows(replicates.n_rows, aux);
            costRet.insert_rows(costRet.n_rows, costAux);
            
    }

    return replicates;

}

int CVOA::selectBestIndividual(arma::mat pob){

    double mejorCoste = INT_MAX;
    double actual_cost;
    int bestIndividual;
    vector<double> aux(10);

    for(int i=0; i<pob.n_rows; i++){
        actual_cost = pob.row(i).max() - pob.row(i).min();
        contador++;

        if(actual_cost < mejorCoste || i ==1){
            bestIndividual = i;
            mejorCoste = actual_cost;
        }
    }
    return bestIndividual;
}

int CVOA::inPoblation(arma::mat pob, arma::rowvec individual){

    int ind_individual = -1;
    bool inPob = false;

    for(int i=0; i<pob.n_rows && !inPob; i++){
        if(arma::approx_equal(pob.row(i), individual, "absdiff", 0.9)){
            ind_individual = i;
            inPob = true;
        }
    }

    return ind_individual;

}

double CVOA::cost(arma::rowvec row){
    arma::rowvec ret(row.n_cols, arma::fill::zeros);
    for(int i=0; i< row.n_cols; i++)
        for(int j=0; j<row.n_cols; j++)
            ret(i)+=datos(row(i), row(j));


    return ret.max() - ret.min();
        
}


/*---------------------------------------*/
/*--------------Métodos BLM--------------*/
/*---------------------------------------*/


arma::rowvec CVOA::blm(arma::rowvec cromosoma){

    //Crear sel original y calcular su coste
    arma::mat sel = arma::mat(m, 2);

    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++){
            sel(i,0) = cromosoma(i);
            sel(i,1) += datos(sel(i,0), sel(j,0));
        }
    
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
    vector<pair<int,int>>::iterator it_index;   //Iterador del vecindario 

    
    //Iteramos mientras no recorramos el vecindario al completo sin encontrar una solución mejor
    //o mientras no pasemos de las 100000 iteraciones en total
    
    while(contador < 400 && seIntercambia){

        mezclaIndex(index);  //"Creamos" un nuevo vecindario.
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
    return sel.col(0).as_row();
}

bool CVOA::Int(arma::mat &  sel, double * costeOriginal, pair <int,int> intercambio, arma::mat datos){
    
    bool seIntercambia = false;
    double distancia_v = 0;
    arma::Mat<double> sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<sel.n_rows; i++){ 
        if(sel_int(i,0) != intercambio.first) distancia_v += datos(intercambio.second, (int)sel_int(i,0));
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
        for(int i=0; i<sel.n_rows && !salir; i++){
            if(sel(i,0) == intercambio.first){
                sel.row(i) = arma::rowvec({(double)intercambio.second, distancia_v});
                salir = true;
            }
        }    
    }
    return seIntercambia;

} 

void CVOA::mezclaIndex(vector<pair<int,int>> &index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size()-1);
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }
}


/*---------------------------------------*/
/*--------------Métodos ES--------------*/
/*---------------------------------------*/


arma::rowvec CVOA::BES(arma::rowvec cromosoma){
       //Establecemos los valores iniciales
    int max_vecinos = 10 * datos.n_rows;
    int max_exitos = 0.1 * max_vecinos;

    //Crear sel original y calcular su coste
    arma::mat sel = arma::mat(m, 2);

    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++){
            sel(i,0) = cromosoma(i);
            sel(i,1) += datos(sel(i,0), sel(j,0));
        }

    //Calculo la temperatura de sel
    double T = (0.3 * (sel.col(1).max() - sel.col(1).min())) / (-log(0.3));  
    double Tf = Tf = pow(10, -3);

    //Calculo valor beta para las actualizanes de temperatura
    double beta = (T - Tf) / ((100000/max_vecinos)*T*Tf);

    //Almacenamos la "mejor solucion"
    arma::mat mejor_solucion = sel;
    double coste_mejor_sel = mejor_solucion.col(1).max() - mejor_solucion.col(1).min();

    int numEnfriamientos = 0; //contador de enfriamientos
    arma::rowvec sel_prima; //Sel tras el intercambio
    
    //Definición de variables auxiliares para los bucles
    int interA, interB;     //Indices para el intercambio
    bool inSel, outSel;     //Booleanos auxiliares para el intercambio
    bool seIntercambia = false, seHaceUnIntercambio = true;
    int exitos = 0;
    double coste_sel = coste_mejor_sel;

    do{
        
        seHaceUnIntercambio = false;
        exitos = 0;
        
        for(int i=0; i<max_vecinos && exitos < max_exitos; i++){
            
            //Aumentamos numero de iteraciones totales
            contador++;
            
            //Generamos el intercambio aleatorio
            //El primer número debe estar en sel
            //El segundo número no puede estar en sel
            do{
                interA = Random::get(0, (int)(datos.n_rows-1));
                interB = Random::get(0, (int)(datos.n_rows-1));
                inSel = arma::any(sel.col(0) == interA);
                outSel = arma::all(sel.col(0) != interB);
                
            }while(!inSel || !outSel || interA==interB);
            

            //Hacemos el intercambio (si procede)
            seIntercambia = IntES(sel, T * numEnfriamientos, coste_sel, interA, interB, datos);

            //Comprobamos si se hace aunque sea un intercambio en todo el bucle interno
            seHaceUnIntercambio = seHaceUnIntercambio || seIntercambia;

            //Si se produce el intercambio actualizamos max_exitos
            if(seIntercambia){
                exitos++;

                //Si se produce el intercambio y la nueva solución es mejor se actualiza la mejor solución
                if(coste_sel < coste_mejor_sel){
                    mejor_solucion = sel;
                    coste_mejor_sel = coste_sel;
                }
            }
            
        }

        //Actualizamos la temperatura
        numEnfriamientos++;
        T = T / (1 + (beta*T));

    }while(Tf < T && contador < 400 && exitos!=0);


    //Devolvemos la mejor solución
    return mejor_solucion.col(0).as_row();

}


bool CVOA::IntES(arma::mat &sel, double delta, double & coste_sel, int interA, int interB, arma::mat datos){
    
   bool seIntercambia = false;
    double distancia_v = 0;
    arma::mat sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<sel.n_rows; i++){ 
        if(sel_int(i,0) != interA) distancia_v  += datos(interB, (int)sel_int(i,0));
        sel_int(i,1) -= datos((int)sel_int(i,0), interA);
        sel_int(i,1) += datos((int)sel_int(i,0), interB);
    }

    float diferencia_coste = (std::max(distancia_v, sel_int.col(1).max()) - std::min(distancia_v, sel_int.col(1).min())) - coste_sel;

    bool pasaProbAceptacion = Random::get(0.0, 1.0) <= exp(((-1)*diferencia_coste)/delta);;

    //Condicion de aceptacion de nueva solución
    if(diferencia_coste < 0 || pasaProbAceptacion){
        seIntercambia = true;
        sel = sel_int;  //Actualizamos los costes del sel original
        coste_sel += diferencia_coste;
        contador++;
        bool salir = false;

        //Hacemos el intercambio
        for(int i=0; i<sel.n_rows && !salir; i++){
            if(sel(i,0) == interA){
                sel.row(i) = arma::rowvec({(double)interB, distancia_v});
                salir = true;
            }
        }         
    }
    
    return seIntercambia;
}
