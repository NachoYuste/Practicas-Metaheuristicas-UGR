#include "solucionCVOA_cec17.h"

CVOA_CEC::CVOA_CEC(int semilla, int dimension, bool meme){
    //Almacenar datos
    dim = dimension;
    memetico = meme;

    seed = semilla;
}

void CVOA_CEC::cvoa_cec(){

    for(int funcid = 1; funcid <=30; funcid++){

        newInfectedPopulation.clear();
        infectedPopulation.clear();
        recovered.clear();
        dead.clear();

        arma::rowvec PZ(dim);

        arma::rowvec bestSolution, currentBestSolution;

        double bestCost, currentBestCost;
        int indexCurrentBestSolution;

        cec17_init("CVOA_CEC", funcid, dim);

        int appended = 0;

        std::uniform_real_distribution<> dis(-100.0, 100.0);
        std:mt19937 gen(seed);

        //Patient zero generation
        while(appended < dim){

            int rdn = dis(gen);

            if(arma::all(PZ != rdn)){
                PZ(appended) = rdn;
                appended++;
            }
        }
        
        vector<double> sol(dim);

        infectedPopulation.insert_rows(infectedPopulation.n_rows, PZ);

        bestSolution = PZ;
        sol = arma::conv_to< vector<double> >::from(PZ);
        bestCost = cec17_fitness(&sol[0]);

        time = 0;
        contador = 0;
        
        while (time < PANDEMIC_DURATION && !infectedPopulation.empty() && contador < 10000*dim){
            die();
            newInfectedPopulation.clear();

            for(int i=0; i<infectedPopulation.n_rows; i++)
                infect(infectedPopulation.row(i));

            if (!newInfectedPopulation.empty()){
                indexCurrentBestSolution = selectBestIndividual(newInfectedPopulation);
                currentBestSolution = newInfectedPopulation.row(indexCurrentBestSolution);
                sol = arma::conv_to< vector<double> >::from(currentBestSolution);
                currentBestCost = cec17_fitness(&sol[0]);
                
                if(currentBestCost < bestCost){
                    bestSolution = currentBestSolution;
                    bestCost = currentBestCost;
                }
            }

            //Actualizar lista de individuos recuperados
            recovered.insert_rows(recovered.n_rows, infectedPopulation);

            infectedPopulation = newInfectedPopulation;
            time++;
            cout << "time " << time << endl;
            //cout << "cont " << contador << endl;
        }

        cout <<"Terminada funcion: " << funcid << endl;
    }

    

}

void CVOA_CEC::infect(arma::rowvec infected){

    int rdn1, rdn2;

    rdn1 = Random::get(0, 101);
    rdn2 = Random::get(0, 101);

    if(rdn1 < P_TRAVEL){

        if(rdn2 < P_SUPERSPREADER)
            newInfection(infected, true, true);
        
        else
            newInfection(infected, false, true);

    }

    else{

        if(rdn2 < P_SUPERSPREADER)
            newInfection(infected, true, false);
        
        else
            newInfection(infected, false, false);
    }


}

void CVOA_CEC::newInfection(arma::rowvec infected, bool superspreader, bool traveler){

    arma::mat aux, costAux;

    int spread_rate;

    if(superspreader)   spread_rate = Random::get(0, SUPERSPREADER_RATE);
    else                spread_rate = Random::get(0, ORDINARY_RATE);

    int rdn3 = Random::get(0, 101);
    int rdn4 = Random::get(0, 101);

    aux = replicate(infected, spread_rate, traveler);

    for(int i=0; i < aux.n_rows; i++){

        if((inPoblation(dead, aux.row(i)) == -1) && (inPoblation(newInfectedPopulation, aux.row(i)) == -1)){

            int ind_aux = inPoblation(recovered, aux.row(i));

            if(ind_aux == -1){
                
                if((rdn4 < P_ISOLATION || time < SOCIAL_DISTANCING_WEEK) && (inPoblation(infectedPopulation, aux.row(i)) == -1 ))
                    newInfectedPopulation.insert_rows(newInfectedPopulation.n_rows, aux.row(i));
                
                else
                    recovered.insert_rows(recovered.n_rows, aux.row(i));
            }

            else if(rdn3 < P_REINFECTION){
                newInfectedPopulation.insert_rows(newInfectedPopulation.n_rows, aux.row(i));
                recovered.shed_row(ind_aux);
            }
        }
    }
}

void CVOA_CEC::die(){
    int rdn5;

    for (int i=0; i<infectedPopulation.n_rows; i++){
        rdn5 = Random::get(0, 101);

        if(rdn5 < P_DIE){
            dead.insert_rows(dead.n_rows, infectedPopulation.row(i));
            infectedPopulation.shed_row(i);    
        }
    }
}


arma::mat CVOA_CEC::replicate(arma::rowvec infected, int spread_rate, bool traveler){

    arma::mat replicates;

    //Definicion de variables auxiliares para bucle
    arma::rowvec aux;
    int interA, interB, travel_distance;
    bool outSel;

    if(traveler)    travel_distance = Random::get(1, dim);
    else            travel_distance = 1;


    std::uniform_real_distribution<> dis(-100.0, 100.0);
    std:mt19937 gen(seed);
    
    //Contagia a spread_rate individuos
    for(int i=0; i<spread_rate; i++){
            aux = infected;

            for(int j=0; j<travel_distance; j++){
                
                //Generamos el intercambio aleatorio
                //El primer número debe estar en sel
                outSel = false;
                //El segundo número no puede estar en sel
                interA = Random::get(0, dim-1);
                do{
                    interB = dis(gen);
                    outSel = arma::all(aux != interB);
                }while(!outSel);

                //Intercambiar genes
                aux(interA) = interB;

            }

            //Introducir a individuos replicados
            replicates.insert_rows(replicates.n_rows, aux);            
    }

    return replicates;

}

int CVOA_CEC::selectBestIndividual(arma::mat pob){

    double mejorCoste = -1;
    double actual_cost;
    int bestIndividual;
    vector<double> aux(dim);

    for(int i=0; i<pob.n_rows; i++){
        aux = arma::conv_to< vector<double> >::from(pob.row(i));
        actual_cost = cec17_fitness(&aux[0]);
        contador++;

        if(actual_cost < mejorCoste || i == 0){
            bestIndividual = i;
            mejorCoste = actual_cost;
        }
    }
    return bestIndividual;
}

int CVOA_CEC::inPoblation(arma::mat pob, arma::rowvec individual){

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


/*---------------------------------------*/
/*--------------Métodos BLM--------------*/
/*---------------------------------------*/


arma::rowvec CVOA_CEC::blm(arma::rowvec cromosoma){

    //Crear sel original y calcular su coste
    arma::mat sel = arma::mat(dim, 2);

    for(int i=0; i<dim; i++)
        for (int j=0;j<dim; j++){
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
    int contador = 0;                           //Contador de iteraciones
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

bool CVOA_CEC::Int(arma::mat &  sel, double * costeOriginal, pair <int,int> intercambio, arma::mat datos){
    
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

void CVOA_CEC::mezclaIndex(vector<pair<int,int>> &index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size()-1);
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }
}
