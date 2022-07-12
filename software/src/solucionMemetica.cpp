#include "solucionMemetica.h"


Memetico::Memetico(int semilla, int seleccion, arma::mat dat, tipo_ejecucion tipo){
    //Almacenar datos y fijar tipo de algoritmo y cruce
    datos = dat;
    m = seleccion;
    t = 0;
    tipoEjecucion = tipo;

    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    //Crear población inicial aleatoria y calcular su coste
    poblacion = arma::mat(numCromosomas, datos.n_rows);
    costePoblacion = arma::rowvec(numCromosomas);
    
    for(int i=0; i<numCromosomas; i++){

        while(arma::sum(poblacion.row(i)) < m){
            int genRandom = Random::get(0, (int)(datos.n_rows-1));
            poblacion(i, genRandom) = 1;
        }

        costePoblacion(i) = calculaCoste(poblacion.row(i));
    }
}



/*---------------------------------------*/
/*--------------Métodos AGG--------------*/
/*---------------------------------------*/

arma::irowvec Memetico::generacional(){

    //Variables del algoritmo
    int numCruces = 0.7 * numCromosomas / 2;
    int numMutaciones = 0.1 * (numCromosomas * datos.n_rows);
    arma::mat pob_prima(numCromosomas, datos.n_rows);
    arma::rowvec mejorSolucion = poblacion.row(costePoblacion.index_min());
    double mejorCoste = costePoblacion.min();

    //Definición de variables auxiliares de los bucles
    int rdn1, rdn2, rdn3;
    arma::rowvec hijo1, hijo2, auxBLM;
    arma::uvec genes1, genes0;
    arma::irowvec genesBLM;
    arma::mat hijos;
    int generaciones = 0;

    while (t<100000){

        //Operación de selección
        for(int i=0; i<numCromosomas; i++){
            rdn1 = Random::get(0, numCromosomas-1);
            rdn2 = Random::get(0, numCromosomas-1);
            
            if(costePoblacion(rdn1) < costePoblacion(rdn2)){
                pob_prima.row(i) = poblacion.row(rdn1);
            }
            else{
                pob_prima.row(i) = poblacion.row(rdn2);
            }
        }
        
        
        //Operación de cruce uniforme
        for (int i=0; i<numCruces-1; i+=2){
            //Obtenemos los dos hijos de dos padres
            hijo1 = cruceUniforme(pob_prima.row(i), pob_prima.row(i+1));
            hijo2 = cruceUniforme(pob_prima.row(i), pob_prima.row(i+1));
            
            //Actualizamos la poblacion
            pob_prima.row(i) = hijo1;
            pob_prima.row(i+1) = hijo2;
        }

        //Operación de mutación
        for(int i=0; i<numMutaciones; i++){
            
            //Cromosoma aleatorio
            rdn1 = Random::get(0, (int)(numCromosomas-1));

            //Genes aleatorios a intercambiar
            genes1 = arma::find(pob_prima.row(rdn1) == 1);
            genes0 = arma::find(pob_prima.row(rdn1) == 0);

            rdn2 = Random::get(0, (int)(genes1.size()-1));
            rdn3 = Random::get(0, (int)(genes0.size()-1));

            //Intercambiar genes
            pob_prima.swap_cols(genes1(rdn2), genes0(rdn3));
        }
        
        //OPTIMIZACION CON BÚSQUEDA LOCAL

        if(generaciones%10 == 0){

            //AM-(10,1.0) - AM1
            if(tipoEjecucion == tipo_ejecucion::AM1){
                for(int i=0; i<pob_prima.n_rows; i++){
                    genesBLM = blm(pob_prima.row(i));
                    pob_prima.row(i) = arma::rowvec(pob_prima.n_cols, arma::fill::zeros);
                    for(int j=0; j<genesBLM.size(); j++)
                        pob_prima(i, genesBLM(j)) = 1;
                    
                }
            }

            //AM -(10, 0.1) - AM2
            if(tipoEjecucion == tipo_ejecucion::AM2){
                for(int i=0; i<pob_prima.n_rows; i++){
                    rdn1 = Random::get(0,100);
                    
                    if(rdn1 <=10){
                        genesBLM = blm(pob_prima.row(i));
                        pob_prima.row(i) = arma::rowvec(pob_prima.n_cols, arma::fill::zeros);
                        for(int j=0; j<genesBLM.size(); j++)
                            pob_prima(i, genesBLM(j)) = 1;
                    }
                }
            }

            //AM - (10, 0,1mej) -AM3
            if(tipoEjecucion == tipo_ejecucion::AM3){
                //Como nuestra poblacion solo tiene 10 cromosomas sólo aplicamos BLM en uno de ellos
                genesBLM= blm(pob_prima.row(costePoblacion.index_min()));

                pob_prima.row(costePoblacion.index_min()) = arma::rowvec(pob_prima.n_cols, arma::fill::zeros);
                for(int j=0; j<genesBLM.size(); j++)
                    pob_prima(costePoblacion.index_min(), genesBLM(j)) = 1;
                
            }
        }
        
        
        //Evaluar las soluciones
        for(int i=0; i<costePoblacion.n_cols; i++){
            costePoblacion(i) = calculaCoste(pob_prima.row(i));
            t++;
        }


        //Reemplazar poblacion
        poblacion = pob_prima;
        
        //Comprobar si hay que actualizar la mejor solución
        if(costePoblacion.max() > mejorCoste){
            //Sustituyo la peor solución por la mejor de la anterior generación
            poblacion.row(costePoblacion.index_max()) = mejorSolucion;
            costePoblacion(costePoblacion.index_max()) = mejorCoste;

            //Actualizo la mejor solución con la nueva generación
            mejorCoste = costePoblacion(costePoblacion.index_min());
            mejorSolucion = poblacion.row(costePoblacion.index_min());
        }

        generaciones++;
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::find(mejorSolucion == 1));    
    return ret;
}

arma::rowvec Memetico::cruceUniforme(arma::rowvec padre1, arma::rowvec padre2){

    arma::rowvec hijo(padre1.size());
    int rdn;

    for(int i=0; i<hijo.size(); i++){
        rdn = Random::get(1, 2);

        if(rdn == 1)    hijo(i) = padre1(i);
        else            hijo(i) = padre2(i);
    }

    //Reparación

    //Si sobran genes
    if(arma::sum(hijo) > m){
        double diff, diffMediaMax;
        double media = calculaDistanciaTotal(hijo)/(arma::sum(hijo));
        int gen;
        arma::uvec indices;

        while(arma::sum(hijo) > m){
            indices = arma::find(hijo == 1);

            //Buscamos el gen más alejado de la media
            diffMediaMax = 0;

            for(int j=0; j<indices.size(); j++){
                diff = 0;
                for(int i=0; i<hijo.size(); i++){
                    diff += abs(datos(i, indices(j)) - media);
                }

                //Si es el más alejado de la media lo guardamos
                if(diff > diffMediaMax){
                    diffMediaMax = diff;
                    gen = indices(j);
                }
            }
            //Eliminamos el gen del cromosoma
            hijo(gen) = 0;

            //Actualizamos la media del coste
            media = calculaDistanciaTotal(hijo) / (arma::sum(hijo));
            
        }
    }

    //Si faltan genes
    else if(arma::sum(hijo) < m){
        double diff, diffMediaMax;

        int genes = arma::sum(hijo);
        if(genes == 0)  genes = 1;

        double media = calculaDistanciaTotal(hijo)/genes;
        int gen;
        arma::uvec indices;
        while(arma::sum(hijo) < m){
            
            //Cogemos los genes que no están en uso y los que sí
            indices = arma::find(hijo == 0);

            //Buscamos el que esté más cerca de la media
            diffMediaMax = 1e8;
            for(int j=0; j<indices.size(); j++){
                diff = 0;
                for(int i=0; i<indices.size();i++){
                    diff += abs(datos(i, indices(j)) - media);
                }

                //Si es el más cercano a la media lo guardamos
                if(diff < diffMediaMax){
                    diffMediaMax = diff;
                    gen = indices(j);
                }
            }
            //Añadimos el gen al cromosoma
            hijo(gen) = 1;


            //Actualizamos la media del coste
            media = calculaDistanciaTotal(hijo) / (arma::sum(hijo));
        }
    }
    return hijo;
}

double Memetico::calculaCoste(arma::rowvec cromosoma){

    //Obtener los genes que son 1
    arma::uvec genes = arma::find(cromosoma==1);
    //Calcular la distancia entre todos los puntos
    arma::rowvec distancias(genes.size());
    double distanciaAcu;
    double distMax = -1, distMin = 1e8, distAux = 0;
    for(int i=0; i<genes.size(); i++){
        distanciaAcu = 0;
        for(int j=0; j<genes.size(); j++)
            distanciaAcu += datos(genes(i),genes(j));
        
        distancias(i) = distanciaAcu;
    }

    return distancias.max() - distancias.min();
}

double Memetico::calculaDistanciaTotal(arma::rowvec cromosoma){
    double dist = 0;
    arma::uvec indices = arma::find(cromosoma==1);

    for(int i=0; i<indices.size(); i++)
        for(int j=0; j<indices.size(); j++)
            dist += datos(indices(i), indices(j));

    return dist;
}



/*---------------------------------------*/
/*--------------Métodos BLM--------------*/
/*---------------------------------------*/


arma::irowvec Memetico::blm(arma::rowvec cromosoma){

    //Crear sel original y calcular su coste
    sel = arma::mat(m, 2);

    arma::uvec genes = arma::find(cromosoma==1);
    for(int i=0; i<genes.size(); i++){
        sel(i,0) = genes(i);
    }
    
    costeOriginal = 0;
    for(int i=0; i<m; i++)
        for (int j=0;j<m; j++)
            if(i!=j)    costeOriginal += datos(sel(i,0), sel(j,0));

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

            for (int j=0; j<sel.n_rows; j++){    
                if(!inSel && sel(j,0) == (*it_index).first)  inSel = true;    //Comprobamos si el primer índice está en sel
                if(outSel && sel(j,0) == (*it_index).second) outSel = false;  //Comprobamos si el segundo índice no está en sel
            }

            if(inSel && outSel){
                contador++;
                t++;
                seIntercambia = Int((*it_index));  //Calculamos si es una solucion aceptable y hacemos el intercambio
                if(seIntercambia){
                    index.erase(it_index);      //Borramos el intercambio del espacio de búsqueda
                }

            }
        }        
    }


    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(sel.col(0)));
    return ret;
}

bool Memetico::Int(pair <int,int> intercambio){
    
    bool seIntercambia = false;
    double distancia_v = 0;
    arma::Mat<double> sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<sel.n_rows; i++){ 
        distancia_v  += datos(intercambio.second, (int)sel_int(i,0));
        sel_int(i,1) -= datos((int)sel_int(i,0), intercambio.first);
        sel_int(i,1) += datos((int)sel_int(i,0), intercambio.second);
    }


    float diferencia_coste = (std::max(distancia_v, sel_int.col(1).max()) - std::min(distancia_v, sel_int.col(1).min())) - costeOriginal;
   
    //Condicion de aceptacion de nueva solución
    if(diferencia_coste < 0){
        seIntercambia = true;
        sel = sel_int;  //Actualizamos los costes del sel original
        costeOriginal += diferencia_coste;  //Actualizamos el coste 
        bool salir = false;

        //Hacemos el intercambio
        for(int i=0; i<sel.n_rows && !salir; i++){
            if(sel(i,0) == intercambio.first){
                sel(i) = intercambio.second, distancia_v;
                salir = true;
            }
        }    
    }
    return seIntercambia;

} 

void Memetico::mezclaIndex(vector<pair<int,int>> &index){
   
    pair<int, int> aux;
    int rdn;
   
    for(int i=0; i<index.size(); i++){
        rdn = Random::get(0, (int)index.size()-1);
        aux = index[i];
        index[i] = index[rdn];
        index[rdn] = aux;
    }
}