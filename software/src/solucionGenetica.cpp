#include "solucionGenetica.h"


Genetico::Genetico(int semilla, int seleccion, arma::mat dat, tipo_cruce c){
    //Almacenar datos y fijar tipo de algoritmo y cruce
    datos = dat;
    cruce = c;
    m = seleccion;

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


arma::irowvec Genetico::generacional(){

    //Variables del algoritmo
    int t=0;
    int numCruces = 0.7 * m / 2;
    int numMutaciones = 0.1 * (numCromosomas * datos.n_rows);
    arma::mat pob_prima(numCromosomas, datos.n_rows);
    arma::rowvec mejorSolucion = poblacion.row(costePoblacion.index_min());
    double mejorCoste = costePoblacion.min();

    //Definición de variables auxiliares de los bucles
    int rdn1, rdn2, rdn3;
    arma::rowvec hijo1, hijo2;
    arma::uvec genes1, genes0;
    arma::mat hijos;
    

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
        

        //Operación de cruce
        if(cruce == tipo_cruce::UNIFORME){

            for (int i=0; i<numCruces-1; i+=2){
                //Obtenemos los dos hijos de dos padres
                hijo1 = cruceUniforme(pob_prima.row(i), pob_prima.row(i+1));
                hijo2 = cruceUniforme(pob_prima.row(i), pob_prima.row(i+1));
                
                //Actualizamos la poblacion
                pob_prima.row(i) = hijo1;
                pob_prima.row(i+1) = hijo2;
            }

        }
        
        else if(cruce == tipo_cruce::POSICION){

           for(int i=0; i<numCruces-1; i+=2){
                //Obtenemos los hijos
                hijos = crucePosicion(pob_prima.row(i), pob_prima.row(i+1));

                //Actualizamos la población
                pob_prima.row(i) = hijos.row(0);
                pob_prima.row(i+1) = hijos.row(1);
           }

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
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::find(mejorSolucion == 1));    
    return ret;
    

}

arma::irowvec Genetico::estacionario(){

    //Variables del algoritmo
    int t=0;

    //Variables auxiliares para los bucels
    int rdn1, rdn2, rdn3, rdn4, indice;
    double costeHijo1, costeHijo2;
    arma::rowvec padre1, padre2, hijo1, hijo2;
    arma::uvec genes1, genes0;
    arma::mat hijos;
    

    while (t<100000){

        //Operador de selección
        
        //Padre 1
        do{
            rdn1 = Random::get(0, numCromosomas-1);
            rdn2 = Random::get(0, numCromosomas-1);
        }while(rdn1 == rdn2);

        //Torneo binario
        if(costePoblacion(rdn1) < costePoblacion(rdn2))
            padre1 = poblacion.row(rdn1);
        
        else
            padre1 = poblacion.row(rdn2);


        //Padre 2
        do{
            rdn1 = Random::get(0, numCromosomas-1);
            rdn2 = Random::get(0, numCromosomas-1);
        }while(rdn1 == rdn2);

        //Torneo binario
        if(costePoblacion(rdn1) < costePoblacion(rdn2))
            padre2 = poblacion.row(rdn1);
        
        else
            padre2 = poblacion.row(rdn2);
        
        hijo1 = padre1;
        hijo2 = padre2;

        rdn1 = Random::get(0, 100);
        rdn2 = Random::get(0,100);

        //Operador de cruce con probabilidad del 70%
        if(cruce == tipo_cruce::UNIFORME){
            if(rdn1 <= 70) hijo1 = cruceUniforme(padre1, padre2);
            if(rdn2 <= 70) hijo2 = cruceUniforme(padre1, padre2);
        }

        else if(cruce == tipo_cruce::POSICION){
            hijos = crucePosicion(padre1, padre2);
            if(rdn1 <= 70) hijo1 = hijos.row(0);
            if(rdn2 <= 70) hijo2 = hijos.row(1);
        }


        //Operador de mutacion con probabilidad del 10%
        rdn1 = Random::get(0, 100);
        rdn2 = Random::get(0, 100);

        if(rdn1 <= 10){
            //Genes aleatorios a intercambiar
            genes1 = arma::find(hijo1 == 1);
            genes0 = arma::find(hijo1 == 0);
            
            rdn3 = Random::get(0, (int)(genes1.n_cols));
            rdn4 = Random::get(0, (int)(genes0.n_cols));

            //Intercambiar genes
            hijo1.swap_cols(genes1(rdn3), genes0(rdn4));
        }

        if(rdn2 <= 10){
            //Genes aleatorios a intercambiar
            genes1 = arma::find(hijo2 == 1);
            genes0 = arma::find(hijo2 == 0);
            
            rdn3 = Random::get(0, (int)(genes1.n_cols));
            rdn4 = Random::get(0, (int)(genes0.n_cols));

            //Intercambiar genes
            hijo2.swap_cols(genes1(rdn3), genes0(rdn4));
        }

        //Operador de reemplazo

        //Sustituyo los dos peores hijos con los nuevos
        costeHijo1 = calculaCoste(hijo1);
        costeHijo2 = calculaCoste(hijo2);
        
        if(costePoblacion.max() > costeHijo1){
            indice = costePoblacion.index_max();
            costePoblacion(indice) = costeHijo1;
            poblacion.row(indice) = hijo1;
        }

        if(costePoblacion.max() > costeHijo2){
            indice = costePoblacion.index_max();
            costePoblacion(indice) = costeHijo2;
            poblacion.row(indice) = hijo2;
        }       
        

        //Aumento las evaluaciones
        t+=2;
    }

    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::find(poblacion.row(costePoblacion.index_min()) == 1));    
    return ret;


}

arma::rowvec Genetico::cruceUniforme(arma::rowvec padre1, arma::rowvec padre2){

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

arma::mat Genetico::crucePosicion(arma::rowvec padre1, arma::rowvec padre2){

    //Creamos un hijo que preserve los valores comunes de los padres
    arma::rowvec hijo(padre1.size(), arma::fill::zeros); //Inicializar el hijo con todo a cero
    arma::uvec comunes = arma::find(padre1 == padre2);  //indices de los elementos comunes de ambos padres
    hijo.elem(comunes) = padre1.elem(comunes);          //asignacion en las mismas posiciones de los elementos comunes

    //Asignacion aleatoria del resto del padre 1
    arma::uvec no_comunes = arma::find(padre1 != padre2);   //indices de los elementos que no coincide, el resto de elementos
    arma::rowvec padre_des1 = padre1.elem(no_comunes).as_row();  
    arma::rowvec padre_des2 = padre1.elem(no_comunes).as_row();       

    //Desordeno los restos del padre
    int rdn;

    for(int i=0; i<padre_des1.size(); i++){
        rdn = Random::get(0, (int)(padre_des1.size()-1));
        padre_des1.swap_cols(rdn, i);

        rdn = Random::get(0, (int)(padre_des2.size()-1));
        padre2.swap_cols(rdn, i);
    }

    //Asignación a los dos hijos
    arma::rowvec hijo1 = hijo;
    arma::rowvec hijo2 = hijo;

    hijo1.elem(no_comunes) = padre_des1;
    hijo2.elem(no_comunes) = padre_des2;


    return arma::join_cols(hijo1, hijo2);
}


double Genetico::calculaCoste(arma::rowvec cromosoma){

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

double Genetico::calculaDistanciaTotal(arma::rowvec cromosoma){
    double dist = 0;
    arma::uvec indices = arma::find(cromosoma==1);

    for(int i=0; i<indices.size(); i++)
        for(int j=0; j<indices.size(); j++)
            dist += datos(indices(i), indices(j));

    return dist;
}
