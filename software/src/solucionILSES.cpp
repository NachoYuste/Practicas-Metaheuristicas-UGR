
#include "solucionILSES.h"

SolucionILSES::SolucionILSES(int semilla, int m_sel, arma::mat dat){
    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    m = m_sel;
    datos = dat;
}

arma::irowvec SolucionILSES::ILSES(){
    
    //Generación de solución inicial
    arma::mat s0(generaSelAleatoria());
    arma::mat mejor_sel(BES(s0));

    //Declaración de variables para el bucle
    arma::mat sprima, sprima2;
    double coste_mejor_sel = mejor_sel.col(1).max() - mejor_sel.col(1).min();
    double coste_sprima2;

    for(int i=0; i<10; i++){
        //Operador de mutación
        sprima = mutacion(mejor_sel);   

        //Aplico la búsqueda local al sel mutado
        sprima2 = BES(sprima);

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

arma::mat SolucionILSES::BES(arma::mat sel){
       //Establecemos los valores iniciales
    int max_vecinos = 10 * datos.n_rows;
    int max_exitos = 0.1 * max_vecinos;

    //Calculo la temperatura de sel
    T = (0.3 * (sel.col(1).max() - sel.col(1).min())) / (-log(0.3));  

    //Calculo valor beta para las actualizanes de temperatura
    double beta = (T - Tf) / ((100000/max_vecinos)*T*Tf);

    //Almacenamos la "mejor solucion"
    arma::mat mejor_solucion = sel;
    double coste_mejor_sel = mejor_solucion.col(1).max() - mejor_solucion.col(1).min();
    
    int contador = 0; //contador de iteraciones
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
            seIntercambia = Int(sel, T * numEnfriamientos, coste_sel, interA, interB, datos);

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

    }while(Tf < T && contador < 10000 && exitos!=0);


    //Devolvemos la mejor solución
    return mejor_solucion;

}


bool SolucionILSES::Int(arma::mat &sel, double delta, double & coste_sel, int interA, int interB, arma::mat datos){
    
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

arma::mat SolucionILSES::generaSelAleatoria(){

    arma::mat sel;
    double coste_inicial;
    //Iteramos hasta que T-inicial sea mayor que T-final
    do{
        //Obtenemos una solución inicial aleatoria y calculamos su coste
        sel = arma::mat(1,2);
        sel(0) = Random::get(0, (int)(datos.n_rows-1)), 0;
        
        while(sel.n_rows < m){
            
            int rdn = Random::get(0, (int)(datos.n_rows-1));

            if(arma::all(sel.col(0)!=rdn))
                sel.insert_rows(sel.n_rows, arma::rowvec({(double) rdn, 0}));

        }

        for(int i=0; i<m; i++)
            for (int j=0;j<m; j++)
                sel(i,1) += datos(sel(i,0), sel(j,0));


        //Fijamos el valor inicial de la temperatura
        coste_inicial = (sel.col(1).max() - sel.col(1).min());
        
        //En el caso de que la función objetivo sea de valor 0, asignamos 1 al coste inicial para que T no sea nulo
        //Esto sirve principalmente para m=2
        if(coste_inicial == 0) coste_inicial = 1;   

        T = (0.3 * coste_inicial) / (-log(0.3));  
    }while (Tf > T);

    return sel;
}

arma::mat SolucionILSES::mutacion(arma::mat sel){
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
            if(sel_int(i,0) != interA) distancia_v  += datos(interB, (int)sel_int(i,0));
            sel_int(i,1) -= datos((int)sel_int(i,0), interA);
            sel_int(i,1) += datos((int)sel_int(i,0), interB);
        }

        //Actualizamos los costes del sel original
        sel = sel_int;

        //Hacemos el intercambio
        salir = false;
        for(int i=0; i<m && !salir; i++){
            if(sel(i,0) == interA){
                sel.row(i) = arma::rowvec({(double)interB, distancia_v});
                salir = true;
            }
        }
    }

    return sel;
}
