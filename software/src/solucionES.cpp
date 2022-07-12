#include "solucionES.h"

arma::irowvec ES(int semilla, int m, arma::mat datos){

    //Fijar la semilla para los números aleatorios
    long int seed = semilla;
    Random::seed(seed);

    //Establecemos los valores iniciales
    int max_vecinos = 10 * datos.n_rows;
    int max_exitos = 0.1 * max_vecinos;

    double T;
    double Tf = pow(10, -3);
    double coste_inicial;

    //Iteramos hasta que T-inicial sea mayor que T-final
    arma::mat sel;
    do{
        //Obtenemos una solución inicial aleatoria y calculamos su coste
        sel.clear();
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

    //Calculo valor beta para las actualizanes de temperatura
    double beta = (T - Tf) / ((100000/max_vecinos)*T*Tf);

    //Almacenamos la "mejor solucion"
    arma::mat mejor_solucion = sel;
    double coste_mejor_sel = mejor_solucion.col(1).max() - mejor_solucion.col(1).min();
    
    int contador = 0; //contador de iteraciones
    int k = 1; //contador de enfriamientos
    arma::rowvec sel_prima; //Sel tras el intercambio
    
    //Definición de variables auxiliares para los bucles
    int interA, interB;     //Indices para el intercambio
    bool inSel, outSel;     //Booleanos auxiliares para el intercambio
    bool seIntercambia = false;
    int exitos = 0;
    double coste_sel = coste_mejor_sel;
    double * ref_coste_sel = &coste_sel;
    

    do{
        
        exitos = 0;
        
        //Condicion de enfriamiento
        //Cuando se hayan generado max_vecinos o se hayan alcanzado max_exitos
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
            seIntercambia = Int(sel, T * k, ref_coste_sel, interA, interB, datos);
            
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
        k++;
        T = T / (1 + (beta*T));

    }while(Tf < T && contador < 100000 && exitos!=0);


    //Devolvemos la mejor solución
    arma::irowvec ret = arma::conv_to<arma::irowvec>::from(arma::floor(mejor_solucion.col(0)));
    return ret;
}

bool Int(arma::mat &sel, double delta, double * coste_sel, int interA, int interB, arma::mat datos){

    bool seIntercambia = false;
    double distancia_v = 0;
    arma::mat sel_int = sel;

    //Cálculo factorizado de coste de intercambio
    for(int i=0; i<sel.n_rows; i++){ 
        if(sel_int(i,0) != interA)  distancia_v  += datos(interB, (int)sel_int(i,0));
        sel_int(i,1) -= datos((int)sel_int(i,0), interA);
        sel_int(i,1) += datos((int)sel_int(i,0), interB);
    }

    float diferencia_coste = (std::max(distancia_v, sel_int.col(1).max()) - std::min(distancia_v, sel_int.col(1).min())) - *coste_sel;
    float rdn = Random::get(0.0, 1.0);
    float prob = exp(((-1)*diferencia_coste)/delta);
    bool pasaProbAceptacion = rdn <= prob;


    //Condicion de aceptacion de nueva solución
    if(diferencia_coste < 0 || pasaProbAceptacion){
        seIntercambia = true;
        sel = sel_int;  //Actualizamos los costes del sel original
        *coste_sel += diferencia_coste;
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

