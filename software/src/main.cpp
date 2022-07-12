#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <cstdlib>

#include "lectorDatos.h"
#include "solucionGreedy.h"
#include "solucionBLM.h"
#include "solucionGenetica.h"
#include "solucionMemetica.h"
#include "solucionES.h"
#include "solucionBMB.h"
#include "solucionILS.h"
#include "solucionILSES.h"
#include "solucionCVOA.h"
#include "solucionCVOA_cec17.h"

void seRepiten (arma::irowvec solucion){
    bool ret = false;

    for (int i=0; i<solucion.size()&& !ret; i++)
        for(int j=0; j<solucion.size(); j++)
            if(solucion(i) == solucion(j) && i!=j)  ret = true;


    if(ret)
        cout << "SE REPITEN" << endl;

}

int main(int argc, char ** argv){

    char opcion;

    do {
        cout << "Experimento parcial de un algoritmo -> P" << endl;
        cout << "Experimento completo                -> C" << endl;
        cin >> opcion;
    } while (opcion != 'P' && opcion != 'C');
    

    string alg, arch, sem;

    if(opcion == 'P'){

       /*do{
        cout << "Introduce el algoritmo que desea probar, el número del archivo de datos y la semilla separados por un espacio" << endl;
        cout << "BLM -> B\tGreedy -> G\tGenetico -> GEN\tMemetico -> MEM\tEnfriamiento Simulado -> ES\tBusqueda multiarranque -> BMB\t ILS -> ILS\tILS-ES -> ILSES" << endl;
        cout << "| Num archivo -> [1-50]\t| Semilla->  Numero entero |" << endl;
        cin >> alg >> arch >> sem;
        } while(alg!="B" && alg!="G" && alg!="GEN" && alg!="MEM" && alg!="ES" && alg!="BMB" && alg!="ILS" && stoi(arch)<1 && stoi(arch)>50);

        arma::irowvec solucion;
        int numSeleccion;
        arma::Mat<double> datos = lectorDatos(stoi(arch), &numSeleccion);
        
        if(alg == "B") solucion = blm(stoi(sem), numSeleccion, datos);

        else if(alg == "G") solucion = greedy(stoi(sem), numSeleccion, datos);

        else if(alg == "GEN"){
            string cruce, algoritmo;
            tipo_cruce tipoCruce;
            do{
                cout << "Introduce el tipo de algoritmo genetico y el tipo de cruce separados por un espacio:" <<endl;
                cout << "Estacionario -> E\tGeneracional -> G\t| Uniforme -> U\tPosicion -> P"<<endl;
                cin >> algoritmo >> cruce;
            }while(algoritmo!="E" && algoritmo!="G" && cruce!="U" && cruce!="P");
        
            if(cruce == "U") tipoCruce = tipo_cruce::UNIFORME;
            if(cruce == "P") tipoCruce = tipo_cruce::POSICION;
            
            Genetico gen(stoi(sem), numSeleccion, datos, tipoCruce);
            
            if(algoritmo == "E") solucion = gen.estacionario();
            if(algoritmo == "G") solucion = gen.generacional();
            
            solucion = gen.estacionario();
        }

        else if(alg == "MEM"){
            string tipoAlg;
            do{
                cout << "Introduce el tipo de algoritmo genetico y el tipo de cruce separados por un espacio:" <<endl;
                cout << "AM(10, 1.0) -> A\tAM(10, 0.1) -> B\tAM(10, 0.1mej) -> C"<<endl;
                cin >> tipoAlg;
            }while(tipoAlg!="A" && tipoAlg!="B" && tipoAlg!="C");

            tipo_ejecucion tipoEje;
            if(tipoAlg == "A")  tipoEje = tipo_ejecucion::AM1;
            if(tipoAlg == "B")  tipoEje = tipo_ejecucion::AM2;
            if(tipoAlg == "C")  tipoEje = tipo_ejecucion::AM3;

            Memetico mem(stoi(sem), numSeleccion, datos, tipoEje);
            solucion = mem.generacional();
        }
        else if(alg == "ES") solucion = ES(stoi(sem), numSeleccion, datos);

        else if(alg == "BMB"){
            SolucionBMB bmb(stoi(sem), numSeleccion, datos);
            solucion = bmb.BMB();
        }

        else if(alg == "ILS"){
            SolucionILS ils(stoi(sem), numSeleccion, datos);
            solucion = ils.ILS();
        }

        else if(alg == "ILSES"){
            SolucionILSES ils(stoi(sem), numSeleccion, datos);
            solucion = ils.ILSES();
        }

        cout << endl << endl;
        cout << "Solucion:" << endl << solucion <<endl;
        cout << "Coste: " << calculaDiff(solucion, datos) << endl;
        */
        cout << "Introduce el número del archivo [1-50]: ";
        string arch;
        int numSeleccion;
        cin >> arch;
        arma::Mat<double> datos = lectorDatos(stoi(arch), &numSeleccion);
        arma::irowvec solucion;
        CVOA obj(1502123, numSeleccion, datos, tipo_memetico::NONE);
        solucion = obj.cvoa();
        seRepiten(solucion);
        //131
        cout << endl << endl;
        cout << "Solucion:" << endl << solucion <<endl;
        cout << "Coste: " << calculaDiff(solucion, datos) << endl;

        
       /*int seed;
       cin >> seed;
       CVOA_CEC cvoacec(seed, 10, false);
       cvoacec.cvoa_cec();*/

    }

    else if (opcion == 'C'){
/*    
        string datos_blm_name = "../datos/datos_blm.csv";
        ofstream datos_blm;
        datos_blm.open(datos_blm_name, fstream::out);
        
        if(!datos_blm.is_open()){
            cout << "Error al abrir el fichero datos_blm.csv" << endl;
            exit(0);
        }
    
        string datos_greedy_name = "../datos/datos_greedy.csv";
        ofstream datos_greedy;
        datos_greedy.open(datos_greedy_name, fstream::out);

    
        if(!datos_greedy.is_open()){
            cout << "Error al abrir el fichero datos_greedy.csv" << endl;
            exit(0);
        }

        //Archivo de salida para Genétio-Estacionario-Uniforme
        string datos_gen_est_uni_name = "../datos/datos_gen_est_uni.csv";
        ofstream datos_gen_est_uni;
        datos_gen_est_uni.open(datos_gen_est_uni_name, fstream::out);

        if(!datos_gen_est_uni.is_open()){
            cout << "Error al abrir el fichero" + datos_gen_est_uni_name << endl;
            exit(0);
        }

        //Archivo de salida para Genétio-Estacionario-Posicion
        string datos_gen_est_pos_name = "../datos/datos_gen_est_pos.csv";
        ofstream datos_gen_est_pos;
        datos_gen_est_pos.open(datos_gen_est_pos_name, fstream::out);

        if(!datos_gen_est_pos.is_open()){
            cout << "Error al abrir el fichero" + datos_gen_est_pos_name << endl;
            exit(0);
        }
/*
        //Archivo de salida para Genétio-Generacional-Uniforme
        string datos_gen_gen_uni_name = "../datos/datos_gen_gen_uni.csv";
        ofstream datos_gen_gen_uni;
        datos_gen_gen_uni.open(datos_gen_gen_uni_name, fstream::out);

        if(!datos_gen_gen_uni.is_open()){
            cout << "Error al abrir el fichero" + datos_gen_gen_uni_name << endl;
            exit(0);
        }

        //Archivo de salida para Genétio-Generacional-Posicion
        string datos_gen_gen_pos_name = "../datos/datos_gen_gen_pos.csv";
        ofstream datos_gen_gen_pos;
        datos_gen_gen_pos.open(datos_gen_gen_pos_name, fstream::out);

        if(!datos_gen_gen_pos.is_open()){
            cout << "Error al abrir el fichero" + datos_gen_gen_pos_name << endl;
            exit(0);
        }        


        //Archivo de salida para AM(10, 1.0)
        string datos_AM1_name = "../datos/datos_AM1.csv";
        ofstream datos_AM1;
        datos_AM1.open(datos_AM1_name, fstream::out);

        if(!datos_AM1.is_open()){
            cout << "Error al abrir el fichero" + datos_AM1_name << endl;
            exit(0);
        }  

        //Archivo de salida para AM(10, 0.1)
        string datos_AM2_name = "../datos/datos_AM2.csv";
        ofstream datos_AM2;
        datos_AM2.open(datos_AM2_name, fstream::out);

        if(!datos_AM2.is_open()){
            cout << "Error al abrir el fichero" + datos_AM2_name << endl;
            exit(0);
        }              

        //Archivo de salida para AM(10, 0.1mej)
        string datos_AM3_name = "../datos/datos_AM3.csv";
        ofstream datos_AM3;
        datos_AM3.open(datos_AM3_name, fstream::out);

        if(!datos_AM3.is_open()){
            cout << "Error al abrir el fichero" + datos_AM3_name << endl;
            exit(0);
        }          
        

        //Archivo de salida para ES
        string datos_ES_name = "../datos/datos_ES.csv";
        ofstream datos_ES;
        datos_ES.open(datos_ES_name, fstream::out);

        if(!datos_ES.is_open()){
            cout << "Error al abrir el fichero" + datos_ES_name << endl;
            exit(0);
        }

        //Archivo de salida para BMB
        string datos_BMB_name = "../datos/datos_BMB.csv";
        ofstream datos_BMB;
        datos_BMB.open(datos_BMB_name, fstream::out);

        if(!datos_BMB.is_open()){
            cout << "Error al abrir el fichero" + datos_BMB_name << endl;
            exit(0);
        }

        //Archivo de salida para ILS
        string datos_ILS_name = "../datos/datos_ILS.csv";
        ofstream datos_ILS;
        datos_ILS.open(datos_ILS_name, fstream::out);

        if(!datos_ILS.is_open()){
            cout << "Error al abrir el fichero" + datos_ILS_name << endl;
            exit(0);
        }

        //Archivo de salida para ILS-ES
        string datos_ILSES_name = "../datos/datos_ILSES.csv";
        ofstream datos_ILSES;
        datos_ILSES.open(datos_ILSES_name, fstream::out);

        if(!datos_ILSES.is_open()){
            cout << "Error al abrir el fichero" + datos_ILSES_name << endl;
            exit(0);
        }*/

        //Archivo de salida para CVOA
        string datos_CVOA_name = "../datos/datos_CVOA.csv";
        ofstream datos_CVOA;
        datos_CVOA.open(datos_CVOA_name, fstream::out);

        if(!datos_CVOA.is_open()){
            cout << "Error al abrir el fichero" + datos_CVOA_name << endl;
            exit(0);
        }

        //Archivo de salida para CVOA-BLM
        string datos_CVOA_BLM_name = "../datos/datos_CVOA_BLM.csv";
        ofstream datos_CVOA_BLM;
        datos_CVOA_BLM.open(datos_CVOA_BLM_name, fstream::out);

        if(!datos_CVOA_BLM.is_open()){
            cout << "Error al abrir el fichero" + datos_CVOA_BLM_name << endl;
            exit(0);
        }

        //Archivo de salida para CVOA-ES
        string datos_CVOA_ES_name = "../datos/datos_CVOA_ES.csv";
        ofstream datos_CVOA_ES;
        datos_CVOA_ES.open(datos_CVOA_ES_name, fstream::out);

        if(!datos_CVOA_ES.is_open()){
            cout << "Error al abrir el fichero" + datos_CVOA_ES_name << endl;
            exit(0);
        }

        

        vector <int> semillas = {23,32,46,69,72};

        for(int i=1; i<=50; i++){    //Todos los archivos de datos

            int numArchivo = i;

            int numSeleccion;
            arma::Mat<double> datos = lectorDatos(numArchivo, &numSeleccion);
            
            /*
            //EJECUCION COMPLETA PARA BLM
            arma::irowvec solucion_blm;

            double tiempoTotal_blm = 0, costeTotal_blm = 0;

            for(int j=0; j<semillas.size(); j++){
                auto inicioBLM = chrono::high_resolution_clock::now();
                solucion_blm = blm(semillas[j], numSeleccion, datos);
                auto finBLM = chrono::high_resolution_clock::now();

                chrono::milliseconds tiempoBLM =  chrono::duration_cast<std::chrono::milliseconds   >(finBLM - inicioBLM);

                costeTotal_blm += calculaDiff(solucion_blm, datos);
                tiempoTotal_blm +=tiempoBLM.count();
            }

            
            datos_blm << costeTotal_blm/semillas.size() << "," << tiempoTotal_blm/semillas.size() << endl;

            //EJECUCION COMPLETA PARA GREEDY

            arma::irowvec solucion_greedy;
            double tiempoTotal_greedy = 0, costeTotal_greedy = 0;

            for(int j=0; j<semillas.size(); j++){
                auto inicioGreedy = chrono::high_resolution_clock::now();
                solucion_greedy = greedy(semillas[j], numSeleccion, datos);
                auto finGreedy = chrono::high_resolution_clock::now();

                chrono::milliseconds tiempoGreedy = chrono::duration_cast<std::chrono::milliseconds>(finGreedy - inicioGreedy);

                costeTotal_greedy += calculaDiff(solucion_greedy, datos);
                tiempoTotal_greedy += tiempoGreedy.count();
            }

            datos_greedy << costeTotal_greedy/semillas.size() << "," << tiempoTotal_greedy/semillas.size() << endl;
            
            arma::irowvec solucion;
            double coste, tiempo;

           /*
            //Genetico-Estacionario-Uniforme
            auto inicio = chrono::high_resolution_clock::now();
            Genetico estacionario_uniforme(23, numSeleccion, datos, tipo_cruce::UNIFORME);
            solucion = estacionario_uniforme.estacionario();
            auto fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_gen_est_uni << coste << "," << tiempo << endl;


            //Genetico-Estacionario-Posicion
            inicio = chrono::high_resolution_clock::now();
            Genetico estacionario_posicion(23, numSeleccion, datos, tipo_cruce::POSICION);
            solucion = estacionario_posicion.estacionario();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_gen_est_pos << coste << "," << tiempo << endl;

            //Genetico-Generacional-Uniforme
            inicio = chrono::high_resolution_clock::now();
            Genetico generacional_uniforme(23, numSeleccion, datos, tipo_cruce::UNIFORME);
            solucion = generacional_uniforme.generacional();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_gen_gen_uni << coste << "," << tiempo << endl;

            //Genetico-Generacional-Posicion
            inicio = chrono::high_resolution_clock::now();
            Genetico generacional_posicion(23, numSeleccion, datos, tipo_cruce::POSICION);
            solucion = generacional_posicion.generacional();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_gen_gen_pos << coste << "," << tiempo << endl;
            
            //AM1
            auto inicio = chrono::high_resolution_clock::now();
            Memetico memeticoAM1(23, numSeleccion, datos, tipo_ejecucion::AM1);
            solucion = memeticoAM1.generacional();
            auto fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_AM1 << coste << "," << tiempo << endl;


            //AM2
            inicio = chrono::high_resolution_clock::now();
            Memetico memeticoAM2(23, numSeleccion, datos, tipo_ejecucion::AM2);
            solucion = memeticoAM2.generacional();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_AM2 << coste << "," << tiempo << endl;


            //AM3
            inicio = chrono::high_resolution_clock::now();
            Memetico memeticoAM3(23, numSeleccion, datos, tipo_ejecucion::AM3);
            solucion = memeticoAM3.generacional();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_AM3 << coste << "," << tiempo << endl;

            

            //ES
            auto inicio = chrono::high_resolution_clock::now();
            solucion = ES(23, numSeleccion, datos);
            auto fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_ES << coste << "," << tiempo << endl;

            //BMB
            inicio = chrono::high_resolution_clock::now();
            SolucionBMB bmb (23, numSeleccion, datos);
            solucion = bmb.BMB();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_BMB << coste << "," << tiempo << endl;

            //ILS
            inicio = chrono::high_resolution_clock::now();
            SolucionILS ils (23, numSeleccion, datos);
            solucion = ils.ILS();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_ILS << coste << "," << tiempo << endl;
            

            //ILS-ES
            inicio = chrono::high_resolution_clock::now();
            SolucionILSES ilses (23, numSeleccion, datos);
            solucion = ilses.ILSES();
            fin = chrono::high_resolution_clock::now();

            coste = calculaDiff(solucion, datos);
            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_ILSES << coste << "," << tiempo << endl;

            */

            //EJECUCION COMPLETA PARA CVOA
            arma::irowvec solucion;


            auto inicio = chrono::high_resolution_clock::now();
            CVOA cvoa(1502123, numSeleccion, datos, tipo_memetico::NONE);
            solucion = cvoa.cvoa();
            auto fin = chrono::high_resolution_clock::now();

            double tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_CVOA << calculaDiff(solucion, datos) << "," << tiempo << endl;       

            //EJECUCION COMPLETA PARA CVOA-BLM
            inicio = chrono::high_resolution_clock::now();
            CVOA cvoablm(1502123, numSeleccion, datos, tipo_memetico::BLM);
            solucion = cvoablm.cvoa();
            fin = chrono::high_resolution_clock::now();

            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_CVOA_BLM<< calculaDiff(solucion, datos) << "," << tiempo << endl;    

            //EJECUCION COMPLETA PARA CVOA-ES
            inicio = chrono::high_resolution_clock::now();
            CVOA cvoaes(1502123, numSeleccion, datos, tipo_memetico::ES);
            solucion = cvoaes.cvoa();
            fin = chrono::high_resolution_clock::now();

            tiempo = (chrono::duration_cast<std::chrono::milliseconds>(fin - inicio)).count();

            datos_CVOA_ES<< calculaDiff(solucion, datos) << "," << tiempo << endl;    
            
            
            cout << "Archivo: " << i << " terminado." << endl;
       }
       
        
        /*datos_blm.close();
        datos_greedy.close();
        datos_gen_est_pos.close();
        datos_gen_est_uni.close();
        datos_gen_gen_pos.close();
        datos_gen_gen_uni.close();
        datos_AM1.close();
        datos_AM2.close();
        datos_AM3.close();
        datos_ES.close();
        datos_BMB.close();
        datos_ILS.close();
        datos_ILSES.close();*/
        datos_CVOA.close();
    }
}

