#include "lectorDatos.h"

arma::Mat<double>  lectorDatos(int numFichero, int * numSeleccion){
    
    //Apertura de fichero
    string archivo = getNombreArchivo(numFichero);
    ifstream file(archivo);
    string cadena;

    
    //Comprobación apertura de fichero
    if (!file)    {
        cout << "El achivo " << archivo << " no se ha podido abrir" << endl; 
        exit(EXIT_FAILURE);
    }

    getline(file,cadena);

    vector<string> descripcion = split_iterator(cadena);

    int numElementos = stoi(descripcion[0]);
    *numSeleccion = stoi(descripcion[1]);

    //Matriz de datos
    arma::Mat<double>  datos (numElementos, numElementos);
    
    //Variables auxiliares para creación de matriz
    vector<string> fila;
    int posX, posY;
    double valor;

    while (getline(file, cadena)){
        fila = split_iterator(cadena);
        posX = stoi(fila[0]);
        posY = stoi(fila[1]);
        valor = stod(fila[2]);

        //Completamos la matriz por conveniencia posterior
        datos(posX, posY) = valor;
        datos(posY, posX) = valor;
    }


    file.close();

    
    return datos;
}



vector<string> split_iterator(const string& str) {
  vector<string> resultado;

  string::const_iterator itBegin = str.begin();
  string::const_iterator itEnd   = str.end();

  int numItems = 1;
  for( string::const_iterator it = itBegin; it!=itEnd; ++it )
    numItems += *it==' ';

  resultado.reserve(numItems);

  for( string::const_iterator it = itBegin; it!=itEnd; ++it )
  {
    if( *it == ' ' )
    {
      resultado.push_back(string(itBegin,it));
      itBegin = it+1;
    }
  }

  if( itBegin != itEnd )
    resultado.push_back(string(itBegin,itEnd));

  return resultado;
}


string getNombreArchivo(int numFichero){
    
    string nombre;

    if(numFichero<=5)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n25_m2.txt";
    
    else if(numFichero>5 && numFichero<=10)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n25_m7.txt";

    else if(numFichero>10 && numFichero<=15)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n50_m5.txt";

    else if(numFichero>15 && numFichero<=20)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n50_m15.txt";

    else if(numFichero>20 && numFichero<=25)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n100_m10.txt";

    else if(numFichero>25 && numFichero<=30)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n100_m30.txt";

    else if(numFichero>30 && numFichero<=35)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n125_m12.txt";

    else if(numFichero>35 && numFichero<=40)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n125_m37.txt";

    else if(numFichero>40 && numFichero<=45)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n150_m15.txt";

    else if(numFichero>45 && numFichero<=50)
        nombre="../datos/GKD-b_"+to_string(numFichero)+"_n150_m45.txt";

 
    return nombre;
}

double calculaDiff(arma::irowvec sol, arma::Mat<double>  datos){

    double diff;
    arma::rowvec  coste = arma::rowvec (sol.size(), arma::fill::zeros);
    for(int i=0; i<sol.size(); i++)
        for(int j=0; j<sol.size();j++)
            coste(i)+=datos(sol(i), sol(j));
    
    return (coste.max() - coste.min());
}
