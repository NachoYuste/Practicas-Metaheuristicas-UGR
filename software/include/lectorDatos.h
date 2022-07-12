#ifndef LECTORDATOS_H_INCLUDED
#define LECTORDATOS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

using namespace std;

//Devuelve la dirección completa del archivo a partir de su número identificador
string getNombreArchivo(int numFichero);

//Separa las palabras de un string
vector<string> split_iterator(const string& str);

//Devuelve la matriz de los datos del fichero
arma::Mat<double> lectorDatos(int numFichero, int * numSeleccion);

double calculaDiff(arma::irowvec sol, arma::Mat<double>  datos);

#endif