#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

///////////////// PROTOTIPOS DE FUNCIONES /////////////////
void obtenerBins(int[], double[]);
void obtenerKerma(int, int, double[][179], double[][179]);
void obtenerAirKermaStrength(double &, double &);
void calcularFuncionGeometrica(int, double[], double[], double[][179]);
void calcularFuncionRadial(int, double[][179], double[][179], double, double, double[]);
void calcularFuncionAnisotropia(int, double[][179], double[][179], double[][179]);
void escribirFicheroConstantes(double, double, double, double, double);
void escribirFicheroGeometrica(int, double[], double[], double[][179], double);
void escribirFicheroRadial(int, double[], double[]);
void escribirFicheroAnisotropia(int, double[], double[], double[][179], double);
void escribirFicheroAnisotropia2(double, double[], double[][179], double);
void reconstruirDosis(double[][179], int, double, double, double, double[][179], double[], double[][179]);
void escribirFicheroDosisReconstruida(int, double[], double[], double[][179]);
double roundit(double, double);

const double PI = 3.1415926535897;

int main(){

    /**************************** LECTURA DE DATOS DE LA SIMULACIÓN EN AGUA ****************************/

    int numBins[2]; //Contiene el número de bins de r numBins[0] y de theta numBins[1]
    double minWidth[4]; //Valor mínimo de r minWidth[0], anchura de bin de r minWidth[1], valor mínimo de th minWidth[2], anchura de bin de th minWidth[3]

    obtenerBins(numBins, minWidth);

    int rBins = numBins[0];
    int thBins = numBins[1];
    double rMin = minWidth[0];
    double rWidth = minWidth[1];
    double thMin = minWidth[2];
    double thWidth = minWidth[3];

    double rValores[rBins];
    double thValores[thBins];
    double kerma[rBins][179];
    double sigma[rBins][179];

    for(int i = 0; i < rBins; i++){
        rValores[i] = (2*rMin + rWidth)/2 + i*rWidth;
    }

    for(int j = 0; j < thBins; j++){
        thValores[j] = (2*thMin + thWidth)/2 + j*thWidth;
    }

    obtenerKerma(rBins, thBins, kerma, sigma);

    /**************************** LECTURA DE DATOS DE LA SIMULACIÓN EN AIRE ****************************/
    double kermaAire, sigmaAire;
    obtenerAirKermaStrength(kermaAire, sigmaAire);

    /********************************* CÁLCULO DE FUNCIONES DEL TG-43 **********************************/

    // Cálculo del air kerma strength con su error
    double S_k = roundit(100 * kermaAire, 3);
    double errorS_k = roundit(100 * sigmaAire, 1);

    // Cálculo de la constante de dosis
    double Lambda = roundit((kerma[9][89])/S_k, 4);
    double errorLambda = roundit(Lambda * ((sigma[9][89])/(kerma[9][89])+errorS_k/S_k), 1);

    // Cálculo de la función geométrica
    double G_L[rBins][179];
    calcularFuncionGeometrica(rBins, rValores, thValores, G_L);

    // Cálculo de la función radial
    double K_0 = kerma[9][89];
    double G_L_0 = G_L[9][89];
    double g_L[thBins];
    calcularFuncionRadial(rBins, kerma, G_L, K_0, G_L_0, g_L);

    // Cálculo de la función de anisotropía
    double F[rBins][179];
    calcularFuncionAnisotropia(rBins, kerma, G_L, F);

    // Cálculo de la dosis a partir de las funciones del TG-43
    double dosisReconstruida[rBins][179];
    reconstruirDosis(dosisReconstruida, rBins, S_k, Lambda, G_L_0, G_L, g_L, F);

    /**************************** ESCRITURA DE RESULTADOS EN FICHEROS ****************************/
    escribirFicheroConstantes(S_k, errorS_k, Lambda,errorLambda, G_L_0);
    escribirFicheroGeometrica(rBins, rValores, thValores, G_L, rWidth);
    escribirFicheroRadial(rBins, rValores, g_L);
    escribirFicheroAnisotropia(rBins, rValores, thValores, F, rWidth);
    escribirFicheroDosisReconstruida(rBins, rValores, thValores, dosisReconstruida);

    int opcion;
    double r;
    int th;

    do{
        cout << "1. Obtener K(r,th) en un punto" << endl;
        cout << "2. Obtener G_L(r,th) en un punto" << endl;
        cout << "3. Obtener g_L(r) en un punto" << endl;
        cout << "4. Obtener F(r,th) en un punto" << endl;
        cout << "5. Obtener F(th) para una distancia determinada" << endl;
        cout << "0. Salir del programa" << endl;

        cout << "\nSelecciona una opcion: ";
        cin >> opcion;
        cout << "\n\n";

        switch(opcion){
            case 1:
                cout << "\n(r,th) = ";
                cin >> r >> th;
                cout << "\n K(" << r << ", " << th << ") = " << kerma[(int)(r/rWidth+0.5) - 1][th - 1]
                    << " +- " << sigma[(int)(r/rWidth+0.5) - 1][th - 1] << "\n\n";
                break;
            case 2:
                cout << "\n(r,th) = ";
                cin >> r >> th;
                cout << "\n G_L(" << r << ", " << th << ") = " << G_L[(int)(r/rWidth+0.5) - 1][th - 1] << "\n\n";
                break;
            case 3:
                cout << "\n r = ";
                cin >> r;
                cout << "\n g_L(" << r << ") = " << g_L[(int)(r/rWidth+0.5) - 1] << "\n\n";
                break;
            case 4:
                cout << "\n(r,th) = ";
                cin >> r >> th;
                cout << "\n F(" << r << ", " << th << ") = " << F[(int)(r/rWidth+0.5) - 1][th - 1] << "\n\n";
                break;
            case 5:
                cout << "\n r = ";
                cin >> r;
                cout << endl;
                escribirFicheroAnisotropia2(r, thValores, F, rWidth);
                break;
        }
    }
    while(opcion != 0);

    return 0;
}

void obtenerBins(int numBins[], double minWidth[]){ //Lee del fichero los datos relacionados con los bins de la simulación

    ifstream ficheroLectura("tallyTrackLengthEsfericas.dat");

    //Sale del programa si ifstream no pudo abrir el archivo
    if(!ficheroLectura){
        cerr << "No se pudo abrir el archivo en modo lectura" << endl;
        exit(1);
    }

    for(int i = 1; i <= 4; i++){
        ficheroLectura.ignore(500, '\n');
    }

    ficheroLectura.ignore(3);
    ficheroLectura >> numBins[0] >> numBins[1];
    ficheroLectura.ignore(1);
    ficheroLectura.ignore(500, '\n');
    ficheroLectura.ignore(5);
    ficheroLectura >> minWidth[0] >> minWidth[1] >> minWidth[2] >> minWidth[3];

    ficheroLectura.close();
}

void obtenerKerma(int rBins, int thBins, double kerma[][179], double sigma[][179]){

    ifstream ficheroLectura("tallyTrackLengthEsfericas.dat");

    //Sale del programa si ifstream no pudo abrir el archivo
    if(!ficheroLectura){
        cerr << "No se pudo abrir el archivo en modo lectura" << endl;
        exit(1);
    }

    for(int i = 1; i <= 12; i++){
        ficheroLectura.ignore(500, '\n');
    }


    for(int j = 0; j < thBins; j++){

        for(int k = 0; k < rBins; k++){

            ficheroLectura.ignore(66); //Ignora los espacios en blanco antes del kerma
            ficheroLectura >> kerma[k][j] >> sigma[k][j]; //Lee y almacena kerma y su error
            ficheroLectura.ignore(1,'\n'); //Salta de linea
        }
    }

    ficheroLectura.close();
}

void obtenerAirKermaStrength(double &kermaAire, double &sigmaAire){

    ifstream ficheroLectura("tallyTrackLengthCilindricas.dat");

    //Sale del programa si ifstream no pudo abrir el archivo
    if(!ficheroLectura){
        cerr << "No se pudo abrir el archivo en modo lectura" << endl;
        exit(1);
    }

    for(int i = 1; i <= 12; i++){
        ficheroLectura.ignore(500, '\n');
    }

    ficheroLectura.ignore(66);

    ficheroLectura >> kermaAire >> sigmaAire;
}

void calcularFuncionGeometrica(int rBins, double rValores[], double thValores[], double G_L[][179]){

    double L;

    cout << "Introduce la longitud de la fuente en cm: ";
    cin >> L;
    cout << endl;

    for(int i = 0; i < 179; i++){

        for(int j = 0; j < rBins; j++){

            G_L[j][i]=(acos((rValores[j]*cos(thValores[i]*PI/180)-L/2)/(sqrt(pow(rValores[j],2)+pow(L/2,2)-L*rValores[j]*cos(thValores[i]*PI/180))))
            - acos((rValores[j]*cos(thValores[i]*PI/180)+L/2)/(sqrt(pow(rValores[j],2)+pow(L/2,2)+L*rValores[j]*cos(thValores[i]*PI/180)))))
            /(L * rValores[j] * sin(thValores[i]*PI/180));
        }
    }
}

void calcularFuncionRadial(int rBins, double kerma[][179], double G_L[][179], double K_0, double G_L_0, double g_L[]){

    for(int i = 0; i < rBins; i++){
        g_L[i] = (kerma[i][89] * G_L_0) / (K_0 * G_L[i][89]);
    }
}

void calcularFuncionAnisotropia(int rBins, double kerma[][179], double G_L[][179], double F[][179]){

    for(int i = 0; i < 179; i++){
            for(int j = 0; j < rBins; j++){
                F[j][i] = (kerma[j][i] * G_L[j][89]) / (kerma[j][89] * G_L[j][i]);
            }
    }
}

void escribirFicheroConstantes(double S_k, double errorS_k, double Lambda, double errorLambda, double G_L_0){

    ofstream ficheroEscritura("constantes.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "************************ CONSTANTES DEL TG-43 *************************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;

    ficheroEscritura << "Air-kerma strength:\t" << S_k << " +- " << errorS_k << "\n\n";

    ficheroEscritura << "Dose-rate constant:\t" << Lambda << " +- " << errorLambda << "\n\n";

    ficheroEscritura << "G_L(r = 1cm, th = 90º):\t" << G_L_0 << "\n\n";

    ficheroEscritura.close();
}

void escribirFicheroGeometrica(int rBins, double rValores[], double thValores[], double G_L[][179], double rWidth){
    ofstream ficheroEscritura("funcionGeometrica.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "************************* FUNCIÓN GEOMÉTRICA **************************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;
    ficheroEscritura << setw(10) << left << "r (cm)" << setw(10) << left << "th (º)" << setw(10) << left << "G_L(r,th)" << endl;
    for(int i = 0; i < 179; i++){
        for(int j = 0; j < rBins; j++){
            ficheroEscritura << left << setw(10) << rValores[j] << left << setw(10) << thValores[i] << left << setw(10) << G_L[j][i] << endl;
        }
    }

    ficheroEscritura.close();

    ofstream ficheroEscrituraMath1("funcionGeometrica1.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath1){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath1 << "***********************************************************************" << endl;
    ficheroEscrituraMath1 << "*********************** FUNCIÓN GEOMÉTRICA r=1 ************************" << endl;
    ficheroEscrituraMath1 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath1 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath1 << "{" << thValores[i] << "," << G_L[(int)(1/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath1 << ",";
        }
    }
    ficheroEscrituraMath1 << "}";

    ficheroEscrituraMath1.close();

    ofstream ficheroEscrituraMath2("funcionGeometrica2.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath2){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath2 << "***********************************************************************" << endl;
    ficheroEscrituraMath2 << "*********************** FUNCIÓN GEOMÉTRICA r=5 ************************" << endl;
    ficheroEscrituraMath2 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath2 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath2 << "{" << thValores[i] << "," << G_L[(int)(5/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath2 << ",";
        }
    }
    ficheroEscrituraMath2 << "}";

    ficheroEscrituraMath2.close();

    ofstream ficheroEscrituraMath3("funcionGeometrica3.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath3){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath3 << "***********************************************************************" << endl;
    ficheroEscrituraMath3 << "*********************** FUNCIÓN GEOMÉTRICA r=10 ************************" << endl;
    ficheroEscrituraMath3 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath3 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath3 << "{" << thValores[i] << "," << G_L[(int)(10/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath3 << ",";
        }
    }
    ficheroEscrituraMath3 << "}";

    ficheroEscrituraMath3.close();

    ofstream ficheroEscrituraMath4("funcionGeometrica4.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath4){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath4 << "***********************************************************************" << endl;
    ficheroEscrituraMath4 << "*********************** FUNCIÓN GEOMÉTRICA th=1 ************************" << endl;
    ficheroEscrituraMath4 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath4 << "{";
    for(int i = 0; i < rBins; i++){
        ficheroEscrituraMath4 << "{" << rValores[i] << "," << G_L[i][0]*pow(rValores[i],2) << "}";
        if(i != rBins-1){
            ficheroEscrituraMath4 << ",";
        }
    }
    ficheroEscrituraMath4 << "}";

    ficheroEscrituraMath4.close();

    ofstream ficheroEscrituraMath5("funcionGeometrica5.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath5){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath5 << "***********************************************************************" << endl;
    ficheroEscrituraMath5 << "*********************** FUNCIÓN GEOMÉTRICA th=45 ************************" << endl;
    ficheroEscrituraMath5 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath5 << "{";
    for(int i = 0; i < rBins; i++){
        ficheroEscrituraMath5 << "{" << rValores[i] << "," << G_L[i][44]*pow(rValores[i],2) << "}";
        if(i != rBins-1){
            ficheroEscrituraMath5 << ",";
        }
    }
    ficheroEscrituraMath5 << "}";

    ficheroEscrituraMath5.close();

    ofstream ficheroEscrituraMath6("funcionGeometrica6.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath6){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath6 << "***********************************************************************" << endl;
    ficheroEscrituraMath6 << "*********************** FUNCIÓN GEOMÉTRICA th=90 ************************" << endl;
    ficheroEscrituraMath6 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath6 << "{";
    for(int i = 0; i < rBins; i++){
        ficheroEscrituraMath6 << "{" << rValores[i] << "," << G_L[i][89]*pow(rValores[i],2) << "}";
        if(i != rBins-1){
            ficheroEscrituraMath6 << ",";
        }
    }
    ficheroEscrituraMath6 << "}";

    ficheroEscrituraMath6.close();
}

void escribirFicheroRadial(int rBins, double rValores[], double g_L[]){
    ofstream ficheroEscritura("funcionRadial.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "*************************** FUNCIÓN RADIAL ****************************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;
    ficheroEscritura << setw(10) << left << "r (cm)" << setw(10)  << left << "g_L(r)" << endl;
    for(int j = 0; j < rBins; j++){
        ficheroEscritura << left << setw(10) << rValores[j] << left << setw(10) << g_L[j] << endl;
    }

    ficheroEscritura.close();

    ofstream ficheroEscrituraMath("funcionRadialMathematica.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath << "***********************************************************************" << endl;
    ficheroEscrituraMath << "*************************** FUNCIÓN RADIAL ****************************" << endl;
    ficheroEscrituraMath << "***********************************************************************\n" << endl;
    ficheroEscrituraMath << "{";
    for(int j = 0; j < rBins; j++){
        ficheroEscrituraMath << "{" << rValores[j] << "," << g_L[j] << "}";
        if(j != (rBins-1)){
            ficheroEscrituraMath << ",";
        }
    }
    ficheroEscrituraMath << "}";

    ficheroEscrituraMath.close();

}

void escribirFicheroAnisotropia(int rBins, double rValores[], double thValores[], double F[][179], double rWidth){
    ofstream ficheroEscritura("funcionAnisotropia.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "*********************** FUNCIÓN DE ANISOTROPÍA ************************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;
    ficheroEscritura << setw(10) << left << "r (cm)" << setw(10) << left << "th (º)" << setw(10) << left << "F(r,th)" << endl;
    for(int i = 0; i < 179; i++){
        for(int j = 0; j < rBins; j++){
            ficheroEscritura << left << setw(10) << rValores[j] << left << setw(10) << thValores[i] << left << setw(10) << F[j][i] << endl;
        }
    }

    ficheroEscritura.close();

    ofstream ficheroEscrituraMath1("funcionAnisotropia1.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath1){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath1 << "***********************************************************************" << endl;
    ficheroEscrituraMath1 << "*********************** FUNCIÓN DE ANISOTROPÍA r=1 ************************" << endl;
    ficheroEscrituraMath1 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath1 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath1 << "{" << thValores[i] << "," << F[(int)(1/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath1 << ",";
        }
    }
    ficheroEscrituraMath1 << "}";

    ficheroEscrituraMath1.close();

    ofstream ficheroEscrituraMath2("funcionAnisotropia2.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath2){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath2 << "***********************************************************************" << endl;
    ficheroEscrituraMath2 << "*********************** FUNCIÓN DE ANISOTROPÍA r=5 ************************" << endl;
    ficheroEscrituraMath2 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath2 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath2 << "{" << thValores[i] << "," << F[(int)(5/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath2 << ",";
        }
    }
    ficheroEscrituraMath2 << "}";

    ficheroEscrituraMath2.close();

    ofstream ficheroEscrituraMath3("funcionAnisotropia3.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscrituraMath3){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscrituraMath3 << "***********************************************************************" << endl;
    ficheroEscrituraMath3 << "*********************** FUNCIÓN DE ANISOTROPÍA r=10 ************************" << endl;
    ficheroEscrituraMath3 << "***********************************************************************\n" << endl;

    ficheroEscrituraMath3 << "{";
    for(int i = 0; i < 179; i++){
        ficheroEscrituraMath3 << "{" << thValores[i] << "," << F[(int)(10/rWidth+0.5) - 1][i] << "}";
        if(i != 178){
            ficheroEscrituraMath3 << ",";
        }
    }
    ficheroEscrituraMath3 << "}";

    ficheroEscrituraMath3.close();
}

void escribirFicheroAnisotropia2(double r, double thValores[], double F[][179], double rWidth){
    ofstream ficheroEscritura("funcionAnisotropiaExtra.dat");

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "******************** FUNCIÓN DE ANISOTROPÍA (r=" << r << ") *********************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;
    ficheroEscritura << setw(10) << left << "th (º)" << setw(10) << left << "F(th)" << endl;
    cout << setw(10) << left << "th (º)" << setw(10) << left << "F(th)" << endl;

    for(int i = 0; i < 179; i++){
        ficheroEscritura << setw(10) << left << thValores[i] << setw(10) << left << F[(int)(r/rWidth+0.5) - 1][i] << endl;
        cout << setw(10) << left << thValores[i] << setw(10) << left << F[(int)(r/rWidth+0.5) - 1][i] << endl;
    }

    ficheroEscritura.close();

    cout << "\nSe ha guardado este resultado en el archivo funcionAnisotropia2.dat\n" << endl;
}

void reconstruirDosis(double dosisReconstruida[][179], int rBins, double S_k, double Lambda,
                      double G_L_0, double G_L[][179], double g_L[], double F[][179]){

    for(int i = 0; i < 180; i++){
        for(int j = 0; j < rBins; j++){
            dosisReconstruida[j][i] = S_k * Lambda * G_L[j][i] * g_L[j] * F[j][i] / G_L_0;
        }
    }
}

void escribirFicheroDosisReconstruida(int rBins, double rValores[], double thValores[], double dosisReconstruida[][179]){
    ofstream ficheroEscritura("dosisReconstruida.dat");

    double th;

    // Sale del programa si no puede crear el archivo
    if(!ficheroEscritura){
        cerr << "No se pudo crear el archivo de resultados" << endl;
        exit(1);
    }

    cout << "\n Reconstruir dosis para theta = ";
    cin >> th;

    ficheroEscritura << "***********************************************************************" << endl;
    ficheroEscritura << "******************** DOSIS RECONSTRUIDA (th=" << th << ") *********************" << endl;
    ficheroEscritura << "***********************************************************************\n" << endl;
    ficheroEscritura << setw(10) << left << "r (cm)" << setw(10) << left << "D(r)" << endl;

    for(int i = 0; i < rBins; i++){
        ficheroEscritura << setw(10) << left << rValores[i] << setw(10) << left << dosisReconstruida[i][(int)th - 1] << endl;
    }

    ficheroEscritura.close();
}

double roundit(double num, double N){
    double d = log10(num);
    double power;
    if (num > 0)
    {
        d = ceil(d);
        power = -(d-N);
    }
    else
    {
        d = floor(d);
        power = -(d-N);
    }

    return (int)(num * pow(10.0, power) + 0.5) * pow(10.0, -power);
}
