#include <iostream>
#include <vector>
#include "pbPlots/pbPlots.cpp"
#include "pbPlots/supportLib.cpp"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

double f(double x){
    return x*x*x + x;
}

MatrixXd macierzA(long long n, double h){
    MatrixXd matrix(n,n);

    for( int i = 0; i < n; i++ )
    {
        for( int j = 0; j < n; j++){
            matrix(i,j) = 0.0;
        }
    }

    matrix(0,0) = 2.0/pow(h,2)+2.0;
    matrix(0,1) = -2.0/pow(h,2);

    for( int i = 1; i < n-1; i++){
        matrix(i,i-1) = -1.0/pow(h,2)-1.0/(2.0*h);
        matrix(i,i) = 2.0/pow(h,2)+2.0;
        matrix(i,i+1) = -1.0/pow(h,2)+1.0/(2.0*h);
    }

    matrix(n-1,n-2) = -2.0/pow(h,2);
    matrix(n-1,n-1) = 2.0/pow(h,2)+2.0/h+1.0;

    return matrix.inverse();
}

MatrixXd macierzF(long long n, double h, vector<double> t){
    MatrixXd matrix(n,1);
    matrix(0,0) = -1.0-(2.0/h);

    for( int i = 1; i < n; i++ ){
        matrix(i, 0) = f(t[i]);
    }

    return matrix;
}

vector<double> rozw_dokladne(long long n, vector<double>t){
    vector<double> u;

    for(int i = 0; i < n; i++){
        u.push_back((-25.0/24.0/exp(2.0))*exp(2*t[i])+(-25.0/12.0/exp(2.0)+7.0/4.0)*exp(-t[i])+ pow(t[i],3)/2.0-3.0* pow(t[i],2)/4.0+11.0*t[i]/4.0-17.0/8.0);
    }

    return u;
}

vector<double> rozw_przyblizone(long long n, double h, vector<double> t){
    MatrixXd matrix1 = macierzA(n, h);
    MatrixXd matrix2 = macierzF(n, h, t);
    MatrixXd matrix3 = matrix1*matrix2;
    vector<double> wynik;

    for(int i = 0; i < n; i++ ){
        wynik.push_back(matrix3(i,0));
        //cout << wynik[i] << endl;
    }

    return wynik;
}

double blad_globalny(vector<double> rozw_dokladne, vector<double> rozw_przyblizone){
    long n = rozw_przyblizone.size();
    double blad_gl = 0.0;

    for( int i = 0; i < n; i++ ){
        if(abs(rozw_przyblizone[i] - rozw_dokladne[i]) > blad_gl){
            blad_gl = abs(rozw_przyblizone[i] - rozw_dokladne[i]);
        }
    }

    return blad_gl;
}

void wykres(vector<double> rozw_dokladne, vector<double> rozw_przyblizone, vector<double> t, string name){
    bool success;
    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();

    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();

    series->xs = &t;
    series->ys = &rozw_przyblizone;
    series->linearInterpolation = true;
    series->lineType = toVector(L"dashed");
    series->lineThickness = 2;
    series->color = CreateRGBColor(1, 0, 0);

    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();

    series2->xs = &t;
    series2->ys = &rozw_dokladne;
    series2->linearInterpolation = true;
    series2->lineType = toVector(L"dashed");
    series2->lineThickness = 2;
    series2->color = CreateRGBColor(0, 1, 0);

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
    settings->width = 800;
    settings->height = 480;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = toVector(L"Zielony - wynik dokladny, Czerwony - wynik przyblizony");
    settings->xLabel = toVector(L"x");
    settings->yLabel = toVector(L"u(x)");
    settings->scatterPlotSeries->push_back(series);
    settings->scatterPlotSeries->push_back(series2);

    success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

    if(success){
        vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, name);
        DeleteImage(imageReference->image);
    }else{
        cerr << "Error: ";
        for(wchar_t c : *errorMessage->string){
            wcerr << c;
        }
        cerr << endl;
    }

}

int main(){

    long long n1;
    long long n2;
    vector <double> t1;
    vector <double> t2;
    double h1;
    double h2;
    vector <double> dokladne_n_10;
    vector <double> przyblizone_n_10;
    vector <double> dokladne_n_50;
    vector <double> przyblizone_n_50;
    double blad_gl_n_10;
    double blad_gl_n_50;

    cout << "Rownanie -u'' + u' + 2u = x^3 + x, u'(0) = 1, u(1)+u'(1) = 0 rozwiazne metoda MRS" << endl;

    for(int i = 0; i < 11; i++){
        t1.push_back(0.0 + (1.0-0.0)*(double)i/10.0);
    }

    n1 = 11;
    h1 = (t1[9]-t1[0])/11.0;

    dokladne_n_10 = rozw_dokladne(n1, t1);
    przyblizone_n_10 = rozw_przyblizone(n1, h1, t1);
    blad_gl_n_10 = blad_globalny(dokladne_n_10, przyblizone_n_10);
    wykres(dokladne_n_10, przyblizone_n_10, t1, "mrs_n_10.png");

    cout << "Wykres dla 10 punktow znajduje sie w pliku mrs_n_10.png" << endl;
    cout << "Blad globalny dla n = 10 wynosi: " << blad_gl_n_10 << endl;

    for(int i = 0; i < 51; i++){
        t2.push_back(0.0 + (1.0-0.0)*(double)i/50.0);
    }

    n2 = 51;
    h2 = (t2[49]-t2[0])/51.0;

    dokladne_n_50 = rozw_dokladne(n2, t2);
    przyblizone_n_50 = rozw_przyblizone(n2, h2, t2);
    blad_gl_n_50 = blad_globalny(dokladne_n_50, przyblizone_n_50);
    wykres(dokladne_n_50, przyblizone_n_50, t2, "mrs_n_50.png");

    cout << "Wykres dla 50 punktow znajduje sie w pliku mrs_n_50.png" << endl;
    cout << "Blad globalny dla n = 50 wynosi: " << blad_gl_n_50 << endl << endl;

    cout << "Nacisnij enter by zakonczyc program" << endl;

    cin.get();

    return 0;
}
