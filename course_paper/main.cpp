#include <iostream>
#include <cmath>
#include <iomanip>
#include "FORSYTHE.H"

double alpha = 0.0;

void firstFunction(double t, double* y, double* dX);
double secondFunction(double value);
double auxiliaryFunction(double y);

int main() {
    double ax = 0.7;
    double bx = 1.5;
    alpha = ZEROIN(ax, bx, secondFunction, 0.0000001);
    std::cout << "alpha:" << alpha << '\n';

    double t = 0.0;
    double tOut = 0.0;
    double h = 0.1;
    double X[2] = {0.0,sqrt(2 * alpha)};
    double work[3 + 6 * 1];
    double relerr = 0.0000001;
    double abserr = 0.0000001;
    int flag = 1;

    std::cout << "---------------RKF45---------------\n";
    for (int i = 0; i <= 1 / h; i++)
    {
        RKF45(firstFunction, 1, X, t, tOut, relerr, abserr, work, flag);
        std::cout << tOut << std::setw(15) << X[0] << std::setw(15) << X[1] << '\n';
        tOut += h;
    }
}

void firstFunction(double t, double* y, double* dX) {
    dX[0] = sqrt(2 * ((alpha + pow(y[0], 3) / 3 - y[0])));
}

double secondFunction(double value) {
    alpha = value;
    double result, errest, flag;
    int nofun;
    QUANC8(auxiliaryFunction, 0.0, 1.0, 0.0000001, 0.0000001, result, errest, nofun, flag);
    return result - 1;
}

double auxiliaryFunction(double y) {
    return 1 / sqrt(2 * (alpha + pow(y, 3) / 3 - y));
}
