
#include <math.h>
#include <iostream>
#include <iomanip>

#include "GaussLegendreQuadrature.h"


int main() {
    Rosetta::GaussLegendreQuadrature<5> gl5;
 
    std::cout << std::setprecision(10);
 
    gl5.print_roots_and_weights(std::cout);

    double rho = 8.18;
    double delta = 0.14;
    std::cout << "Integrating Exp(X) over [0, 50]: " << gl5.integrate(0., 50., Rosetta::RosettaExp, delta, 0.1) << '\n';
    std::cout << "Actual value:                    " << Rosetta::RosettaExp(delta,50,0.1) - Rosetta::RosettaExp(delta,0,0.1) << '\n'; 

}