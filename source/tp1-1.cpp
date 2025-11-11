#include <cmath>
#include <femtool.hpp>

int main(){
    // Instantiation of a 2D domain
    Mesh2D Omega;

    // Loading a 2D mesh
    Read(Omega,"tp1-1.mesh");

    // Assembly of a finite element space over Omega
    auto Vh = FeSpace(Omega);

    // Frequency parameter
    double k = 10*M_PI;

    // Function x = (x1,x2) -> cos(k*(x1+x2))
    auto F = [&k](const R3& x){return std::cos(k*(x[0]+x[1]));};

    // Evaluating nodal values of f at the degrees of freedom of Vh
    auto u = Vh(F);

    // Plotting F with vizir4
    Plot(Vh,u,"tp1-1-output");
}