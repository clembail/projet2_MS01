#include <cmath>
#include <femtool.hpp>

int main(){

    // QUESTION 2

    // // Instantiation of a 2D domain
    // Mesh2D Omega;

    // // Loading a 2D mesh
    // Read(Omega,"tp1-exo1.mesh");

    // // Assembly of a finite element space over Omega
    // auto Vh = FeSpace(Omega);

    // // Frequency parameter
    // double k = 10*M_PI;

    // // Function x = (x1,x2) -> cos(k*(x1+x2))
    // auto F = [&k](const R3& x){return std::cos(k*(x[0]+x[1]));};

    // // Evaluating nodal values of f at the degrees of freedom of Vh
    // auto u = Vh(F);

    // // Plotting F with vizir4
    // Plot(Vh,u,"tp1-exo1-output");

    // QUESTION 3

    // Instantiation of a 3D domain
    Mesh3D Omega;

    // Loading a 3D mesh
    Read(Omega,"tp1-2.mesh");

    // Assembly of a finite element space over Omega
    FeSpace<3> Vh(Omega);

    // Assembly of a finite element space over the boundary of Omega
    auto [Wh,B] = Boundary(Vh);

    // Frequency parameter
    double k = 5*M_PI;

    // Function x = (x1,x2) -> cos(k*(x1+x2+x3))
    auto F = [&k](const R3& x){return std::cos(k*(x[0]+x[1]+x[2]));};

    // Evaluating nodal values of f at the degrees of freedom of Vh
    auto u = Wh(F);

    // Plotting F with vizir4
    Plot(Wh,u,"tp1-2-output");
}