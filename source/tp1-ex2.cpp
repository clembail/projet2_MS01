#include <cmath>
#include <femtool.hpp>

int main(){
  // Instantiation of a 2D domain
  Mesh2D Omega;

  // Loading a 2D mesh
  Read(Omega,"tp1-ex2.mesh");

  // Assembly of a finite element space over Omega
  auto Vh   = FeSpace(Omega);

  // Frequency parameter
  double k = 10.*M_PI;
  
  // Function x = (x1,x2) -> cos(omega*x1)
  auto F    = [&k](const R3& x){return std::cos(k*x[0]);};

  // Evaluating nodal values of f at the degrees of freedom of Vh
  auto ue = Vh(F);

  // Assembly of the finite element matrix of the
  // boundary value problem
  // -Delta u + u = rhs on Omega
  // \partial_n u = 0   on the boundary  
  auto M = Mass(Vh);

  auto A = Stiffness(Vh)+M;

  // Computing a corresponding  right-hand side
  auto b = A(ue);

  // Lazy invertion of A (sparse LU factorization) 
  auto InvA = Inv(A);  

  // Computing A^{-1}b
  auto uh = InvA(b);

  // Plotter la solution calculée u_h
  Plot(Vh, uh, "tp1-ex2-output");
  Plot(Vh, uh-ue, "tp1-ex2-output-erreur");

  // Calculer et afficher l'erreur relative
  std::cout << "Erreur relative Norm(u_h - u_ex) / Norm(u_ex) = ";
  std::cout <<  Norm(uh - ue) / Norm(ue) << "\n";

  return 0;
}
