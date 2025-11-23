#include <cmath>
#include <vector>
#include <iostream>     // Pour std::cout, std::endl
#include <iomanip>      // Pour std::setw
#include <fstream>      // Pour std::ofstream
#include <sstream>      // Pour std::ostringstream
#include <string>
#include "femtool.hpp"  // Assurez-vous que ce fichier est accessible


int main(int argc,char* argv[]){

  // QUESTION 2 ET 3

  //// Instantiation of a 2D domain
  // Mesh2D Omega;

  // // Loading a 2D mesh
  // Read(Omega,"tp1-ex2.mesh");

  // // Assembly of a finite element space over Omega
  // auto Vh   = FeSpace(Omega);

  // // Frequency parameter
  // double k = 10.*M_PI;
  
  // // Function x = (x1,x2) -> cos(omega*x1)
  // auto F    = [&k](const R3& x){return std::cos(k*x[0]);};

  // // Evaluating nodal values of f at the degrees of freedom of Vh
  // auto ue = Vh(F);

  // // Assembly of the finite element matrix of the
  // // boundary value problem
  // // -Delta u + u = rhs on Omega
  // // \partial_n u = 0   on the boundary  
  // auto M = Mass(Vh);

  // auto A = Stiffness(Vh)+M;

  // // Computing a corresponding  right-hand side
  // auto b = A(ue);

  // // Lazy invertion of A (sparse LU factorization) 
  // auto InvA = Inv(A);  

  // // Computing A^{-1}b
  // auto uh = InvA(b);

  // // Plotter la solution calculée u_h
  // Plot(Vh, uh, "tp1-ex2-output");
  // Plot(Vh, uh-ue, "tp1-ex2-output-erreur");

  // // Calculer et afficher l'erreur relative
  // std::cout << "Erreur relative Norm(u_h - u_ex) / Norm(u_ex) = ";
  // std::cout <<  Norm(uh - ue) / Norm(ue) << "\n";

  // QUESTION 4 À 6

  if (argc != 3) {
    std::cerr << "Erreur: Mauvais nombre d'arguments." << std::endl;
    std::cerr << "Usage: ./mon_programme <fichier_mesh> <fichier_csv_sortie>" << std::endl;
    return 1;
  }

  std::string mesh_filename = argv[1];
  std::string csv_filename  = argv[2];

  // Instantiation of a 2D domain
  Mesh2D Omega;
  Read(Omega,mesh_filename);
  auto Vh   = FeSpace(Omega);
  std::cout << Omega.size() << std::endl;

  // Frequency parameter
  double k = 10.*M_PI; // k élevé pour tester
  
  auto F    = [&k](const R3& x){return std::cos(k*x[0]);};
  auto ue = Vh(F); // ue = Solution exacte (discrète)
  
  auto A = Stiffness(Vh)+Mass(Vh);
  auto b = A(ue); // b = A * u_exact

  std::vector<double> uh; // Le vecteur solution uh (sera rempli par cgsolve)
  
  // Appel du solveur
  std::cout << "test execution code" << std::endl;
  auto errors_per_iter = cgsolve(A, b, ue, uh);
 
  std::cout << "===================================================" << std::endl;
  std::cout << "La méthode a convergé en " << errors_per_iter.size() - 1 << " itérations." << std::endl;

  std::ofstream outfile(csv_filename);
  if (outfile.is_open()) {
    outfile << "iteration,relative_error\n";
    for (size_t i = 0; i < errors_per_iter.size(); ++i) {
      outfile << i << "," << errors_per_iter[i] << "\n";
    }
    outfile.close();
    std::cout << "Résultats écrits dans : " << csv_filename << std::endl;
  } else {
    std::cerr << "Erreur: Impossible d'ouvrir le fichier " << csv_filename << std::endl;
  }
  std::cout << "===================================================" << std::endl;

  // Afficher erreur finale
  if (!errors_per_iter.empty()) {
    std::cout << "Erreur relative finale Norm(u_h - u_ex) / Norm(u_ex) = ";
    std::cout <<  Norm(uh - ue) / Norm(ue) << "\n";
    std::cout << "Erreur finale (depuis tableau) : " << errors_per_iter.back() << std::endl;
  } else {
    std::cerr << "Erreur: Le tableau des erreurs est vide." << std::endl;
  }

  return 0;
}
