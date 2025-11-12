#include <cmath>
#include <vector>
#include <iostream>     // Pour std::cout, std::endl
#include <iomanip>      // Pour std::setw
#include <fstream>      // Pour std::ofstream
#include <sstream>      // Pour std::ostringstream
#include <string>
#include "femtool.hpp"  // Assurez-vous que ce fichier est accessible

template <typename PrecondType>
std::vector<double> pcgsolve(const CooMatrix<double>&   A,
                            const PrecondType&          P,
                            const std::vector<double>&  b,
                            const std::vector<double>& u_exact)
{
    
  // Vérification des dimensions
  assert((NbCol(A)==NbRow(A)) && 
        (NbCol(P)==NbRow(P))  &&
        (b.size()==NbCol(A))  && (b.size()==u_exact.size()) );

  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    z   = P*r;
  auto    p   = z;
  auto   Ap   = A*p;
  double rz   = std::real((r|z));
  double r2 = std::pow(Norm(r), 2); // r_0^T * r_0
  double eps2 = 1.e-12 * r2; 
  eps2       *= std::abs((b|b));
  double alpha,beta,pAp;
  std::vector<double> relative_errors;
  double norm_ue = Norm(u_exact);

  // Erreur à l'itération 0 (pour x=0)
  relative_errors.push_back(Norm(x - u_exact) / norm_ue);
    
  std::size_t niter = 0;    

  // Affichage de débogage AVANT la boucle
  std::cout << "--- Début du CG ---" << std::endl;
  std::cout << "Norme résidu initial au carré (r2) : " << r2 << std::endl;
  std::cout << "Seuil de tolérance au carré (eps2) : " << eps2 << std::endl;
  std::cout << "Condition (r2 > eps2) : " 
            << (r2 > eps2 ? "VRAI (entre dans la boucle)" : "FAUX (saute la boucle)") 
            << std::endl;
  std::cout << "-------------------" << std::endl;

  // 2. Boucle du Gradient Conjugué
  while( r2 > eps2 && niter++ < 2000 ){
      
    // Calculs de l'itération k
    auto Ap    = A*p;
    pAp   = std::real((Ap|p));
    
    // Si pAp est nul ou trop petit, on arrête pour éviter une division par zéro
    if (std::abs(pAp) < 1.e-30) {
        std::cerr << "AVERTISSEMENT: Division par zéro évitée (pAp proche de 0). Arrêt." << std::endl;
        break;
    }

    Ap    = A*p;
    pAp   = std::real((Ap|p));
    alpha = rz/pAp;
    beta  = 1/rz;
    x    += alpha*p;
    r    -= alpha*Ap;
    r2    = std::pow(Norm(r),2);
    z     = P*r;
    rz    = std::real((r|z));
    beta *= rz;
    p     = beta*p+z;

    // Stockage de l'erreur
    relative_errors.push_back(Norm(x - u_exact) / norm_ue);

    // Affichage (plus fréquent pour voir ce qui se passe)
    if((niter % 50)==0 || niter == 1){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << "Residual Norm: " << std::sqrt(r2) << "\t";
      std::cout << "Relative Error: " << relative_errors.back() << std::endl;
    }
  }

  std::cout << "--- Fin du PCG ---" << std::endl;
  std::cout << "PCG terminé en " << niter << " itérations." << std::endl;
  std::cout << "Résidu final (sqrt(r2)) : " << std::sqrt(r2) << std::endl;
  
  return relative_errors;
}

// --- Fonction Main ---
int main(int argc,char* argv[]){
  // Vérification simple des arguments
  if (argc != 3) {
    std::cerr << "Erreur: Mauvais nombre d'arguments." << std::endl;
    std::cerr << "Usage: ./mon_programme <fichier_mesh> <fichier_csv_sortie>" << std::endl;
    return 1; // Quitte avec une erreur
  }

  // 1. Récupérer les arguments
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
  auto P = CholeskyPrec(A);
  auto b = A(ue); // b = A * u_exact

  std::vector<double> uh; // Le vecteur solution uh (sera rempli par pcgsolve)
  
  // Appel du solveur
  std::cout << "test execution code" << std::endl;
  auto errors_per_iter = pcgsolve(A, P, b, ue);
 
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

  // Plotter
  // Plot(Vh, uh, "tp1-ex2-output");
  // Plot(Vh, uh-ue, "tp1-ex2-output-erreur");

  // Afficher erreur finale
  if (!errors_per_iter.empty()) {
    std::cout << "Erreur relative finale Norm(u_h - u_ex) / Norm(u_ex) = ";
    // std::cout <<  Norm(uh - ue) / Norm(ue) << "\n";
    std::cout << "Erreur finale (depuis tableau) : " << errors_per_iter.back() << std::endl;
  } else {
    std::cerr << "Erreur: Le tableau des erreurs est vide." << std::endl;
  }

  return 0;
}
