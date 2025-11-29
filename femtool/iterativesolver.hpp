#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>
#include "smallvector.hpp"
#include "coomatrix.hpp"

std::vector<double>
cgsolve(const CooMatrix<double>&   A,
	const std::vector<double>& b) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    p   = r;
  auto   Ap   = A*p;    
  double r2   = std::pow(Norm(r),2);  
  double eps2 = (1e-8)*r2;
  eps2       *= std::abs((b|b));      
  double alpha,beta,pAp;
    
  std::size_t niter = 0;    
  while( r2>eps2 && niter++<1000 ){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));    
    alpha = r2/pAp;
    x    += alpha*p;
    r    -= alpha*Ap;    
    r2    = std::pow(Norm(r),2);
    beta  = r2/(alpha*pAp);
    p     = beta*p+r;
      
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }
      
  }
    
  return x;
    
}

std::vector<double> cgsolve(const CooMatrix<double>&   A,
                            const std::vector<double>& b,
                            const std::vector<double>& u_exact,
                            std::vector<double>&       x)
{
    
  assert((NbCol(A)==NbRow(A)) && 
        (b.size()==NbCol(A)) && (b.size()==u_exact.size()) );
 
  x.assign(b.size(), 0.0);
  
  auto    r   = b;
  auto    p   = r;
  
  double r2_old = std::pow(Norm(r), 2);
  double r2_new = r2_old;
  
  double eps2 = 1.e-12 * r2_old; 
  
  double alpha, beta, pAp;
  std::vector<double> relative_errors;
  double norm_ue = Norm(u_exact);

  relative_errors.push_back(Norm(x - u_exact) / norm_ue);
    
  std::size_t niter = 0;    

  // Affichage de débogage AVANT la boucle
  std::cout << "--- Début du CG ---" << std::endl;
  std::cout << "Norme résidu initial au carré (r2_new) : " << r2_new << std::endl;
  std::cout << "Seuil de tolérance au carré (eps2)     : " << eps2 << std::endl;
  std::cout << "Condition (r2_new > eps2)            : " 
            << (r2_new > eps2 ? "VRAI (entre dans la boucle)" : "FAUX (saute la boucle)") 
            << std::endl;
  std::cout << "-------------------" << std::endl;

  while( r2_new > eps2 && niter++ < 2000 ){
      
    auto Ap    = A*p;
    pAp   = std::real((Ap|p));
    
    // Si pAp est nul ou trop petit, on arrête pour éviter une division par zéro
    if (std::abs(pAp) < 1.e-30) {
        std::cerr << "AVERTISSEMENT: Division par zéro évitée (pAp proche de 0). Arrêt." << std::endl;
        break;
    }
        
    alpha = r2_old / pAp;
    
    x    += alpha*p;
    r    -= alpha*Ap;
    
    r2_new = std::pow(Norm(r), 2);
    
    beta   = r2_new / r2_old;
    p      = r + beta*p;
      
    r2_old = r2_new;

    relative_errors.push_back(Norm(x - u_exact) / norm_ue);

    if((niter % 50)==0 || niter == 1){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << "Residual Norm: " << std::sqrt(r2_new) << "\t";
      std::cout << "Relative Error: " << relative_errors.back() << std::endl;
    }
  }

  std::cout << "--- Fin du CG ---" << std::endl;
  std::cout << "CG terminé en " << niter << " itérations." << std::endl;
  std::cout << "Résidu final (sqrt(r2_new)) : " << std::sqrt(r2_new) << std::endl;
  
  return relative_errors;
}


template <typename PrecondType>
std::vector<double>
pcgsolve(const CooMatrix<double>&   A,
  const PrecondType&     P,
	const std::vector<double>& b) {
    
  assert((NbCol(A)==NbRow(A)) &&
	 (b.size()==NbCol(A)) );
 
  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    z   = P*r;
  auto    p   = z;
  auto   Ap   = A*p;
  double rz   = std::real((r|z));
  double r2   = std::pow(Norm(r),2);
  double eps2 = (1e-8)*r2;
  eps2       *= std::abs((b|b));
  double alpha,beta,pAp;
    
  std::size_t niter = 0;
  while( r2>eps2 && niter++<1000 ){
      
    Ap    = A*p;
    pAp   = std::real((Ap|p));
    alpha = rz/pAp;
    beta  = 1/rz;
    x    += alpha*p;
    r    -= alpha*Ap;
    z     = P*r;
    r2    = std::pow(Norm(r),2);
    beta *= rz;
    p     = beta*p+r;
      
    if((niter%50)==0){
      std::cout << std::left << std::setw(7) << niter << "\t";
      std::cout << Norm(r) << std::endl;
    }

  }
    
  return x;
    
}

template <typename PrecondType>
std::vector<double> pcgsolve(const CooMatrix<double>&   A,
                            const PrecondType&          P,
                            const std::vector<double>&  b,
                            const std::vector<double>& u_exact)
{
    
  assert((NbCol(A)==NbRow(A)) && 
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

  while( r2 > eps2 && niter++ < 2000 ){
      
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

    relative_errors.push_back(Norm(x - u_exact) / norm_ue);

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

template <typename PrecondType>
std::vector<double> pcgsolve_history(const CooMatrix<double>& A,
                                     const PrecondType& P,
                                     const std::vector<double>& b,
                                     double tol = 1e-10,
                                     int max_iter = 1000) {
  assert((NbCol(A)==NbRow(A)) && 
        (b.size()==NbCol(A)) );

  auto    x   = std::vector<double>(b.size(),0.);
  auto    r   = b-A*x;
  auto    z   = P*r;
  auto    p   = z;
  auto   Ap   = A*p;
  double rz   = std::real((r|z));
  double r2 = std::pow(Norm(r), 2); // r_0^T * r_0
  double eps2 = tol * tol * r2; 
  eps2       *= std::abs((b|b));
  double alpha,beta,pAp;
  std::vector<double> relative_errors;

  relative_errors.push_back(Norm(r) / Norm(b));
    
  std::size_t niter = 0;    

  // Affichage de débogage AVANT la boucle
  std::cout << "--- Début du CG ---" << std::endl;
  std::cout << "Norme résidu initial au carré (r2) : " << r2 << std::endl;
  std::cout << "Seuil de tolérance au carré (eps2) : " << eps2 << std::endl;
  std::cout << "Condition (r2 > eps2) : " 
            << (r2 > eps2 ? "VRAI (entre dans la boucle)" : "FAUX (saute la boucle)") 
            << std::endl;
  std::cout << "-------------------" << std::endl;

  while( r2 > eps2 && niter++ < max_iter ){
      
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

    relative_errors.push_back(Norm(r) / Norm(b));

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

#endif
