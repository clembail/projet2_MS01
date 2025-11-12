#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include <functional>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>

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


#endif
