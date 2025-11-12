#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "coomatrix.hpp"
#include "mesh.hpp"

using Mesh2DPart = std::vector<Mesh2D>;

std::pair< Mesh2DPart,CooMatrix<double> > Partition4(const Mesh2D& Omega);

void Plot(const std::vector<Mesh2D>& Sigma,
    const std::string& filename){

  assert( mesh.size()==u.size() );

  std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};
  
  //###############//
  //  Fichier mesh //  
  Write(mesh,filename);

  //###############//
  //  Fichier sol  //

  // Ouverture
  filename.replace_extension(".sol");
  std::ofstream f;
  f.open(filename.c_str());

  // Preambule
  f << "MeshVersionFormatted 3\n\n";
  f << "Dimension\n3\n\n";
  
  // Donnees
  f << "SolAt";
  f << tag[DIM] << "\n";
  f << mesh.size() << "\n";
  f << "1\t1\n";
  for(std::size_t j=0; j<u.size(); ++j){
    f << u[j] << "\n";}

  // Fermeture
  f << "\nEnd";
  f.close();
  
}

#endif