#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "coomatrix.hpp"
#include "element.hpp"
#include "mesh.hpp"
#include <cstddef>

using Mesh2DPart = std::vector<Mesh2D>;

std::pair< Mesh2DPart,CooMatrix<double> > Partition4(const Mesh2D& Omega){
    auto nodes = Omega.nodes();
    Mesh2D subMesh01 = nodes;
    Mesh2D subMesh02 = nodes;
    Mesh2D subMesh03 = nodes;
    Mesh2D subMesh04 = nodes;
    Mesh2DPart Sigma = {subMesh01, subMesh02, subMesh03, subMesh04};
    int nb_elements = Omega.size();
    CooMatrix<double> Q(nb_elements,4);

    // calcul des coordonnées min/max et du centre
    double minX = 0 ; double maxX = 0 ;
    double minY = 0 ; double maxY = 0 ;

    for (const auto& node : nodes){
        if(node[0] < minX) minX = node[0];
        if(node[0] > maxX) maxX = node[0];
        if(node[1] < minY) minY = node[1];
        if(node[1] > maxY) maxY = node[1];
    }

    double centreX = (minX + maxX)/2;
    double centreY = (minY + maxY)/2;

    // calcul des éléments dans chaque sous-maillage en calculant leurs barycentres
    // bas-gauche : 0 ; bas-droite : 1 ; haut-gauche : 2 ; haut-droit : 3
    for (std::size_t j = 0; j < nb_elements ; j++){
        const Element<2>& elt = Omega[j];
        double ex = (elt[0][0] + elt[1][0] + elt[2][0])/3;
        double ey = (elt[0][1] + elt[1][1] + elt[2][1])/3;

        std::size_t p = 0;
        if (ex > centreX) p += 1;
        if (ey > centreY) p += 2;
        
        Sigma[p].push_back(elt);
        Q.push_back(j,p,1.0);
    }
            
    Q.sort();

    return {Sigma, Q};
};

std::pair< Mesh2DPart,CooMatrix<double> >
Partition4(const Mesh2D& Omega, const std::size_t& nl){
    Mesh2DPart Sigma;
    CooMatrix<double> R;
    return {Sigma, R}
};

void Plot(const std::vector<Mesh2D>& Sigma,
    std::filesystem::path    filename){
        
    if (Sigma.empty()) return;

    std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};

    // Ouverture fichier
    filename.replace_extension(".mesh");
    std::ofstream f;
    f.open(filename.c_str());

    // Préambule
    f << "MeshVersionFormatted 3\n\n";
    f << "Dimension\n3\n\n";
    
    int nb_domains = Sigma.size();
    std::size_t nb_nodes = (Sigma[0].nodes()).size();
    std::size_t nb_elements = 0;
    std::size_t DIM = 2;

    auto v = Sigma[0].nodes();
    auto& v0 = v[0];

    for (std::size_t i = 0; i < nb_domains; i++){
        // nb_nodes += (Sigma[i].nodes()).size();
        std::cout << Sigma[i].size() << std::endl;
        nb_elements += Sigma[i].size();
    }

    // Section noeuds
    // auto v = m.nodes();
    // const auto& v0 = v[0];
    f << "Vertices\n";
    f << nb_nodes << "\n";
    for(const auto& x:v){
        f << x << "\t1\n";
    }
    f << "\n";
    // for (std::size_t i = 0; i < nb_domains; i++){   
    //     auto v = Sigma[i].nodes();
    //     for(const auto& x:v){
    //         f << x << "\t1\n"; 
    //     }
    // }
    f << "\n";

    // Section elements
    f << tag[DIM] << "\n";
    f << nb_elements << "\n";
    for (std::size_t j = 0; j < nb_domains; j++){
        auto m = Sigma[j];
        for(const auto& e:m){
            for(std::size_t k=0; k<DIM+1; k++){
                f << 1 + int(&e[k]-&v0) << "\t";
            }
            f << (j+1) << "\n";
        }
    }

    // Fermeture
    f << "\nEnd";
    f.close();
  
}

#endif