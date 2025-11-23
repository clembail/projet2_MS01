#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "coomatrix.hpp"
#include "element.hpp"
#include "fespace.hpp"
#include "mesh.hpp"
#include <cstddef>
#include <tuple>

using Mesh2DPart = std::vector<Mesh2D>;
using FeSpace2D = FeSpace<2>;
using FeSpace2DxCoo = std::pair<FeSpace2D,CooMatrix<double>>;

std::pair< Mesh2DPart,CooMatrix<double> > 
Partition4(const Mesh2D& Omega){
    auto nodes = Omega.nodes();
    Mesh2D subMesh01 = nodes;
    Mesh2D subMesh02 = nodes;
    Mesh2D subMesh03 = nodes;
    Mesh2D subMesh04 = nodes;
    Mesh2DPart Sigma = {subMesh01, subMesh02, subMesh03, subMesh04};
    int nb_elements = Omega.size();
    CooMatrix<double> Q(nb_elements,4);

    // calcul des coordonnées min/max et du centre
    double minX = 1.e10 ; double maxX = -1.e10 ;
    double minY = 1.e10 ; double maxY = -1.e10 ;

    for(const auto& elt : Omega) {
        for(std::size_t k=0; k<3; ++k) {
            const auto& node = elt[k];
            
            if(node[0] < minX) minX = node[0];
            if(node[0] > maxX) maxX = node[0];
            if(node[1] < minY) minY = node[1];
            if(node[1] > maxY) maxY = node[1];
        }
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

void Plot(const std::vector<Mesh2D>& Sigma,
    std::filesystem::path    filename){
        
    if (Sigma.empty()) return;

    std::vector<std::string> tag =
    {"Vertices", "Edges", "Triangles", "Tetrahedra"};
    
    filename.replace_extension(".mesh");
    std::ofstream f;
    f.open(filename.c_str());

    f << "MeshVersionFormatted 3\n\n";
    f << "Dimension\n3\n\n";
    
    int nb_domains = Sigma.size();
    std::size_t nb_nodes = (Sigma[0].nodes()).size();
    std::size_t nb_elements = 0;
    std::size_t DIM = 2;

    auto v = Sigma[0].nodes();
    auto& v0 = v[0];

    for (std::size_t i = 0; i < nb_domains; i++){
        std::cout << Sigma[i].size() << std::endl;
        nb_elements += Sigma[i].size();
    }

    f << "Vertices\n";
    f << nb_nodes << "\n";
    for(const auto& x:v){
        f << x << "\t1\n";
    }
    f << "\n";

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

    f << "\nEnd";
    f.close();

    filename.replace_extension(".sol");
    std::ofstream f_sol;
    f_sol.open(filename.c_str());

    f_sol << "MeshVersionFormatted 3\n\nDimension\n3\n\n";

    f_sol << "SolAtTriangles\n";
    f_sol << nb_elements << "\n";
    f_sol << "1 1\n";

    for (std::size_t j = 0; j < Sigma.size(); ++j) {
        double partition_value = (double)(j + 1);
        
        for (std::size_t i = 0; i < Sigma[j].size(); ++i) {
            f_sol << partition_value << "\n";
        }
    }

    f_sol << "\nEnd";
    f_sol.close();
  
}

// Surcharge de la fonction Partition4
std::pair<std::vector<Mesh2D>, CooMatrix<double>> 
Partition4(const Mesh2D& Omega, const std::size_t& nl) {

    auto [Sigma, Q] = Partition4(Omega);

    std::vector<Mesh2D> Gamma;
    Gamma.reserve(4);
    CooMatrix<double> R(Omega.size(), 4);

    // Pré-traitement de Q pour récupérer les listes d'éléments initiaux
    std::vector<std::vector<std::size_t>> elements_in_p(4);
    
    for(const auto& [row, col, val] : GetData(Q)) {
        if (val > 0.0) {
            elements_in_p[col].push_back(row);
        }
    }

    // Boucle principale sur les 4 sous-domaines
    for(std::size_t p = 0; p < 4; ++p) {
        
        // Initialisation avec les données de Q
        std::vector<std::size_t> current_elements = elements_in_p[p];

        std::vector<bool> is_in_gamma(Omega.size(), false);
        for(auto idx : current_elements) {
            is_in_gamma[idx] = true;
        }

        // Ajout des couches (Overlap)
        auto nodes = Omega.nodes();
        
        for(std::size_t layer = 0; layer < nl; ++layer) {
            std::vector<bool> active_nodes(nodes.size(), false);
            for(auto elt_idx : current_elements) {
                const auto& e = Omega[elt_idx];
                active_nodes[int(&e[0] - &nodes[0])] = true;
                active_nodes[int(&e[1] - &nodes[0])] = true;
                active_nodes[int(&e[2] - &nodes[0])] = true;
            }

            for(std::size_t k = 0; k < Omega.size(); ++k) {
                if(is_in_gamma[k]) continue; 

                const auto& e = Omega[k];
                if( active_nodes[int(&e[0] - &nodes[0])] || 
                    active_nodes[int(&e[1] - &nodes[0])] || 
                    active_nodes[int(&e[2] - &nodes[0])] ) {
                    
                    is_in_gamma[k] = true;
                    current_elements.push_back(k);
                }
            }
        }

        // Construction de Gamma[p] et remplissage de R
        Mesh2D mesh_p(nodes);
        
        for(auto global_idx : current_elements) {
            mesh_p.push_back(Omega[global_idx]);
            
            double local_idx = (double)(mesh_p.size() - 1);
            R.push_back(global_idx, p, local_idx);
        }
        
        Gamma.push_back(mesh_p);
    }
    
    R.sort();
    return {Gamma, R};
}

std::pair< Mesh2DPart,CooMatrix<double> > 
Partition16(const Mesh2D& Omega){
    auto [Sigma,Q] = Partition4(Omega);
    std::size_t nb_elements = Omega.size();

    Mesh2DPart Sigma_16;
    Sigma_16.reserve(16);

    CooMatrix<double> Q_16(nb_elements,16);
    Q_16.reserve(Omega.size());

    std::vector<std::vector<std::size_t>> global_indices(4);
    for(const auto& [row, col, val] : GetData(Q)) {
        if (val > 0.0) {
            global_indices[col].push_back(row);
        }
    }

    for (std::size_t p = 0; p < 4; p++){
        auto [Sigma_local,Q_local] = Partition4(Sigma[p]);
        Sigma_16.insert(Sigma_16.end(), Sigma_local.begin(), Sigma_local.end());

        const auto& map_p = global_indices[p]; // Le dictionnaire pour p

        for(const auto& [local_row, local_col, val] : GetData(Q_local)) {
            if (val > 0.0) {
                std::size_t global_row = map_p[local_row];
                
                std::size_t final_col = local_col + (p * 4);

                Q_16.push_back(global_row, final_col, 1.0);
            }
        }
    }

    Q_16.sort();

    return {Sigma_16, Q_16};
}

// Surcharge de la fonction Partition16
std::pair<std::vector<Mesh2D>, CooMatrix<double>> 
Partition16(const Mesh2D& Omega, const std::size_t& nl) {

    // Récupérer la partition de base (sans recouvrement)
    auto [Sigma, Q] = Partition16(Omega);

    std::vector<Mesh2D> Gamma;
    Gamma.reserve(16);
    CooMatrix<double> R(Omega.size(), 16);

    // Pré-traitement de Q pour récupérer les listes d'éléments initiaux
    std::vector<std::vector<std::size_t>> elements_in_p(16);
    
    for(const auto& [row, col, val] : GetData(Q)) {
        if (val > 0.0) {
            elements_in_p[col].push_back(row);
        }
    }

    // Boucle principale sur les 4 sous-domaines
    for(std::size_t p = 0; p < 16; ++p) {
        
        // Initialisation avec les données de Q
        std::vector<std::size_t> current_elements = elements_in_p[p];

        std::vector<bool> is_in_gamma(Omega.size(), false);
        for(auto idx : current_elements) {
            is_in_gamma[idx] = true;
        }

        // Ajout des couches (Overlap)
        auto nodes = Omega.nodes();
        
        for(std::size_t layer = 0; layer < nl; ++layer) {
            std::vector<bool> active_nodes(nodes.size(), false);
            for(auto elt_idx : current_elements) {
                const auto& e = Omega[elt_idx];
                active_nodes[int(&e[0] - &nodes[0])] = true;
                active_nodes[int(&e[1] - &nodes[0])] = true;
                active_nodes[int(&e[2] - &nodes[0])] = true;
            }

            for(std::size_t k = 0; k < Omega.size(); ++k) {
                if(is_in_gamma[k]) continue; 

                const auto& e = Omega[k];
                if( active_nodes[int(&e[0] - &nodes[0])] || 
                    active_nodes[int(&e[1] - &nodes[0])] || 
                    active_nodes[int(&e[2] - &nodes[0])] ) {
                    
                    is_in_gamma[k] = true;
                    current_elements.push_back(k);
                }
            }
        }

        // Construction de Gamma[p] et remplissage de R
        Mesh2D mesh_p(nodes);
        
        for(auto global_idx : current_elements) {
            mesh_p.push_back(Omega[global_idx]);
            
            double local_idx = (double)(mesh_p.size() - 1);
            R.push_back(global_idx, p, local_idx);
        }
        
        Gamma.push_back(mesh_p);
    }
    
    R.sort();
    return {Gamma, R};
}

/////////
// TP3 //
/////////

// EXERCICE 1

FeSpace2DxCoo Restrict(const FeSpace2D& Vh, 
    const Mesh2D& Gamma, 
    const std::vector<std::size_t>& tbl){
        auto Uh = FeSpace(Gamma);
        CooMatrix<double> P(dim(Uh),dim(Vh));
        
        for (std::size_t j= 0; j < Gamma.size(); j++){
            std::size_t k = tbl[j];
            
            auto cell_Uh = Uh[j];
            auto cell_Vh = Vh[k];

            for (std::size_t i=0; i < cell_Uh.size(); i++){
                std::size_t row = cell_Uh[i];
                std::size_t col = cell_Vh[i];
                
                P.push_back(row,col,1.0);
            }
        }
        P.sort();
        for (auto& [r,c,v] : GetData(P)){
            v = 1.0;
        }
        return {Uh,P};
    }

std::vector<FeSpace2DxCoo>
Partition4(const FeSpace2D& Vh, const std::size_t& nl){

    auto Omega = Vh.mesh();
    auto [Gamma,R] = Partition4(Omega,nl);

    std::vector<FeSpace2DxCoo> result;
    result.reserve(Gamma.size());

    std::vector<std::vector<std::size_t>> tbl(Gamma.size());

    for (int p = 0; p < Gamma.size();p++){
        tbl[p].resize(Gamma[p].size());
    }

    for (auto& [j,p,vDouble] : GetData(R)){
        std::size_t v = static_cast<std::size_t>(vDouble); //conversion du double en size_t
        tbl[p][v] = j;
    }

    for (int p=0; p < Gamma.size(); p++){
        result.push_back(Restrict(Vh, Gamma[p], tbl[p]));
    }

    return result;
}

std::vector<FeSpace2DxCoo>
Partition16(const FeSpace2D& Vh, const std::size_t& nl){

    auto Omega = Vh.mesh();
    auto [Gamma,R] = Partition16(Omega,nl);

    std::vector<FeSpace2DxCoo> result;
    result.reserve(Gamma.size());

    std::vector<std::vector<std::size_t>> tbl(Gamma.size());

    for (int p = 0; p < Gamma.size();p++){
        tbl[p].resize(Gamma[p].size());
    }

    for (auto& [j,p,vDouble] : GetData(R)){
        std::size_t v = static_cast<std::size_t>(vDouble); //conversion du double en size_t
        tbl[p][v] = j;
    }

    for (int p=0; p < Gamma.size(); p++){
        result.push_back(Restrict(Vh, Gamma[p], tbl[p]));
    }

    return result;
}

#endif