#include <iostream>
#include "femtool.hpp"


int main(int argc,char* argv[]){

    std::string input = "./solutions/mesh_h_0.01.mesh";
    std::string output4 = "./solutions/mesh_h_0.01_partition4.mesh";
    std::string output4o = "./solutions/mesh_h_0.01_partition4o.mesh";
    std::string output16 = "./solutions/mesh_h_0.01_partition16.mesh";
    std::string output16o = "./solutions/mesh_h_0.01_partition16o.mesh";

    Mesh2D Omega;
    Read(Omega, input);

    ////////////////////////////
    // Test partitionnement en 4
    ////////////////////////////
    std::cout << "=========================================================================" << std::endl;

    std::cout << "Partitionnement du maillage en 4" << std::endl;
    auto [Sigma4, Q4] = Partition4(Omega);

    std::cout << "Création du maillage prenant en compte le partionnement" << std::endl;
    Plot(Sigma4, output4);

    //////////////////////////////////////////////
    // TEST partitionnement en 4 avec recouvrement
    //////////////////////////////////////////////
    std::cout << "=========================================================================" << std::endl;

    std::cout << "Partitionnement du maillage en 4 avec recouvrement (nl==2)" << std::endl;
    auto [Gamma4, R4] = Partition4(Omega, 2);

    std::cout << "Création du maillage prenant en compte le partionnement" << std::endl;
    Plot(Gamma4, output4o);

    // 2.1 Vérification Géométrique
    // Comme ça se chevauche, un Plot unique serait illisible.
    // On sauvegarde chaque sous-domaine étendu séparément.
    for(int p=0; p<4; ++p) {
        std::string name = "./solutions/partition/gamma4_" + std::to_string(p);
        Write(Gamma4[p], name); // Utilise la fonction Write de mesh.hpp
        std::cout << "-> " << name << ".mesh généré (" << Gamma4[p].size() << " elts)" << std::endl;
    }

    std::cout << "Matrice R : " << NbRow(R4) << "x" << NbCol(R4) << std::endl;

    // Le nombre de non-zéros dans R doit être égal à la somme des tailles des maillages Gamma
    std::size_t total_size_gamma4 = 0;
    for(const auto& m : Gamma4) total_size_gamma4 += m.size();
    std::cout << "Somme des tailles Gamma : " << total_size_gamma4 << std::endl;
    std::cout << "Non-Zéros dans R : " << Nnz(R4) << std::endl;

    assert(total_size_gamma4 == Nnz(R4)); // Doit être strictement égal
    assert(total_size_gamma4 > Omega.size()); // Doit être plus grand car recouvrement

    /////////////////////////////
    // TEST partitionnement en 16
    /////////////////////////////
    std::cout << "=========================================================================" << std::endl;

    std::cout << "Partitionnement du maillage en 16" << std::endl;
    auto [Sigma16, Q16] = Partition16(Omega);

    std::cout << "Création du maillage prenant en compte le partionnement" << std::endl;
    Plot(Sigma16, output16);

    std::size_t total_elements_16 = 0;
    for(const auto& m : Sigma16) total_elements_16 += m.size();

    std::cout << "Elements total Omega : " << Omega.size() << std::endl;
    std::cout << "Somme elements Sigma16 : " << total_elements_16 << std::endl;
    assert(Omega.size() == total_elements_16);

    auto sum_rows = Q16 * std::vector<double>(16, 1.0); // Somme des colonnes
    double min_sum = 1e9, max_sum = -1e9;
    for(auto v : sum_rows) {
        if(v < min_sum) min_sum = v;
        if(v > max_sum) max_sum = v;
    }
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Q16 Row Sum -> Min: " << min_sum << " Max: " << max_sum << std::endl;
    assert(std::abs(min_sum - 1.0) < 1e-9 && std::abs(max_sum - 1.0) < 1e-9);

    ///////////////////////////////////////////////
    // TEST partitionnement en 16 avec recouvrement
    ///////////////////////////////////////////////
    std::cout << "=========================================================================" << std::endl;

    std::cout << "Partitionnement du maillage en 16 avec recouvrement (nl == 2)" << std::endl;
    auto [Gamma16, R16] = Partition16(Omega, 2);

    std::cout << "Création du maillage prenant en compte le partionnement" << std::endl;
    Plot(Gamma16, output16o);

    for(int p=0; p<16; ++p) {
        std::string name = "./solutions/partition/gamma16_" + std::to_string(p);
        Write(Gamma16[p], name); // Utilise la fonction Write de mesh.hpp
        std::cout << "-> " << name << ".mesh généré (" << Gamma16[p].size() << " elts)" << std::endl;
    }

    std::cout << "Matrice R : " << NbRow(R16) << "x" << NbCol(R16) << std::endl;

    // Le nombre de non-zéros dans R doit être égal à la somme des tailles des maillages Gamma
    std::size_t total_size_gamma16 = 0;
    for(const auto& m : Gamma16) total_size_gamma16 += m.size();
    std::cout << "Somme des tailles Gamma : " << total_size_gamma16 << std::endl;
    std::cout << "Non-Zéros dans R : " << Nnz(R16) << std::endl;

    assert(total_size_gamma16 == Nnz(R16)); // Doit être strictement égal
    assert(total_size_gamma16 > Omega.size()); // Doit être plus grand car recouvrement

    return 0;
}