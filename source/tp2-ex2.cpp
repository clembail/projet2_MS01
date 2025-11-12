#include <iostream>
#include "femtool.hpp"


int main(int argc,char* argv[]){

    std::string input = "./solutions/mesh_h_0.01.mesh";
    std::string output = "./solutions/mesh_h_0.01_partition4.mesh";

    Mesh2D Omega;
    Read(Omega, input);

    std::cout << "Partitionnement du maillage" << std::endl;
    auto [Sigma, Q] = Partition4(Omega);

    std::cout << "Création du maillage prenant en compte le partionnement" << std::endl;
    Plot(Sigma, output);

    return 0;
}