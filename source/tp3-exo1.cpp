#include <iostream>
#include <cmath>
#include <string>
#include <femtool.hpp>

int main(){
    std::string input = "./tp/solutions/mesh_h_0.01_partition4.mesh";

    Mesh2D Omega;
    Read(Omega, input);

    std::cout << "Partitionnement du maillage en 4" << std::endl;
    auto [Sigma, R] = Partition4(Omega,0);

    std::vector<std::size_t> tbl(Sigma[0].size());

    for (auto& [j,p,vDouble] : GetData(R)){
        if (p == 0){
            std::size_t v = static_cast<std::size_t>(vDouble); //conversion du double en size_t
            tbl[v] = j;
        }
    }

    auto Vh = FeSpace(Omega);
    auto [Uh,P] = Restrict(Vh,Sigma[0], tbl);

    double k = 7*M_PI;
    auto F = [&k](const R3& x){return std::cos(k*(x[0]+x[1]));};

    auto fh = Vh(F);

    auto fhRestricted = P*fh;

    Plot(Uh,fhRestricted, "./tp/solutions/tp3-exo1");

    return 0;
}