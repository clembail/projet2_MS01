#include <cmath>
#include "femtool.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mesh_file> <output_csv>" << std::endl;
        return 1;
    }

    std::string mesh_file = argv[1];
    std::string output_csv = argv[2];

    std::cout << "Lecture du maillage : " << mesh_file << std::endl;
    Mesh2D Omega;
    Read(Omega, mesh_file);
    FeSpace2D Vh(Omega);
    
    auto F = [](const R3& x) { return std::cos(10.0 * M_PI * x[0]); };
    auto ue = Vh(F);

    auto A = Stiffness(Vh) + Mass(Vh);
    auto b = A * ue;

    std::size_t nl = 2; 
    std::cout << "Partitionnement en 4 sous-domaines avec " << nl << " couches..." << std::endl;
    
    std::vector<FeSpace2DxCoo> subdomains = Partition4(Vh, nl);


    std::vector<CooMatrix<double>> vec_R;
    vec_R.reserve(subdomains.size());

    for(const auto& domain_pair : subdomains) {
        vec_R.push_back(domain_pair.second);
    }

    std::cout << "Construction du Preconditionneur Schwarz Additif..." << std::endl;
    AdditiveSchwarz ASM(A, vec_R);

    std::cout << "Lancement du PCG..." << std::endl;
    auto residuals = pcgsolve_history(A, ASM, b, 1e-10, 1000);

    std::ofstream f(output_csv);
    f << "iteration,rel_residual\n";
    for(size_t k=0; k<residuals.size(); ++k) {
        f << k << "," << residuals[k] << "\n";
    }
    f.close();
    std::cout << "Resultats sauvegardes dans " << output_csv << std::endl;

    return 0;
}