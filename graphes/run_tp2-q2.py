import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import os

# 1. Définir les maillages et le k
# (Vous devez avoir généré ces maillages avec gmsh !)
meshes = {
    0.025: "../solutions/mesh_h_0.025.mesh",
    0.01: "../solutions/mesh_h_0.01.mesh",
    0.005: "../solutions/mesh_h_0.005.mesh",
    0.0025: "../solutions/mesh_h_0.0025.mesh"
}
cpp_executable = "../build/tp2-ex1" # Nom de votre C++ compilé

plt.figure(figsize=(10, 7))

# 2. Boucle pour exécuter le C++
for h, mesh_file in meshes.items():
    csv_file = f"results_h_{h}_PCG.csv"
    print(f"--- Lancement du solveur C++ pour h={h} ---")
    
    # Construit la commande
    command = [cpp_executable, mesh_file, csv_file]
    
    # Exécute la commande
    subprocess.run(command, check=True)
    
    # 3. Lire le CSV et tracer (répond à Q4 et Q5)
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        plt.plot(df['iteration'], df['relative_error'], label=f'h = {h}')

print("--- Tracé des résultats (Q2) ---")
plt.yscale('log')
plt.title("Convergence du CG (Erreur H1) pour différents maillages")
plt.xlabel("Itération (k)")
plt.ylabel("Erreur relative $||u_k - u_h||_{H1} / ||u_h||_{H1}$")
plt.legend()
plt.grid(True, which="both")
plt.savefig("tp2_ex1_question_2_plot.png")
print("Graphique 'tp2_ex1_question_2_plot.png' sauvegardé.")
