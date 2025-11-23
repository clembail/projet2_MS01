import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# --- CONFIGURATION ---
results_dir = "../results_PCG-AS"

data_points = []

print("--- Génération du graphe Question 2 (Scalabilité) ---")

# Trouver tous les fichiers CSV correspondant au pattern
csv_files = glob.glob(os.path.join(results_dir, "results_h_*.csv"))

for csv_file in csv_files:
    try:
        # 1. Extraction propre de h depuis le nom de fichier
        filename = os.path.basename(csv_file)
        # On enlève "res_h_" et ".csv" pour garder le nombre
        h_str = filename.replace("results_h_", "").replace(".csv", "")
        h = float(h_str)
        
        # 2. Lecture des données
        df = pd.read_csv(csv_file)
        
        if not df.empty:
            # Le nombre d'itérations total est la dernière valeur de la colonne 'iteration'
            last_iter = df['iteration'].iloc[-1]
            data_points.append((h, last_iter))
            
    except Exception as e:
        print(f"Erreur avec {csv_file}: {e}")

# 3. Tri des données par h décroissant (pour que la ligne soit tracée correctement)
data_points.sort(key=lambda x: x[0], reverse=True)

if not data_points:
    print("Aucune donnée trouvée !")
    exit()

# Séparation x et y
h_vals = [p[0] for p in data_points]
iter_vals = [p[1] for p in data_points]

# 4. Tracé
plt.figure(figsize=(10, 7))

# Semilogx car h varie sur des ordres de grandeur
plt.loglog(h_vals, iter_vals, 'o-', linewidth=2, markersize=8, color='blue', label="ASM (4 domaines)")

# Inversion de l'axe X : on aime voir h diminuer vers la droite (raffinage)
plt.gca().invert_xaxis()

plt.title(r"Scalabilité : Nombre d'itérations vs Finesse $h$ (Tol=$10^{-10}$)")
plt.xlabel(r"Finesse du maillage $h$ (échelle log)")
plt.ylabel("Nombre d'itérations à convergence")
plt.grid(True, which="both", linestyle='--', alpha=0.7)
plt.legend()

output_file = "tp3_exo3_question_2.png"
plt.savefig(output_file)
print(f"Graphique sauvegardé : {output_file}")
plt.show()