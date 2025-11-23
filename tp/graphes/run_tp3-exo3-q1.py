import pandas as pd
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
results_dir = "../results_PCG-AS"
# On choisit quelques valeurs représentatives pour ne pas surcharger le graphe
h_to_plot = [0.025, 0.01, 0.005, 0.0025] 

plt.figure(figsize=(10, 7))
markers = ['o-', 's-', '^-', 'd-', 'x-']

print("--- Génération du graphe Question 1 (Convergence) ---")

plot_count = 0
for i, h in enumerate(h_to_plot):
    csv_file = os.path.join(results_dir, f"results_h_{h}.csv")
    
    if os.path.exists(csv_file):
        try:
            df = pd.read_csv(csv_file)
            
            # Tracé : Itération en X, Résidu Relatif en Y (échelle log)
            label = f"h = {h}"
            marker = markers[i % len(markers)]
            
            # markevery=10 évite d'avoir trop de symboles sur la courbe
            plt.semilogy(df['iteration'], df['rel_residual'], marker, label=label, markevery=5, markersize=6)
            plot_count += 1
        except Exception as e:
            print(f"Erreur lecture {csv_file}: {e}")
    else:
        print(f"Fichier manquant pour h={h} : {csv_file}")

if plot_count > 0:
    plt.title(r"Convergence PCG - Schwarz Additif ($N_{dom}=4, n_l=2$)")
    plt.xlabel("Nombre d'itérations $k$")
    plt.ylabel(r"Résidu relatif $\|b - Ax_k\|_2 / \|b\|_2$")
    plt.grid(True, which="both", linestyle='--', alpha=0.7)
    plt.legend()
    
    output_file = "tp3_exo3_question_1.png"
    plt.savefig(output_file)
    print(f"Graphique sauvegardé : {output_file}")
    plt.show()
else:
    print("Aucune donnée à tracer.")