import pandas as pd
import matplotlib.pyplot as plt
import glob

print("--- Analyse pour la Question 6 ---")

# 1. Trouver tous les fichiers de résultats
csv_files = glob.glob("../results_CG/results_h_*.csv")
results = []

for f in csv_files:
    # Extraire h du nom de fichier 'results_h_0.025.csv'
    try:
        h = float(f.split('_')[-1].replace('.csv', ''))
    except:
        continue
        
    # Lire le CSV
    df = pd.read_csv(f)
    
    # Le nombre d'itérations est juste le nombre de lignes
    # (ou la dernière valeur de 'iteration')
    iterations = df['iteration'].iloc[-1]
    
    results.append((h, iterations))

# 2. Trier les résultats par h pour un joli tracé
results.sort(key=lambda x: x[0])
h_values = [r[0] for r in results]
iter_values = [r[1] for r in results]

# 3. Tracer (répond à Q6)
plt.figure(figsize=(10, 7))
plt.plot(h_values, iter_values, 'o-')
plt.title("Nombre d'itérations CG vs. Finesse du maillage (h)")
plt.xlabel("Finesse du maillage (h)")
plt.ylabel("Nombre d'itérations (pour atteindre 10e-6)")
plt.xscale('log') # Souvent utile de voir h en log
plt.yscale('log') # Souvent utile de voir les itérations en log
plt.grid(True, which="both")
plt.savefig("tp1_exo2_question_6.png")
print("Graphique 'tp1_exo2_question_6.png' sauvegardé.")
