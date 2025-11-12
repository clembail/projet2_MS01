#!/bin/zsh

# =======================================================
# 1. PARAMÈTRES À CONFIGURER
# =======================================================

# La liste des finesses de maillage (h) à générer
H_VALUES=(0.05 0.025 0.01 0.005 0.0025 0.001)
#
CPP_EXECUTABLE="./build/tp2-ex1"


# =======================================================
# 2. VÉRIFICATIONS
# =======================================================

# Vérifier si l'exécutable C++ existe
if [ ! -f "$CPP_EXECUTABLE" ]; then
    echo "Erreur: Exécutable C++ non trouvé à '$CPP_EXECUTABLE'"
    echo "Veuillez compiler votre projet avant de lancer ce script."
    exit 1
fi

echo "--- Lancement du script de simulation (pour Q5 et Q6) ---"
echo "Exécutable C++ : $CPP_EXECUTABLE"


# =======================================================
# 3. BOUCLE DE SIMULATION
# =======================================================

# Boucle sur chaque valeur de h dans la liste
for h in "${H_VALUES[@]}"; do
    
    echo "" # Ligne vide pour la clarté
    echo "=================================================="
    echo "TRAITEMENT POUR h = $h"
    echo "=================================================="

    # Définir les noms de fichiers dynamiquement
    MESH_FILE="./solutions/mesh_h_${h}.mesh"
    CSV_FILE="./results_PCG/results_h_${h}.csv"

    # --- Lancement du solveur C++ ---
    echo "2. Lancement du solveur C++..."
    
    # Commande pour lancer votre C++
    $CPP_EXECUTABLE $MESH_FILE $CSV_FILE
    
    if [ $? -ne 0 ]; then
        echo "   ERREUR: Le solveur C++ a échoué pour h = $h."
    else
        echo "   Succès: Résultats '$CSV_FILE' créés."
    fi
done

echo ""
echo "=================================================="
echo "--- Script terminé. Tous les maillages ont été traités. ---"
