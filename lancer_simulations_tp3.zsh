#!/bin/zsh

# --- CONFIGURATION ---
EXECUTABLE="./build/tp3-exo3"
MESH_DIR="./tp/solutions"
OUTPUT_DIR="./tp/results_PCG-AS"

# Liste des finesses h à tester
H_VALUES=(0.05 0.025 0.01 0.005 0.0025 0.001)

# Création du dossier de résultats s'il n'existe pas
mkdir -p $OUTPUT_DIR

echo "=============================================="
echo "   Lancement des simulations Schwarz Additif"
echo "=============================================="

# Vérification de la présence de l'exécutable
if [[ ! -x "$EXECUTABLE" ]]; then
    echo "Erreur : L'exécutable $EXECUTABLE est introuvable ou non exécutable."
    echo "Essayez : chmod +x $EXECUTABLE"
    exit 1
fi

for h in $H_VALUES; do
    MESH_FILE="${MESH_DIR}/mesh_h_${h}.mesh"
    CSV_OUTPUT="${OUTPUT_DIR}/results_h_${h}.csv"

    echo -n "-> Traitement h = $h ... "

    if [[ -f "$MESH_FILE" ]]; then
        # Exécution du code C++
        # Redirection de la sortie standard vers /dev/null pour ne pas polluer le terminal
        $EXECUTABLE "$MESH_FILE" "$CSV_OUTPUT" > /dev/null
        
        # Vérification rapide si le fichier CSV a été créé
        if [[ -f "$CSV_OUTPUT" ]]; then
            echo "OK (Résultats dans $CSV_OUTPUT)"
        else
            echo "ERREUR (Pas de fichier de sortie généré)"
        fi
    else
        echo "IGNORÉ (Maillage $MESH_FILE introuvable)"
    fi
done

echo "=============================================="
echo "   Simulations terminées."
echo "=============================================="