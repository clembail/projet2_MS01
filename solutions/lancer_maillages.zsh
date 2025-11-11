#!/bin/zsh

# 1. Le fichier .geo source qui est paramétré
GEO_FILE="generer_carre.geo"

# 2. Liste des finesses de maillage (h) à générer (de la Question 5)
H_VALUES=(0.05 0.025 0.01 0.005 0.0025 0.001)

echo "--- Début de la génération des maillages ---"

# 3. Boucle sur chaque valeur de h
for h in $H_VALUES; do
  echo "Génération du maillage pour h = $h ..."
  
  # Définir le nom du fichier de sortie
  OUTPUT_FILE="mesh_h_${h}.mesh"
  
  # 4. La commande Gmsh
  #    -2 : Génère un maillage 2D
  #    -setnumber h $h : Définit la variable 'h' dans le .geo
  #    -format mesh : Spécifie le format de sortie (MEDIT/GMF, celui attendu par votre code)
  #    -o $OUTPUT_FILE : Spécifie le fichier de sortie
  gmsh -2 -setnumber h $h $GEO_FILE -o $OUTPUT_FILE
  
  if [ $? -eq 0 ]; then
    echo "Succès: Maillage sauvegardé dans $OUTPUT_FILE"
  else
    echo "Erreur: Échec de la génération pour h = $h"
  fi
done

echo "--- Génération de tous les maillages terminée. ---"
