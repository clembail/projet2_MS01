SetFactory("OpenCASCADE");

// Taille de maille
h = 1e-2;

// Définit le volume géométrique 3D (le cube)
// Il reçoit le tag 1
Box(1) = {0, 0, 0, 1, 1, 1};

// Applique la taille de maille
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

// LA LIGNE CRUCIALE :
// Dit à Gmsh: "L'entité géométrique 'Box(1)'
// est un 'Physical Volume'. Il faut le mailler."
Physical Volume(1) = {1};