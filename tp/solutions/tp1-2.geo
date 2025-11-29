SetFactory("OpenCASCADE");

// Taille de maille
If (!Exists(h))
    h = 1e-2;
EndIf

Box(1) = {0, 0, 0, 1, 1, 1};

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Physical Volume(1) = {1};
