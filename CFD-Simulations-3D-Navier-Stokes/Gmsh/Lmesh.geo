// Gmsh project created on Sat May  3 20:58:02 2014

Point(1) = {0, -2, 0, 1.0};
Point(2) = {8, -2, 0, 1.0};
Point(3) = {8, 6, 0, 1.0};
Point(4) = {4, 6, 0, 1.0};
Point(5) = {4, 2, 0, 1.0};
Point(6) = {0, 2, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(7) = {1, 2, 3, 4, 5, 6};
Plane Surface(8) = {7};
DefineConstant[ lc = { 0.01, Path "Gmsh/Parameters"}];

