// Gmsh project created on Sun Nov  1 11:36:12 2015
cl__1 = 20;
cl__2 = 10;
Point(1) = {0, -200, 0, 20};
Point(2) = {300, -200, 0, 10};
Point(3) = {508.3778132, -1381.76930361, 0, 20};
Point(4) = {902.300914405, -1312.31003255, 0, 20};
Point(5) = {723.215385232, -296.665527018, 0, 10};
Point(6) = {767.819926237, 0, 0, 10};
Point(7) = {831.659443104, 293.128838748, 0, 10};
Point(8) = {1795.49895998, 1256.96835562, 0, 20};
Point(9) = {1724.78828186, 1327.67903374, 0, 20};
Point(10) = {597.109248121, 200, 0, 10};
Point(11) = {0, 200, 0, 20};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {9, 10};
Line(8) = {10, 11};
Line(9) = {11, 1};
Circle(10) = {5, 6, 7};
Line Loop(16) = {1, 2, 3, 4, 10, 5, 6, 7, 8, 9};
Plane Surface(16) = {16};