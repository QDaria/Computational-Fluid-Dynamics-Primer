cl__1 = 50;
cl__2 = 20;
Point(1) = {0, -200, 0, 50};
Point(2) = {300, -200, 0, 20};
Point(3) = {404.1889066, -790.884651807, 0, 50};
Point(4) = {798.112007805, -721.42538074, 0, 50};
Point(5) = {723.215385232, -296.665527018, 0, 20};
Point(6) = {767.819926237, 0, 0, 20};
Point(7) = {831.659443104, 293.128838748, 0, 20};
Point(8) = {1231.659443106, 693.12883875, 0, 50};
Point(9) = {1160.94876499, 763.839516869, 0, 50};
Point(10) = {597.109248121, 200, 0, 20};
Point(11) = {0, 200, 0, 50};
Point(12) = {0, -200, 185, 50};
Point(13) = {300, -200, 185, 20};
Point(17) = {404.1889066, -790.884651807, 185, 50};
Point(21) = {798.112007805, -721.42538074, 185, 50};
Point(25) = {723.215385232, -296.665527018, 185, 20};
Point(29) = {767.819926237, 0, 185, 20};
Point(30) = {831.659443104, 293.128838748, 185, 20};
Point(34) = {1231.659443106, 693.12883875, 185, 50};
Point(38) = {1160.94876499, 763.839516869, 185, 50};
Point(42) = {597.109248121, 200, 185, 20};
Point(46) = {0, 200, 185, 50};
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
Line(18) = {12, 13};
Line(19) = {13, 17};
Line(20) = {17, 21};
Line(21) = {21, 25};
Circle(22) = {25, 29, 30};
Line(23) = {30, 34};
Line(24) = {34, 38};
Line(25) = {38, 42};
Line(26) = {42, 46};
Line(27) = {46, 12};
Line(29) = {1, 12};
Line(30) = {2, 13};
Line(34) = {3, 17};
Line(38) = {4, 21};
Line(42) = {5, 25};
Line(46) = {7, 30};
Line(50) = {8, 34};
Line(54) = {9, 38};
Line(58) = {10, 42};
Line(62) = {11, 46};
Line Loop(16) = {1, 2, 3, 4, 10, 5, 6, 7, 8, 9};
Plane Surface(16) = {16};
Line Loop(31) = {1, 30, -18, -29};
Ruled Surface(31) = {31};
Line Loop(35) = {2, 34, -19, -30};
Ruled Surface(35) = {35};
Line Loop(39) = {3, 38, -20, -34};
Ruled Surface(39) = {39};
Line Loop(43) = {4, 42, -21, -38};
Ruled Surface(43) = {43};
Line Loop(47) = {10, 46, -22, -42};
Ruled Surface(47) = {47};
Line Loop(51) = {5, 50, -23, -46};
Ruled Surface(51) = {51};
Line Loop(55) = {6, 54, -24, -50};
Ruled Surface(55) = {55};
Line Loop(59) = {7, 58, -25, -54};
Ruled Surface(59) = {59};
Line Loop(63) = {8, 62, -26, -58};
Ruled Surface(63) = {63};
Line Loop(67) = {9, 29, -27, -62};
Ruled Surface(67) = {67};
Line Loop(68) = {18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
Plane Surface(68) = {68};
Surface Loop(1) = {16, 68, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67};
Volume(1) = {1};
