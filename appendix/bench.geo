cl0 = 0.1;
cl1 = 0.001;
cl2=0.005;
cl3=0.01;

r=0.05;
r1=0.06;
r2=0.08;

L=2.2;
W=0.41;

cx=0.2;
cy=0.2;

Point(1) = {0, 0, 0, cl0};
Point(2) = {L, 0, 0, cl0};
Point(3) = {L, W, 0, cl0};
Point(4) = {0, W, 0, cl0};
Point(5) = {cx, cy, 0, cl1};
Point(6) = {cx-r, cy, 0, cl1};
Point(7) = {cx+r, cy, 0, cl1};
Point(8) = {cx, cy-r, 0, cl1};
Point(9) = {cx, cy+r, 0, cl1};
Point(10) = {cx-r1, cy, 0, cl2};
Point(11) = {cx+r1, cy, 0, cl2};
Point(12) = {cx, cy-r1, 0, cl2};
Point(13) = {cx, cy+r1, 0, cl2};
Point(20) = {cx-r2, cy, 0, cl3};
Point(21) = {cx+r2, cy, 0, cl3};
Point(22) = {cx, cy-r2, 0, cl3};
Point(23) = {cx, cy+r2, 0, cl3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 9};
Circle(6) = {9, 5, 7};
Circle(7) = {7, 5, 8};
Circle(8) = {8, 5, 6};
Circle(10) = {10, 5, 13};
Circle(11) = {13, 5, 11};
Circle(12) = {11, 5, 12};
Circle(13) = {12, 5, 10};
Circle(20) = {20, 5, 23};
Circle(21) = {23, 5, 21};
Circle(22) = {21, 5, 22};
Circle(23) = {22, 5, 20};
Line Loop(11) = {1, 2, 3, 4};
Line Loop(12) = {8, 7, 6, 5};
Plane Surface(111) = {11 , 12};

Point {6} In Surface {111};
Point {7} In Surface {111};

Line { 10, 11, 12, 13} In Surface {111};
Line { 20, 21, 22, 23} In Surface {111};

Physical Line(12) = {4};
Physical Line(13) = {2};
Physical Line(14) = {1, 3};
Physical Line(15) = {5, 6, 7, 8};
Physical Surface(16) = {111};

Mesh.Smoothing = 100;
Smoother Surface{111} = 100;
Mesh.Algorithm = 8; // delquad
