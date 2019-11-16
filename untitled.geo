SetFactory("OpenCASCADE");

//+
Point(1) = {-1, -1, -0, 1.0};
//+
Point(2) = {-1, 1, -0, 1.0};
//+
Point(3) = {1, 1, -0, 1.0};
//+
Point(4) = {1, -1, -0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Disk(2) = {0, 0, 0, 0.1, 0.1};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

//+
Extrude {0.0,0.0,1.0}{ 
    Surface{1}; 
    } 
//+


Mesh.MshFileVersion = 1;//+
Circle(2) = {0, 0, 0, 0.1, 0, 2*Pi};
//+
SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
//+
SetFactory("Built-in");
//+
Coherence;
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, 10, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 1, 0.3, 2*Pi};
