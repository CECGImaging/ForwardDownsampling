function [points, elements, nop,noe]=CreateSphere(origin, radius, sub)
%function [points, elements, nop, noe]=CreateSphere(origin, radius, sub)
%sub = 1, 2, or 3 defines the number of points in the mesh
%Builds a spherical mesh with given origin and radius
%Uses meshes generated with Matlab routine trisphere.m written by Finn
%Lindgren, Lund University, Sweden
%(http://www.maths.lth.se/matstat/staff/finn/prog/index.html.en)
if sub==1
    load sphere_sub1
elseif sub==2
    load sphere_sub2
elseif sub==3
    load sphere_sub3
elseif sub==4
    load sphere_sub4
else
    error('Improper sub value')
end
nop=length(p);
noe=length(elements);
p=p.*radius;
trans=ones(nop,1);
trans=trans*origin;
points=p'+trans;
