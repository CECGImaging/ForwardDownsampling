function res=RotatePointset(points,Rx,Ry,Rz)
% function res=RotatePointset(points,Rx,Ry,Rz)
% Ri = rotation around i-axis (radians)
% points = N x 3
%
% Rotation directions:
% x: rotates y-axis towards z
% y: rotates z-axis towards x
% z: rotates x-axis towards y

Rxmat=[1 0 0; 0 cos(Rx) -sin(Rx); 0 sin(Rx) cos(Rx)];
Rymat=[cos(Ry) 0 sin(Ry);0 1 0; -sin(Ry)  0 cos(Ry)];
Rzmat=[cos(Rz) -sin(Rz) 0;sin(Rz) cos(Rz) 0; 0 0 1];
res=Rxmat*points';
res=Rymat*res;
res=(Rzmat*res)';