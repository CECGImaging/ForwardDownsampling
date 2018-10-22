function [points,triangles]=LoadTrimesh(filename)
% function [points,triangles]=LoadTrimesh(filename)
% load triangle mesh from ascii file
% file format: 
% non 
% nodes (non x 3 -matrix of node coordinates)
% noe
% elements (noe x 3 -matrix of indexes to node coordinates)
%
% where non = number of nodes, noe = number of elements

fp=fopen(filename,'r');
nop=fscanf(fp,'%d',1);
points=zeros(nop,3);
for I=1:nop
    points(I,:)=fscanf(fp,'%f',3)';
end
not=fscanf(fp,'%d',1);
triangles=zeros(not,3);
for I=1:not
    triangles(I,:)=fscanf(fp,'%d',3)';
end
fclose(fp);
