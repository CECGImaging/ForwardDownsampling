function SaveTrimesh(filename,points,elements)
% function SaveTrimesh(filename,points,elements)
% saves a triangle mesh to an ascii file
% file format: 
% non 
% nodes (non x 3 -matrix of node coordinates)
% noe
% elements (noe x 3 -matrix of indexes to node coordinates)
%
% where non = number of nodes, noe = number of elements

fp=fopen(filename,'w');
fprintf(fp,'%d\n',length(points));
[r c]=size(points);
if r==3
    points = points';
end
for I=1:length(points)
    fprintf(fp,'%f\t%f\t%f\n',points(I,:));
end

[r c]=size(elements);
if r==3
    elements = elements';
end

fprintf(fp,'%d\n',length(elements));
for I=1:length(elements)
    fprintf(fp,'%d\t%d\t%d\n',elements(I,:));
end
fclose(fp);
