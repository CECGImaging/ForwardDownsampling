function res=MeanNodeDistance(mesh)
% function res=MeanNodeDistance(mesh)
% Calculates mean distance between two neighboring nodes
% Equals to mean element side length.
r12=mesh.p(mesh.e(:,2),:)-mesh.p(mesh.e(:,1),:);
r13=mesh.p(mesh.e(:,3),:)-mesh.p(mesh.e(:,1),:);
r23=mesh.p(mesh.e(:,3),:)-mesh.p(mesh.e(:,2),:);
r12=sqrt(sum(r12.*r12,2));
r13=sqrt(sum(r13.*r13,2));
r23=sqrt(sum(r23.*r23,2));
res=mean([r12;r13;r23]);