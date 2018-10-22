% Calculates an operator matrix for linear interpolation of values
% from one triangle mesh to another using barycentric coordinates
% 
% Syntax:
% barycentricInterp(struct sourceMesh, struct targetMesh, ...
%                   double normalLength, [bool verbose])
% 
% The function creates a matrix 'interpMat' in the workspace.
% This matrix can be used to linearly interpolate values from all points in
% sourceMesh to all points in targetMesh.
% 
% Additionally, two debug files 'barycentricInterp_intersections.vtp' and
% 'barycentricInterp_target.vtp' are created in the current folder.
% 
% sourceMesh and targetMesh must have the following fields
% (as created by readVTK):
%     points:    [numPoints x 3 double]
%     cells:     [numCells x 3 int32]
% 
% normalLength defines the maximum search distance in positive and negative 
% normal direction.
% 
% verbose is false on default
% 
% Written in October 2018 by
% Steffen Schuler
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu

function barycentricInterp(sourceMesh, targetMesh, normalLength, verbose) %#ok
