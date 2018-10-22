% Writes a VTK file (legacy: vtk or XML: vtp, vtu)
% Also supported are: ply, stl, obj, off (only geometry, no data arrays)
%
% Syntax: writeVTK(struct inStruct, char filename, [bool verbose], [char dataMode]);
%
% inStruct must have the following fields:
%     points:    [numPoints x 3 double]
%     cells:     [numCells x maxNumPointsPerCell int32]
%     cellTypes: [numCells x 1 uint8]
%     pointData: [1 x 1 struct of point data arrays] (optional)
%     cellData:  [1 x 1 struct of cell data arrays] (optional)
%     fieldData: [1 x 1 struct of field data arrays] (optional)
%
% verbose is false on default
% 
% dataMode can be chosen from {'ascii', 'binary'} for vtk, ply, stl
%                         and {'ascii', 'binary', 'appended'} for vtp, vtu
%                         and {'ascii'} for obj, off
%
% Written in November 2017 by
% Nils Daub, Steffen Schuler
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu

function writeVTK(inStruct, filename, verbose) %#ok
