% Reads a VTK file (legacy: vtk or XML: vtp, vtu)
% Also supported are: ply, stl, obj, (off not yet)
%
% Syntax: outStruct = readVTK(char filename, [bool verbose])
%
% Written in June 2017 by
% Steffen Schuler
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu

function outStruct = readVTK(filename, verbose) %#ok
