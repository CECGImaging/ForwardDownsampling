function L = laplaceBeltrami(mesh)
% Calculates the operator matrix for the surface Laplacian
% of a triangular mesh (Laplace-Beltrami operator).
% Written by Steffen Schuler in March 2017,
% based on libigl.
% Input:  Mesh in the BEM-library format.
% Output: Laplace operator matrix.

[Gx,Gy,Gz] = gradient(mesh, 0);
G = [Gx; Gy; Gz];
T = spdiags(repmat(mesh.a,3,1), 0, 3*mesh.noe, 3*mesh.noe);
L = -G' * T * G;

end
