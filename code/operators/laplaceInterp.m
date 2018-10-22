function M = laplaceInterp(L, knownInd, method)

% Calculates the operator matrix for Laplace interpolation on a
% triangular mesh.
%
% Input:
%     L: Laplace operator matrix
%     knownInd: Indices of nodes, where the values are known
%     method:
%         'A': Laplacian = 0 at nodes, where the value is unknown
%         'B' (default): Minimizes the Laplacian at all nodes
% Output:
%     M: Interpolation matrix M
%
% This function implements:
% Oostendorp T, Oosterom A, Huiskamp G (1989).
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.
%
% Written by Steffen Schuler in 2018.
% Based on mesh_laplacian_interp by Darren Weber.

switch nargin
    case 2
        method = 'B';
    case 3
    otherwise
        error('Wrong number of input arguments.');
end

k = length(knownInd);
n = length(L);

% find unknown indices
unknownInd = setdiff(1:n, knownInd);

% reorder rows and cols of L
indL = [knownInd, unknownInd];
L = L(indL, :);
L = L(:, indL);

L11 = L(1:k    , 1:k    );
L12 = L(1:k    , (k+1):n);
L21 = L((k+1):n, 1:k    );
L22 = L((k+1):n, (k+1):n);

clear L;

switch method
    case 'A' 
        M = -L22 \ L21;
    case 'B'
        M = -[L12; L22] \ [L11; L21];
    otherwise
        warning('Unknown method "%s". Using default method "B".', method);
        M = -[L12; L22] \ [L11; L21];
end

clear L11 L12 L21 L22;

% append an identity matrix, taking care of the known values
M = [speye(k); M];

% reorder rows of M
[~,indM] = sort(indL);
M = M(indM, :);

end
