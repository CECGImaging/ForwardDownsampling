function [actTimes,precomp] = globalActTimes(sig, sampPeriod, mesh, neighborhoodSize, refNode, gaussSigma, mode, precomp)

% Estimates activation times by
% - Calculating time delays for pairs of nearby nodes using
%   cross-correlation of a derivative signal
% - Solving the equation differenceMatrix * activationTimes = delays
%   in a least-squares sense
%
% This implements the "global activation time" approach described in:
% R. Dubois et al., "Global and Directional Activation Maps for 
% Cardiac Mapping in Electrophysiology, CinC 2012
%
% Input: 
%     sig:              Signal (numNodes x numTimesteps)
%     sampPeriod:       Sampling period
%     mesh:             Mesh in BEM-library format
%     neighborhoodSize: Number of edges between pairs of nodes
%     refNode:          Reference node used to remove the const time offset
%     gaussSigma:       Sigma used for temporal gaussian filtering of sig
%     mode:             Mode of derivative. One of:
%                      'spatial' (default), 'temporal', 'spatiotemporal'
%
% Output:
%     actTimes: Activation times in the same unit as sampPeriod
%     precomp:  Struct with precomputed intermediate results, which can be 
%               used to speed up further executions on the same geometry
%
% Usage without precomputed struct:
% [actTimes,precomp] = globalActTimes(sig, sampPeriod, mesh, neighborhoodSize, refNode, gaussSigma, [mode])
%
% Usage with precomputed struct:
% actTimes = globalActTimes(sig, sampPeriod, [], [], [], gaussSigma, mode, precomp)
%
% Written by Steffen Schuler in 2018.

%%

% addpath('../gptoolbox/mesh');
% addpath('../gptoolbox/matrix');
% addpath('../operators');

%% Check input arguments

noPrecomp = true;
if nargin < 7
    mode = 'spatial';
elseif nargin > 7
    noPrecomp = false;
end

%% Determine node neighborhood pairs

if noPrecomp
    [~,dist] = dijkstra_shepard(mesh.p, double(mesh.e), mesh.p);
    d = mean(mean(edge_lengths(mesh.p, mesh.e)));

    precomp.pairs = NaN(20*mesh.nop,2);
    numPairs = 0;

    for i = 1:mesh.nop
        neighs = find(dist(:,i) > (neighborhoodSize-0.5)*d & dist(:,i) < (neighborhoodSize+0.5)*d);
        numNeighs = numel(neighs);
        precomp.pairs(numPairs+1:numPairs+numNeighs,:) = [repmat(i,numNeighs,1) neighs];
        numPairs = numPairs + numNeighs;
    end

    precomp.pairs(isnan(precomp.pairs(:,1)),:) = [];
    precomp.pairs = unique(sort(precomp.pairs,2), 'rows');
end

%% Define derivative signal (derivSig) used to determine delays

% increase sampling rate using linear interpolation
upsampFactor = 10;
resampPeriod = sampPeriod/upsampFactor;
s = NaN(size(sig,1), upsampFactor*(size(sig,2)-1)+1);
for i = 1:size(sig,1)
    s(i,:) = interp1(0:size(sig,2)-1, sig(i,:), 0:1/upsampFactor:size(sig,2)-1);
end

if gaussSigma > 0
    gaussSigma = gaussSigma/resampPeriod;
    bt = sqrt(log(2))/(2*pi*gaussSigma);
    span = round(6*gaussSigma/2)*2;
    b = gaussdesign(bt, span, 1);
end

if gaussSigma > 0
    s = filtfilt(b,1,s')';
end
% for i = 1:size(s,1)
%     plot(s(i,:));
%     title(i);
%     waitforbuttonpress;
% end

if noPrecomp
    [precomp.Gx,precomp.Gy,precomp.Gz] = Gradient(mesh);
end
gradSig = sqrt((precomp.Gx*s).^2+(precomp.Gy*s).^2+(precomp.Gz*s).^2);
% plot(gradSig(1,:));

diffSig = (s(:,3:end)-s(:,1:end-2))/2;
% plot(diffSig(1,:));

switch mode
    case 'spatial'
        derivSig = gradSig;
    case 'temporal'
        derivSig = diffSig;
    case 'spatiotemporal'
        derivSig = gradSig(:,2:end-1) .* diffSig;
    otherwise
        warning('Unknown derivativeMode ''%s''. Using default mode ''spatial'' instead.', mode);
        derivSig = gradSig;
end
% for i = 1:size(derivSig,1)
%     plot(derivSig(i,:));
%     title(i);
%     waitforbuttonpress;
% end

%% Compute delays

delays = NaN(size(precomp.pairs,1),1);
parfor i = 1:size(precomp.pairs,1)
    [xc,lag] = xcorr(derivSig(precomp.pairs(i,1),:), derivSig(precomp.pairs(i,2),:));
    [~,ind] = max(xc);
    delays(i) = lag(ind);
    % plot(xc)
    % title(lag(ind))
    % waitforbuttonpress
end

%% Fix the activation time of the reference node
%  to resolve the ambiguity in constant time offset, which cannot be
%  determined from delays alone

if noPrecomp
    precomp.refT = [refNode 0];
end

%% Build difference matrix and compute least-squares solution

if noPrecomp
    D = zeros(size(precomp.pairs,1), mesh.nop);
    for i = 1:size(precomp.pairs,1)
        D(i,precomp.pairs(i,1)) = 1;
        D(i,precomp.pairs(i,2)) = -1;
    end

    I = zeros(size(precomp.refT,1), mesh.nop);
    for i = 1:size(precomp.refT,1)
        I(i,precomp.refT(i,1)) = 1;
    end

    M = sparse(cat(1,I,D));
    precomp.invM = (M'*M) \ M';
end

v = cat(1,precomp.refT(:,2),delays);
T = precomp.invM * v;

%% Determine absolute activation time offset
%  by aligning all derivSigs, building a template
%  and detecting the maximum of this template

derivSigAligned = delayseq(derivSig', -T);
[~,offsetT] = max(mean(derivSigAligned,2));

actTimes = resampPeriod*(T+offsetT);

end
