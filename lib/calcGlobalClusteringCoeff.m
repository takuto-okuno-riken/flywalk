% Calculate global clustering coefficient of adjacency matrix.
% returns global clustering coefficient (coeff), triangle num (triangle) and all triplet (allTriplet).
% input:
%  A            adjacency matrix (node x node, logical)
%  isWeight     use weight mat for out degrees (default: true)

function [coeff, triangle, allTriplet] = calcGlobalClusteringCoeff(A, isWeight)
    if nargin < 2, isWeight = true; end

    B = A + A'; % weight type (to be symmetric). this results large triplet
    tri1 = sum((B*B).*B','all') / 2;
%    triplet = tri1;
    C = A | A'; % binary type (to be symmetric). this results small triplet
    tri2 = sum((C*C).*C','all') / 2;
%    triplet = tri2;
    triplet = (tri1 + tri2) / 2; % this close to graph-tool compatible. hmm ...

    if isWeight
        allTriplet = outDegrees(B); % this is graph-tool compatible.
    else
        allTriplet = outDegrees(C);
    end
    coeff = triplet / allTriplet;
    triangle = floor(triplet/3);
end

function c = outDegrees(A)
    K = sum(A,2);
    c = sum(K.*(K-1))/2;
%    c = 0;
%    for i=1:size(A,1)
%        k = sum(A(i,:));
%        c = c + k*(k-1)/2;
%    end
end

