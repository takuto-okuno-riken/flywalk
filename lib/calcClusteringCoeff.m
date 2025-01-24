% Calculate clustering coefficient of adjacency matrix.
% returns clustering coefficient in each node (C), averaged clustering coefficient (gR),
%         all cluster count (count), clustering matrix (S3).
% input:
%  S            adjacency matrix (node x node, logical)

function [C, aC, count, S3] = calcClusteringCoeff(S)
    nlen = size(S,1);
    S2 = S | S'; % in or out degree
    K = sum(S2,2); % node total degree
    S3 = zeros(nlen,nlen,'logical');

    % set pool num. this calculation takes time. we need big pool num.
%    delete(gcp('nocreate')); % shutdown pools
%    parpool(12);

%    parfor i=1:nlen
    for i=1:nlen
%        disp(num2str(i));
        logis = S2(i,:);
%        Si = S3(logis,:); % this one is slow
%        Si(:,~logis) = 0;
%        Si = sum(Si,1);
%        S4b(i,:) = Si;

        Si = S2(logis,logis);
        Si = sum(Si,1);
        S3(i,logis) = Si;
    end
    Tr = sum(S3,2); % triangle degree
    C = Tr ./ (K.*(K-1));
    aC = nanmean(C);
    count = length(find(Tr>0));
end
