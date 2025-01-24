% Calculate reciprocity of adjacency matrix.
% returns reciprocal count in each node (R), global reciprocity (gR), 
%         all reciprocal count (count), reciprocal matrix (S2).
% input:
%  S            adjacency matrix (node x node, logical)

function [R, gR, count, S2] = calcReciprocity(S)
    L = sum(S,'all'); % total number of links
    S2 = S & S';
    R = sum(S2,2); count = length(find(R>0));
    Lr = sum(S2,'all'); % get reciprocal number
    gR = Lr/L;
end
