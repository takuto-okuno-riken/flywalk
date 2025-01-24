% generate (directed) Erdos-Renyi random graph by probability
% input:
%  n            node size
%  prob         probability (density)

function E = generateERgraph(n, prob)
    E = full(logical(sprand(n,n,prob))); % this includes diag elements.
    E = E & ~eye(n,n,'logical'); % remove diag.
end
