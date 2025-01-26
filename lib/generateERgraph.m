% generate (directed) Erdos-Renyi random graph by probability
% input:
%  n            node size
%  prob         probability (density)
%  isSparse     generate sparse matrix (default: false)

function E = generateERgraph(n, prob, isSparse)
    if nargin < 3, isSparse = false; end

    E = logical(sprand(n,n,prob)); % this includes diag elements.
    if isSparse
        E = E & ~speye(n,n); % remove diag.
    else
        E = full(E); % this includes diag elements.
        E = E & ~eye(n,n,'logical'); % remove diag.
    end
end
