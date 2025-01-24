% generate directed configuration (CFG) model graph
% this is equivalent to the Maslov–Sneppen edge-swapping null model
%  References: Maslov, S. & Sneppen, K., Science 296, 910–913 (2002).
% input:
%  S            adjacency matrix (node x node, logical)

function S = generateCFGgraph(S)
    n = size(S,1);   % node num
    [row,col] = find(S);
    L = length(row); % total link num
    
    retryNum = ceil(L/(n-1)); % retry number
    
    for itr=1:L
        for r=1:retryNum
            i=ceil(L*rand());
            j=ceil(L*rand());
            a=row(i); b=col(i);
            c=row(j); d=col(j);
    
            % check random index
            if a==c || a==d || b==c || b==d
                continue;
            end
            % check elements are zero
            if S(a,d) || S(c,b)
                continue;
            end
    
            % swap matrix elements (0 to 1)
            S(a,d)=1; S(a,b)=0;
            S(c,b)=1; S(c,d)=0;
            
            col(i) = d;
            col(j) = b;
            break;
        end
    end
end
