%%
% load SWC neuron tree file
% This function can call Trees toolbox (https://www.treestoolbox.org/index.html), if exist.
% input:
%  fname        SWC file name
%  isTrees      call Trees toolbox load_tree function, if exist (default: false)

function swc = loadSwc(fname, isTrees)
    if nargin < 2, isTrees = false; end
    if isTrees && exist('load_tree','file')
        swc = load_tree(fname); % Trees toolbox
    else
        X = readmatrix(fname,'FileType','text','OutputType','double');
        swc = nan(1,size(X,2)-1);
        idx = X(:,1);
        swc(idx,:) = nan;
        swc(idx,:) = X(:,2:end);
    end
end