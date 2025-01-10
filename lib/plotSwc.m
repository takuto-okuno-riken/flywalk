%%
% Plotting SWC neuron tree.
% This function can call Trees toolbox (https://www.treestoolbox.org/index.html), if exist.
% input:
%  swc          swc matrix (list x items)
%  color        line color (1 x 3 [0 1])
%  lineWidth    line width (default: 0.5)
%  isTrees      call Trees toolbox plot_tree function, if exist (default: false)

function plotSwc(swc, col, lineWidth, isTrees)
    if nargin < 4, isTrees = false; end
    if nargin < 3, lineWidth = 0.5; end
    if isTrees && exist('plot_tree','file')
        swc = plot_tree(swc, col); % Trees toolbox
        return;
    end
    V = swc(:,2:4);
    F = [];
    for i=1:size(swc,1)-1
        if swc(i,end) > 0
            F(end+1,:) = [i swc(i,end)];
        end
    end
    patch('Faces',F,'Vertices',V,'FaceColor','none','EdgeColor',col,'LineWidth',lineWidth);
end