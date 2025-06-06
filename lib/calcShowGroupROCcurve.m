%%
% Calculate and show ROC curve with group data matrix
% returns X, Y (node x elements) and AUC (node x 1)
% input:
%  gmb         Ground Truth binary matrix (node x elements)
%  tm          evaluating group data matrix (node x elements)
%  groupStr    group name for ROC curve plot
%  isShow      flag to plot ROC curve
%  col         color of ROC curve

function [X, Y, AUC] = calcShowGroupROCcurve(gmb, tm, groupStr, isShow, col)
    if nargin < 5, col = []; end
    injNum = size(tm, 1);
    n = size(tm, 2);
    X = zeros(injNum,n);
    Y = zeros(injNum,n);
    AUC = nan(injNum,1);
    [~,I] = sort(tm,2,'descend');
    if isShow, figure; end
    for i=1:injNum
        r = uint32(gmb(i,I(i,:))); % 0 to 4294967295
        tc = sum(gmb(i,:));
        if tc==0
            disp(['no active ground truth i=' num2str(i)]);
            continue;
        end
        fc = n - tc;
        Y(i,:) = double(cumsum(r,2)) / double(tc);
        X(i,:) = double(cumsum(1-r,2)) / double(fc);
        AUC(i) = trapz(X(i,:),Y(i,:));
        if isShow
            hold on;
            if isempty(col), plot(X(i,:),Y(i,:));
            else, plot(X(i,:),Y(i,:),'Color',col); end
            hold off;
        end
    end
    if isShow
        xlim([0 1]); ylim([0 1]); daspect([1 1 1]);
        hold on; plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]); hold off;
        title(['ROC curve of ' groupStr]);
    end
end
