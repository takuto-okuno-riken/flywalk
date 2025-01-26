function testGraph
%%{
    % Erdos-Renyi random directed graph
    p = 0.1;
    n = 1000;
    E = generateERgraph(n,p);

    figure; imagesc(E); daspect([1 1 1]); xlabel('target'); ylabel('source');
    esum = sum(E,'all');
    disp(['density=' num2str(esum/(n*(n-1)))]);

    % histogram of in-degree and out-degree follow Poisson distribution.
    Eid = sum(E,1);
    Eod = sum(E,2);
    bin = [5*n/100:1:15*n/100] ;
    figure; histogram(Eid,bin); hold on; histogram(Eod,bin); hold off;
    mGid = mean(Eid);
    mGod = mean(Eod);
    [hGid,bin2] = histcounts(Eid,bin);
%    [lambdahat,lambdaci] = poissfit(Gid);

    y = poisspdf(bin,mGid) * n; % use mean value as possison lambda
    hold on; plot(bin,y); hold off;

    % fitting version
    n = sum(hGid);
    pd = fitdist(bin(1:end-1)','Poisson','Frequency',hGid'); % this shows same fitting.
    y2 = n * pdf(pd,bin);
    hold on; plot(bin,y2); hold off;

    % show QQ plot
    figure; qqplot(hGid,y);
    figure; qqplot(hGid,y2);

    % show CDF plot
    figure;
    cdfplot(hGid);
    hold on; cdfplot(y); hold off;
    legend('Empirical CDF','Poisson CDF','Location','best')

    % chi2 test
    [h,p,st] = chi2gof(bin(1:end-1),'Ctrs',bin(1:end-1),'Frequency',hGid, ...
        'Expected',y2(1:end-1),'NParams',1);
    disp(['chi2 test. vs. possison distribution p=' num2str(p)]);

%    [h,p,stats] = chi2gof(Gid); % this test for normal distribution

    % directed configuration (CFG) model
    G = generateCFGgraph(E);
    Gid = sum(G,1);
    God = sum(G,2);

    % in-degree and out-degree should be same
    d1 = sum(abs(Gid-Eid));
    d2 = sum(abs(God-Eod));
    disp(['in-degree diff=' num2str(d1) ', out-degree diff=' num2str(d2)]);
%}

    %-----------------------------------------------------------------------
    % generate ER random directed graph for calc cluster coeff.
    ns = [6 8 10 20 50 100 500 1000 10000];
    p = 0.02;
%{
    % save mat file to be checked by graph-tool.
    for i = 1:length(ns)
        n = ns(i);
        for j = 1:10
            E = generateERgraph(n,p,true);
            save(['d:\work\gs\erdir\erdir' num2str(n) '-' num2str(j) '_' num2str(p*100) '.mat'],'E');
        end
    end
%}
    for i = 1:length(ns)
        n = ns(i);
        for j = 1:10
            load(['d:\work\gs\erdir\erdir' num2str(n) '-' num2str(j) '_' num2str(p*100) '.mat']);
            figure; imagesc(E); daspect([1 1 1]); title([num2str(n) '-' num2str(j)]);

            E = full(E);
            nlen = size(E,1);
            L = sum(E,'all'); % total number of links
            dens = L / (nlen*(nlen-1));
            disp([num2str(n) '-' num2str(j) ') ' num2str(nlen) ' neurons, ' num2str(L) ' connections, density=' num2str(dens)]);

            [R, gR, count, S2] = calcReciprocity(E);
            disp([num2str(n) '-' num2str(j) ') ' num2str(gR)]);

            [coeff, tri, alltri] = calcGlobalClusteringCoeff(E, true);
            disp([num2str(n) '-' num2str(j) ') ' '(' num2str(coeff) ', ' num2str(tri) ', ' num2str(alltri) ')']);
        end
    end
end

function c = numAllTriples(A)
    c = 0;
    A = A | A'; % symmetric or
    I = 1:size(A,2); % row index array
    for i=1:size(A,1)
        idx=I(A(i,:)); % logical, a bit faster than find.
        if length(idx)==2, c=c+1;  % one triple
        elseif length(idx)>2, c=c+nchoosek(length(idx),2); % several triples
        end
    end
end

