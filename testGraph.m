function testGraph
    % Erdos-Renyi random graph
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
end
