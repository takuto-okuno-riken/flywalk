% analyze and plot SC and FC relation.

function plotFuncConnectivity
    vslabels = {
        'log10(cell count2) vs. FC(z)', ...
        'log10(synapse weight2) vs. FC(z)', ...
        'log10(synapse count) vs. FC(z)', ...
        'ROI in-neuron weight vs. FC(z)', ...
        'ROI in-synapse weight vs. FC(z)', ...
        'ROI out-neuron weight vs. FC(z)', ...
        'log10(cell count2b) vs. FC(z)', ... %7
        'log10(synapse weight2b) vs. FC(z)', ...
        'log10(synapse count b) vs. FC(z)', ...
        'ROI in-neuron weight b vs. FC(z)', ...
        'ROI in-synapse weight b vs. FC(z)', ...
        'ROI out-neuron weight b vs. FC(z)', ...
        'log10(cell count2) vs. FC T-val', ... %13
        'log10(synapse weight2) vs. FC T-val', ...
        'log10(synapse count) vs. FC T-val', ...
        'ROI in-neuron weight vs. FC T-val', ...
        'ROI in-synapse weight vs. FC T-val', ...
        'ROI out-neuron weight vs. FC T-val', ...
        'log10(cell count2b) vs. FC T-val', ... %19
        'log10(synapse weight2b) vs. FC T-val', ...
        'log10(synapse count b) vs. FC T-val', ...
        'ROI in-neuron weight b vs. FC T-val', ...
        'ROI in-synapse weight b vs. FC T-val', ...
        'ROI out-neuron weight b vs. FC T-val', ...
    };

    % check smoothing result around 50 ROIs
    checkSmoothingResult50(vslabels);

    % hemibrain ROI check other piece (Orphan) body type check.
%{
    load('data/flyemroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    checkOtherPieceSynapse('data/neuprint_connectlist.mat', ids, labelNames); % hemibrain ROI
%}
end

function checkSmoothingResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roitypes = {'flyemroi','hemiBranson7065km50','hemiKm50'}; % flyem ROI (Turner compatible)

    rlabel = {}; ii=1;
    Rm = []; AUC = [];
    for r = 1:length(roitypes)
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    rlabel{ii} = pftype; ii=ii+1;
                    aucmat = ['results/' pftype '-fcauc.mat'];
                    load(aucmat);
                    AUC = cat(3,AUC,A);
                    Rm = [Rm,R(:)];
                end
            end
        end
    end

    % FC-SC correlation (4-type mixed box plot)
    figure; imagescLabel2(Rm,rlabel,vslabels); colorbar; title('FC-SC correlation');
    figure; plot(Rm'); legend(vslabels); title('FC-SC correlation'); setlineColors(6);

%    AA = squeeze(AUC(20,:,:));
%    figure; boxplot(AA,'Labels',rlabel); title([' FC-SC detection (FC T-val vs. ROI in-neuron weight b)']);

    AA = squeeze(nanmean(AUC,2));
  %  figure; boxplot(AA(5:8,:),'Labels',rlabel); title([' FC-SC detection (4-type mixed plot)']);
  %  hold on; plot(AA(5:8,:)'); hold off; legend;
    figure; imagescLabel2(AA,rlabel,vslabels); colorbar; title('FC-SC detection');
    figure; plot(AA'); legend(vslabels); title('FC-SC detection'); setlineColors(6);

    for i=[7 13]
        X = squeeze(AUC(i,:,:));
        figure; plot(X); legend(rlabel); title(['FC-SC detection results by threshold in (' num2str(i) ') ' vslabels{i}]); setlineColors(9);
    end
%{
    for ii=1:length(rlabel) % pattern
        X = squeeze(AUC(:,:,ii));
        figure; plot(X'); legend(vslabels); title(['FC-SC detection results by threshold ' rlabel{ii}]); 
        xlabel('percentile'); ylabel('AUC'); setlineColors(6);
    end
%}
    % both FC-SC correlation & detection
    B = abs(Rm) + abs(AA-0.5)*2;
    figure; imagescLabel2(B,rlabel,vslabels); colorbar; title('FC-SC correlation & detection');
    figure; plot(B'); legend(vslabels); title('FC-SC correlation & detection'); setlineColors(6);
end

function imagescLabel2(mat, xlabel, ylabel)
    imagesc(mat);
    set(gca,'XTick',1:length(xlabel)); set(gca,'XTickLabel',xlabel);
    set(gca,'YTick',1:length(ylabel)); set(gca,'YTickLabel',ylabel);
end

function setlineColors(maxcol)
    cols = zeros(3*maxcol,3);
    for i=1:3
        cols(maxcol*(i-1)+1:maxcol*i,mod(i,3)+1) = 1;
        cols(maxcol*(i-1)+1:maxcol*i,mod(i+1,3)+1) = 0:1/(maxcol-1):1;
    end
    for i=4:6
        cols(maxcol*(i-1)+1:maxcol*i,mod(i,3)+1) = 0.4;
        cols(maxcol*(i-1)+1:maxcol*i,mod(i+1,3)+1) = [0:1/(maxcol-1):1] * 0.5;
        cols(maxcol*(i-1)+1:maxcol*i,mod(i+2,3)+1) = [1/maxcol:1/maxcol:1] * 0.5;
    end
    ax = gca; ax.ColorOrder = cols;
end

function checkOtherPieceSynapse(fname, ids, labelNames)
    % load matrix info
    load([fname(1:end-4) '_cnids.mat']); % connected nid cells
    roimax = size(Cnids,1);

    countMat = nan(roimax,roimax,3,'single');
    for p=2:3
        for i=1:roimax
            for j=1:roimax
                countMat(i,j,p) = length(Cnids{i,j,p});
            end
        end
    end

    CM2 = countMat(ids,ids,2);
    CM3 = countMat(ids,ids,3);
    logi = (CM3>0);
    r = corr(CM2(:),CM3(:));
    disp(['neuron count vs. other count : ' num2str(r)]);

    lsz = length(ids);
    E = eye(lsz,lsz); E=1-E;
    CM2logi = CM2 .* logi .* E;
    CM3 = CM3 .* E;
    figure; imagescLabel(CM2logi.*E, labelNames, 'CM2 logi'); % ignore diag
    figure; imagescLabel(CM3.*E, labelNames, 'other count'); % ignore diag
    figure; scatter(CM2logi(:), CM3(:)); xlabel('CM2 logi'); ylabel('other count')

    % find top 20 other count or other rate.
    [M,idx1] = sort(CM3(:),'descend');
    for i=1:20
        [x,y] = ind2sub([lsz,lsz],idx1(i));
        disp([num2str(i) ') ' labelNames{x} '-' labelNames{y}  ' neurons:' num2str(CM2(x,y)) ' others:' num2str(CM3(x,y))]);
    end
    Gr = CM3 ./ CM2; Gr(isnan(Gr))=0; Gr(CM2<5)=0; % others connection rate compared to neuron connection (neurons>=5)
    [M,idx2] = sort(Gr(:),'descend');
    for i=1:20
        [x,y] = ind2sub([lsz,lsz],idx2(i));
        disp([num2str(i) ') ' labelNames{x} '-' labelNames{y}  ' neurons:' num2str(CM2(x,y)) ' others:' num2str(CM3(x,y)) ' rate:' num2str(Gr(x,y))]);
    end

    % read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    % read synapse info and show other synaptic info
    load('data/hemibrain_v1_2_synapses.mat');

    for i=1:20
        [x,y] = ind2sub([lsz,lsz],idx2(i));
        nids = Cnids{ids(x),ids(y),3};
        logi = ismember(StoN,nids);
        sids = find(logi==1);
        idx = find(Sdir(sids)==1); % get pre-synapse ids
        presids = sids(idx);
        idx = find(Sdir(sids)==2); % get post-synapse ids
        postsids = sids(idx);
        disp([num2str(i) ') ' labelNames{x} '-' labelNames{y}  ' cells:' num2str(length(nids)) ' pre:' num2str(length(presids)) ' post:' num2str(length(postsids))]);

        % show sample locations (sorted by large size)
        ST={'non', 'Traced', 'Orphan', 'Assign', 'Unimportant'};
        prenids = StoN(presids);
        logi = ismember(Nid,prenids);
        prenidx = find(logi==1);
        [M,I] = sort(Nsize(prenidx),'descend');
        for j=1:5
            idx = prenidx(I(j));
            logi = ismember(StoN,Nid(idx));
            S2=find(logi==1); k=find(Sdir(S2)==1);
            disp([' nid:' num2str(Nid(idx)) '  pre sidx:' num2str(S2(k(1))) ') loc: ' num2str(Sloc(S2(k(1)),:)) '   st:' ST{Nstatus(idx)+1} ' sz:' num2str(Nsize(idx))]);
        end
        postnids = StoN(postsids);
        logi = ismember(Nid,postnids);
        postnidx = find(logi==1);
        [M,I] = sort(Nsize(postnidx),'descend');
        for j=1:5
            idx = postnidx(I(j));
            logi = ismember(StoN,Nid(idx));
            S2=find(logi==1); k=find(Sdir(S2)==2);
            disp([' nid:' num2str(Nid(idx))  ' post sidx:' num2str(S2(k(1))) ' loc: ' num2str(Sloc(S2(k(1)),:)) '   st:' ST{Nstatus(idx)+1} ' sz:' num2str(Nsize(idx))]);
        end
    end
end

function imagescLabel(mat, labelNames, titlestr)
    imagesc(mat); colorbar; daspect([1 1 1]); title(titlestr);
    set(gca,'XTick',1:size(mat,1));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end

