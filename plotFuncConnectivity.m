% analyze and plot SC and FC relation.

function plotFuncConnectivity
    vslabels = {
        'log10(neurons f) vs. m-FCz', ...
        'log10(synapse weight f) vs. m-FCz', ...
        'log10(synapses f) vs. m-FCz', ...
        'ROI in-neuron weight f vs. m-FCz', ...
        'ROI in-synapse weight f vs. m-FCz', ...
        'ROI out-neuron weight f vs. m-FCz', ...
        'log10(neurons) vs. m-FCz', ... %7
        'log10(synapse weight) vs. m-FCz', ...
        'log10(synapses) vs. m-FCz', ...
        'ROI in-neuron weight vs. m-FCz', ...
        'ROI in-synapse weight vs. m-FCz', ...
        'ROI out-neuron weight vs. m-FCz', ...
        'log10(neurons f) vs. FC-Tval', ... %13
        'log10(synapse weight f) vs. FC-Tval', ...
        'log10(synapses f) vs. FC-Tval', ...
        'ROI in-neuron weight f vs. FC-Tval', ...
        'ROI in-synapse weight f vs. FC-Tval', ...
        'ROI out-neuron weight f vs. FC-Tval', ...
        'log10(neurons) vs. FC-Tval', ... %19
        'log10(synapse weight) vs. FC-Tval', ...
        'log10(synapses) vs. FC-Tval', ...
        'ROI in-neuron weight vs. FC-Tval', ...
        'ROI in-synapse weight vs. FC-Tval', ...
        'ROI out-neuron weight vs. FC-Tval', ...
    };

    % check smoothing result around 50 ROIs
%    checkSmoothingResult50(vslabels);

    % check correlation result in each ROI num
%    checkSmoothingByRoinum(vslabels)

    % check nuisance result round 50 ROIs
    checkNuisanceResult50(vslabels);

    % hemibrain ROI check other piece (Orphan) body type check.
%{
    load('data/flyemroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    checkOtherPieceSynapse('data/neuprint_connectlist.mat', ids, labelNames); % hemibrain ROI
%}
end

function checkSmoothingByRoinum(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {{'hemiBranson7065km',''},{'hemiCmkm',''},{'hemiCmkm','r1w1'},{'hemiDistKm',''},{'hemiRand',''},{'hemiVrand',''}}; %,{'hemiCmkm','r2w1'}};
    roitypelabels = {'Branson','Cm','CmR1w1','Dist','Rand','Vand'}; %,'CmR2w1'};
    
    ylabels = {}; R3 = []; A3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; ii=1; xlabels = {};
        for rr=1:length(roinums)
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                for k=1:length(smooth)
                    for n=1:length(nuisance)
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(roinums(rr)) roitypes{r}{2}];
                        xlabels{ii} = [smooth{k} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/' pftype '-fcauc.mat'];
                        if exist(aucmat,'file')
                            load(aucmat);
                        else
                            A = nan(size(Am,1),100);
                            R = nan(size(Rm,1),1);
                        end
                        Am = [Am,nanmean(A,2)];
                        Rm = [Rm,R(:)];
                    end
                end
            end
        end

        R3 = [R3;Rm]; A3 = [A3;Am];
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24 48 72 96 120]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48 72 96 120]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120]);
    figure; plot(B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSmoothingResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roitypes = {'flyemroi','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
    roitypelabels = {'FlyEM','Branson','Cm','CmR1w1','Dist','Rand','Vrand'};

    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = [];
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    aucmat = ['results/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(Am,1),100);
                        R = nan(size(Rm,1),1);
                    end
                    AA = cat(3,AA,A);
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end

        R3 = [R3;Rm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA);
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end


    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),smooth,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([0:8]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);
    
    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),smooth,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([0:8]',A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2);

%    AA = squeeze(AA3(20,:,:));
%    figure; boxplot(AA,'Labels',smooth); title(['FC-SC detection ' ylabels{20}]);

    for i=[7 9]
        X = []; slabels = {};
        for j=1:length(roitypelabels)
            Y = squeeze(AA3(i+(j-1)*24,:,:)); Y(Y==0)=nan; % ignore zero result
            X = [X, Y];
            C = cell(9,1); C(1:9) = {[roitypelabels{j} ' ']};
            slabels = [slabels(:); strcat(C(:),smooth(:))];
        end
        figure; plot(X); legend(slabels); title(['FC-SC detection results by threshold in (' num2str(i) ') ' vslabels{i}]); setlineColors(9);
    end

    % both FC-SC correlation & detection
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),smooth,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([0:8]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);
end

function checkNuisanceResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
    roitypes = {'flyemroi','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
    roitypelabels = {'FlyEM','Branson','Cm','CmR1w1','Dist','Rand','Vrand',};

    ylabels = {}; R3 = []; A3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    aucmat = ['results/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(Am,1),100);
                        R = nan(size(Rm,1),1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end

        R3 = [R3;Rm]; A3 = [A3;Am];
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),nuisance,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([1:37]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),nuisance,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([1:37]',A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2);

    % both FC-SC correlation & detection
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),nuisance,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120 144]);
    figure; plot([1:37]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);
end

function I = getR3idx(A,B)
    I = repmat(A',[1 length(B)]) + repmat(B,[length(A) 1]); I=I(:);
end

function imagescLabel2(mat, xlabel, ylabel, range)
    if isempty(range)
        imagesc(mat);
    else
        imagesc(mat, range);
    end
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

function setlineStyles(styles)
    ax = gca; ax.LineStyleOrder = styles; % this is for R2023b
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

