% analyze and plot SC and FC relation.
% this script can run after analyzeStructConnectivity.m

function plotStructConnectivity
    % check ROI volume 
    % roitype: hemiroi52, CmKm50, DistKm50 (figure.1)
    % roitype: CmKm & DistKm 20 to 10000 (figure.2)
%    showSCroiVolume();
    
    % check SC post synapse cloud (figure.3)
    % roitype: hemiroi primary, hemiDistKm500
    showSCpostSynapse();

    % check SC matrix connection similarity FlyEM vs. FlyWire
    % matrix diff is used ext fig.3-1
    % roitype: hemiroi primary, hemiDistKm500
%    showSCmatrixSimilarity();

    % check SC vs. FC matrix connection FlyEM vs. FlyWire (figure.3)
    % roitype: hemiroi primary, hemiDistKm500
%    showSCFCmatrixSimilarity();

    % check neural transmitter type in each neuron (figure.3)
    % roitype: hemiroi primary, hemiDistKm500
    showNeuralTransmitter();

    % check mushroom body huge GABA neuron (APL-R) FlyEM, FlyWire
%    showAPLneuron();

    % check SC matrix connection count diff FlyEM 0sr80 vs. FlyWire 0sr140 (roi 20 to 1000)
    % roitype: Branson,Cm,DistKm
%    showSCdiffConnectionCount(); % no use

    % hemibrain ROI check other piece (Orphan) body type check.
%{
    load('data/hemiroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    checkOtherPieceSynapse('results/sc/hemiroi_connectlist.mat', ids, labelNames); % hemibrain ROI
%}
end

function showSCroiVolume()
    % volume size of one voxel
    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    volsz = prod(info.PixelDimensions*1000);

    % load SC & atlas
%    roitypes = {'hemiroi','hemiCmkm50','hemiDistKm50'};
%    maxroi = 52;
    roitypes = {'hemiCmkm20','hemiCmkm50','hemiCmkm100','hemiCmkm500','hemiCmkm1000','hemiCmkm5000','hemiCmkm10000',...
        'hemiDistKm20','hemiDistKm50','hemiDistKm100','hemiDistKm500','hemiDistKm1000','hemiDistKm5000','hemiDistKm10000'};
    maxroi = 10000;

    csvM = [];
    roiV = nan(maxroi,length(roitypes));
    for i = 1:length(roitypes)
        roitype = roitypes{i};

        % load structural connectivity matrix (from makeStructConnectivity.m)
        switch(roitype)
        case 'hemiroi'
            [roiIdxs, sz] = getHemiroiRoiIdxs();
        otherwise
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' roitype 'atlasCal.nii.gz']);
        end
        fname = ['results/sc/' lower(roitype) '_connectlist.mat'];
        load(fname);

        for k=1:length(primaryIds)
            r = primaryIds(k);
            roiV(k,i) = length(roiIdxs{r}) * volsz;
            csvM = [csvM; i, roiV(k,i)];
        end
    end
    figure; boxplot(roiV,'Labels',roitypes); title('ROI volumes in each ROI type');
    % copy & paste csvM to swarmplot excel files.
end


function [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(fname)
    atlV = niftiread(fname);
    roimax = max(atlV(:));
    sz = size(atlV);

    roiIdxs = {};
    for i=1:roimax
        roiIdxs{i} = find(atlV==i);
    end

    primaryIds = 1:roimax;
    roiNum = length(primaryIds);
end

function [roiIdxs, sz] = getHemiroiRoiIdxs()
    roiIdxs = {};
    listing = dir(['atlas/hemiroi/*.nii.gz']);
    for i=1:length(listing)
        V = niftiread(['atlas/hemiroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
        idx = find(V>0);
        roiIdxs{i} = idx;
        sz = size(V);
    end
end

function showSCpostSynapse()
    % primary, R/L, name order
    load('data/hemiroi.mat');

    % load SC & atlas
    roitypes = {'hemiroi'}; %, 'hemidistkm500'};

    for i = 1:length(roitypes)
        roitype = roitypes{i};

        fname = ['results/sc/' roitype '_postsyncount.mat'];
        load(fname);
        ids = primaryIds;
            
        % plot bar total post-synapse count
        labels = {}; X = []; Y = [];
        for c=1:length(hbsynThs)
            sth = hbsynThs(c);
            for r=1:length(confThs)
                rth = confThs(r);
                labels{end+1} = ['FlyEM' num2str(sth) 'sr' num2str(rth)];
                X(end+1) = sum(hbS(ids,r,c),'all');
                Y = [Y, hbS(ids,r,c)];
            end
        end
        for c=1:length(fwsynThs)
            sth = fwsynThs(c);
            for r=1:length(scoreThs)
                rth = scoreThs(r);
                labels{end+1} = ['FlyWire' num2str(sth) 'sr' num2str(rth)];
                X(end+1) = sum(fwS(ids,r,c),'all');
                Y = [Y, fwS(ids,r,c)];
            end
        end

        cats=categorical(labels,labels);
        figure; bar(cats,X); %set(gca,'yscale','log');
        title(['total post-synapse count in all ROI : ' roitype]);

        % plot bar mean post-synapse count of ROIs (FlyEM80-FlyWire140 case)
        p = {[4 5+5],[1 5+1]};
        figure; boxplot(Y); ylim([0 1.5e6]); %set(gca,'yscale','log');
        pval = ranksum(Y(:,p{1}(1)),Y(:,p{1}(2))); xticklabels(labels);
        title(['boxplot post-synapse count of ROIs : ' roitype ' ' labels{p{1}(1)} '-' labels{p{1}(2)} ' p=' num2str(pval)]);

        % plot bar in each ROI (FlyEM80-FlyWire140 case)
        str = split(roitype,'_');
        switch(str{1})
        case 'hemiroi'
            labelNames = roiname(ids,1);
        otherwise
            labelNames = {};
            for j=1:size(Y,1), labelNames{end+1} = num2str(j); end
        end
        cats=categorical(labelNames, labelNames);
        for j=1:length(p)
            figure; h=bar(cats,Y(:,p{j}(1))); h.FaceAlpha = 0.4;
            hold on; h=bar(cats,Y(:,p{j}(2))); h.FaceAlpha = 0.4; hold off;
            legend({labels{p{j}(1)}, labels{p{j}(2)}}); set(gca,'yscale','log');
            title(['valid post-synapse count in each ROI : ' roitype]);
        end

        % plot scatter plot by each ROI (FlyEM80-FlyWire140 case)
        if i==1
            r1 = corr(Y(:,p{1}(1)),Y(:,p{1}(2)));
            r2 = corr(log10(Y(:,p{1}(1))),log10(Y(:,p{1}(2))));
            r3 = corr(Y(:,p{1}(1)),Y(:,p{1}(2)),'type','Spearman');
            figure; scatter(Y(:,p{1}(1)),Y(:,p{1}(2)),14,'red','filled');
            hold on; plot([0 16e5], [0 16e5],':','Color',[0.5 0.5 0.5]); hold off; daspect([1 1 1]);
            text(Y(:,p{1}(1)),Y(:,p{1}(2)), labelNames, 'Vert','top', 'Horiz','left', 'FontSize',11);
            title(['EM vs. Wire valid post-synapse count in each ROI : ' roitype ' r=' num2str(r1)]);
            xlabel('FlyEM'); ylabel('FlyWire');
        end
    end
end

function showSCmatrixSimilarity()
    rgbs = [107 41 147; 55 41 185; 0 0 0; 192 0 0; 254 254 41];
    gradmap = colormapGen(rgbs,[0,0.25,0.5,0.75,1],256);

    hbsynThs = [0];
    hbrateThs = [50 60 70 80 90];
    fwsynThs = [0];
    fwrateThs = [50 70 100 130 140 150];

    % load SC & atlas
    roitypes = {'hemiroi'}; %, 'hemidistkm500'};

    for i = 1:length(roitypes)
        roitype = roitypes{i};
        fname = ['results/sc/' roitype '_postsyncount.mat'];
        load(fname);

        switch(roitype)
        case 'hemiroi'
            % primary, R/L, name order
            load('data/hemiroi.mat');
        otherwise
            roiname = {};
            for j=1:length(primaryIds), roiname{end+1,1} = num2str(j); end
        end
        ids = primaryIds;

        idlen = length(ids);
        zsz = length(hbrateThs)*length(hbsynThs)+length(fwrateThs)*length(fwsynThs);
        ncountMat = zeros(idlen,idlen,zsz);
        sycountMat = zeros(idlen,idlen,zsz);
        nweightMat = zeros(idlen,idlen,zsz);
        syweightMat = zeros(idlen,idlen,zsz);

        % check connection matrix
        ii = 1; labels = {};
        for c=1:length(hbsynThs)
            sth = hbsynThs(c);
            for r=1:length(hbrateThs)
                rth = hbrateThs(r);
                fname = ['results/sc/' roitype '_hb' num2str(sth) 'sr' num2str(rth) '_connectlist.mat'];
                if exist(fname,'file')
                    t = load(fname);
                    ncountMat(:,:,ii) = t.ncountMat(ids,ids,2);
                    sycountMat(:,:,ii) = t.sycountMat(ids,ids,2);
                    nweightMat(:,:,ii) = t.nweightMat(ids,ids,2);
                    syweightMat(:,:,ii) = t.syweightMat(ids,ids,2);
                    labels{ii} = ['FlyEM' num2str(sth) 'sr' num2str(rth)];
                    ii = ii + 1;
                end
            end
        end
        for c=1:length(fwsynThs)
            sth = fwsynThs(c);
            for r=1:length(fwrateThs)
                rth = fwrateThs(r);
                fname = ['results/sc/' roitype '_fw' num2str(sth) 'sr' num2str(rth) '_connectlist.mat'];
                if exist(fname,'file')
                    t = load(fname);
                    ncountMat(:,:,ii) = t.ncountMat(ids,ids,2);
                    sycountMat(:,:,ii) = t.sycountMat(ids,ids,2);
                    nweightMat(:,:,ii) = t.nweightMat(ids,ids,2);
                    syweightMat(:,:,ii) = t.syweightMat(ids,ids,2);
                    labels{ii} = ['FlyWire' num2str(sth) 'sr' num2str(rth)];
                    ii = ii + 1;
                end
            end
        end
        ii = ii - 1;

        for j=1:ii
            N1 = ncountMat(:,:,j);
            S1 = sycountMat(:,:,j);
            Nw1 = nweightMat(:,:,j); Nw1(isnan(Nw1)) = 0;
            Sw1 = syweightMat(:,:,j); Sw1(isnan(Sw1)) = 0;
            E = logical(eye(size(N1,1)));
            mN(j) = mean(N1(:)); mS(j) = mean(S1(:));
            mNin(j) = mean(N1(E),'all'); mSin(j) = mean(S1(E),'all');
            mNoth(j) = mean(N1(~E),'all'); mSoth(j) = mean(S1(~E),'all');

            for k=j+1:ii
                N2 = ncountMat(:,:,k);
                S2 = sycountMat(:,:,k);
                Nw2 = nweightMat(:,:,k); Nw2(isnan(Nw2)) = 0;
                Sw2 = syweightMat(:,:,k); Sw2(isnan(Sw2)) = 0;
                NSims(j,k) = getCosSimilarity(N1, N2);
                SSims(j,k) = getCosSimilarity(S1, S2);
                NwSims(j,k) = getCosSimilarity(Nw1, Nw2);
                SwSims(j,k) = getCosSimilarity(Sw1, Sw2);
                Nr(j,k) = corr(N1(:),N2(:));
                Sr(j,k) = corr(S1(:),S2(:));
                Nwr(j,k) = corr(Nw1(:),Nw2(:));
                Swr(j,k) = corr(Sw1(:),Sw2(:));

                % plot SC matrix
                if (j==4 && k==5+5) % FlyEM80-FlyWire140 case
                    labelNames = roiname(ids,1);
                    lN1 = log10(N1); lN1(isinf(lN1)) = 0; lN2 = log10(N2); lN2(isinf(lN2)) = 0; 
%                    figure; imagesc(lN1); colorbar; daspect([1 1 1]); title([num2str(j) '-' num2str(k) ' ' roitype ' neuron']);
%                    figure; imagesc(lN2); colorbar; daspect([1 1 1]); title([num2str(j) '-' num2str(k) ' ' roitype '\_fw neuron']);
                    figure; imagescLabel(N1-N2,labelNames,[-1000 1000], [labels{j} '-' labels{k} ' ' roitype ' Hemi-Wire neuron diff']); colormap(gradmap);
                    figure; imagescLabel(S1-S2,labelNames,[-100000 100000], [labels{j} '-' labels{k} ' ' roitype ' Hemi-Wire synapse diff']); colormap(gradmap);
                    figure; imagescLabel(Nw1-Nw2,labelNames,[-0.5 0.5], [labels{j} '-' labels{k} ' ' roitype ' Hemi-Wire nweight diff']); colormap(gradmap);
                    figure; imagescLabel(Sw1-Sw2,labelNames,[-0.5 0.5], [labels{j} '-' labels{k} ' ' roitype ' Hemi-Wire syweight diff']); colormap(gradmap);

                    % scatter plot of connected neuron count
                    m = max([N1(:); N2(:)]);
                    figure; scatter(N1(:),N2(:)); ylim([0 m]); xlim([0 m]); daspect([1 1 1]); set(gca,'xscale','log'); set(gca,'yscale','log');
                    hold on; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected neuron count (FlyEM)'); ylabel('connected neuron count (FlyWire)');
                    title([labels{j} '-' labels{k} ' ' roitype ' r=' num2str(Nr(j,k))]);

                    figure; scatter(Nw1(:),Nw2(:)); ylim([0 1]); xlim([0 1]); daspect([1 1 1]);
                    hold on; plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected neuron weight (FlyEM)'); ylabel('connected neuron weight (FlyWire)');
                    title([labels{j} '-' labels{k} ' ' roitype ' r=' num2str(Nwr(j,k))]);

                    % scatter plot of connected post-synapse count
                    m = max([S1(:); S2(:)]);
                    figure; scatter(S1(:),S2(:)); ylim([0 m]); xlim([0 m]); daspect([1 1 1]); set(gca,'xscale','log'); set(gca,'yscale','log');
                    hold on; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected synapse count (FlyEM)'); ylabel('connected synapse count (FlyWire)');
                    title([labels{j} '-' labels{k} ' ' roitype ' r=' num2str(Sr(j,k))]);

                    figure; scatter(Sw1(:),Sw2(:)); ylim([0 1]); xlim([0 1]); daspect([1 1 1]);
                    hold on; plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected synapse weight (FlyEM)'); ylabel('connected synapse weight (FlyWire)');
                    title([labels{j} '-' labels{k} ' ' roitype ' r=' num2str(Swr(j,k))]);
                end
            end
        end
        figure; imagescLabel(Nr,labels,[0.8 1],'neuron count matrix similarity');
        figure; imagescLabel(Sr,labels,[0.8 1],'synapse count matrix similarity');
        figure; imagescLabel(Nwr,labels,[0.8 1],'neuron weight matrix similarity');
        figure; imagescLabel(Swr,labels,[0.8 1],'synapse weight matrix similarity');
    end
end

function showSCFCmatrixSimilarity()
    preproc = 'ar'; % for move correct, slice time correct
    smoothSet = {{'', 's30', 's80'},{'', 's30', 's80','s150','s230','s300'}};
    nuisanceSet = {{'','poltcomp'},{'','poltcomp'}};
    combinationIdx = [4, 11]; % using combination index for bar plot (4=no smooth & poltcomp)
    barYlimSet = {[0.5 0.8; 0.7 1; 1.0 1.6], [0.5 0.8; 0.6 0.9; 1.0 1.5]};

    % roytypes set of SC & atlas 
    roitypesSet = {{'hemiroi_hb0sr50','hemiroi_hb0sr60','hemiroi_hb0sr70','hemiroi_hb0sr80','hemiroi_hb0sr90', ...
            'hemiroi_fw0sr50','hemiroi_fw0sr70','hemiroi_fw0sr100','hemiroi_fw0sr130','hemiroi_fw0sr140','hemiroi_fw0sr150', ...
            'hemiroi_fw5sr50','hemiroi_fw5sr140'}, ...
            {'hemiDistKm500_hb0sr50','hemiDistKm500_hb0sr60','hemiDistKm500_hb0sr70','hemiDistKm500_hb0sr80','hemiDistKm500_hb0sr90', ...
            'hemiDistKm500_fw0sr50','hemiDistKm500_fw0sr70','hemiDistKm500_fw0sr100','hemiDistKm500_fw0sr130','hemiDistKm500_fw0sr140','hemiDistKm500_fw0sr150', ...
            'hemiDistKm500_fw5sr50','hemiDistKm500_fw5sr140'}};

    for i=1:length(smoothSet)
        showSCFCmatrixSimilaritySet(preproc, smoothSet{i}, nuisanceSet{i}, roitypesSet{i}, combinationIdx(i), barYlimSet{i})
    end
end

function showSCFCmatrixSimilaritySet(preproc, smooth, nuisance, roitypes, cIdx, barYlim)
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

    ylabels = {}; rlabels = {}; R3 = []; A3 = []; roiR3 = [];
    for r = 1:length(roitypes)
        str = split(roitypes{r},'_'); rlabels{r} = str{2};
        Am = []; Rm = []; roiRm = []; ii=1; xlabels = {};
        for n=1:length(nuisance)
            for k=1:length(smooth)
                pftype = [smooth{k} nuisance{n} preproc roitypes{r}];
                xlabels{ii} = [smooth{k} nuisance{n}]; ii=ii+1;
                aucmat = ['results/auc/' pftype '-fcauc.mat'];
                if exist(aucmat,'file')
                    load(aucmat);
                else
                    A = nan(size(Am,1),100);
                    R = nan(size(Rm,1),1);
                    roiR = nan(size(roiRm,1),24);
                end
                Am = [Am,nanmean(A,2)];
                Rm = [Rm,R(:)];
                roiRm = cat(3,roiRm,roiR);
            end
        end

        R3 = [R3;Rm]; A3 = [A3;Am]; roiR3 = cat(2,roiR3,roiRm);
        C = cell(24,1); C(1:24) = {[str{2} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end
    % FC-SC correlation (all)
    T = [0:24:24*12];
    I = getR3idx([7 9],T);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],T);
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9],T);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9],T);
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9],T);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],T);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});

    % plot SC vs. poltcomp m-FC(z) result
    I1 = getR3idx([7],T); I2 = getR3idx([9],T); I3 = getR3idx([10],T); I4 = getR3idx([11],T);
    llabels = {'neuron count','synapse count','neuron weight','synapse weight'};
    cats=categorical(rlabels,rlabels);

    figure; bar(cats,[R3(I1,cIdx)';R3(I2,cIdx)';R3(I3,cIdx)';R3(I4,cIdx)']); ylim(barYlim(1,:));
    title([str{1} ' FC-SC correlation : SC vs. poltcomp m-FC(z)']); legend(llabels);

    figure; bar(cats,[A3(I1,cIdx)';A3(I2,cIdx)';A3(I3,cIdx)';A3(I4,cIdx)']); ylim(barYlim(2,:));
    title([str{1} ' FC-SC detection : SC vs. poltcomp m-FC(z)']); legend(llabels);

    figure; bar(cats,[B(I1,cIdx)';B(I2,cIdx)';B(I3,cIdx)';B(I4,cIdx)']); ylim(barYlim(3,:));
    title([str{1} ' FC-SC correlation & detection : SC vs. poltcomp m-FC(z)']); legend(llabels);

    % plot ROI SC vs. poltcomp m-FC(z) result
    switch(str{1})
    case 'hemiroi'
        fname = ['results/sc/' str{1} '_postsyncount.mat'];
        load(fname);
        ids = primaryIds;
        % primary, R/L, name order
        load('data/hemiroi.mat');
        labelNames = roiname(ids,1);
    otherwise
        labelNames = {};
        for j=1:size(roiR3,1), labelNames{end+1} = num2str(j); end
    end
    cats=categorical(labelNames, labelNames);
    figure; h=bar(cats,roiR3(:,72+7,cIdx)); h.FaceAlpha = 0.4;
    hold on; h=bar(cats,roiR3(:,216+7,cIdx)); h.FaceAlpha = 0.4; hold off;
    title(['FC-SC correlation in each ROI']); legend({ylabels{72+7},ylabels{216+7}});

    figure; scatter(roiR3(:,72+7,cIdx),roiR3(:,216+7,cIdx)); ylim([0 1]); xlim([0 1]); daspect([1 1 1]);
    hold on; plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]); hold off; xlabel(ylabels{72+7}); ylabel(ylabels{216+7});
    r = corr(roiR3(:,72+7,cIdx),roiR3(:,216+7,cIdx));
    title(['correlation of regional FC-SC between FlyEM and FlyWire r=' num2str(r)]);
end


function showNeuralTransmitter()
    llabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    rateThs = [140]; %50 70 100  140 150];
    synThs = [0];
    roitypes = {'hemiroi'}; %, 'hemidistkm500'};

    % use appropriate smooth and nuisance
    preproc = 'ar';
    smooth = {'','s230'};
    nuisance = {'poltcomp','poltcomp'};

    for n = 1:length(roitypes)
        roitype = roitypes{n};
        for rt=1:length(rateThs)
            rth = rateThs(rt);
            for j=1:length(synThs)
                sth = synThs(j);
                idstr = [roitype '_fw' num2str(sth) 'sr' num2str(rth)];
                fname = ['results/sc/' idstr '_transmitter.mat'];
                load(fname);

                fname = ['results/sc/' idstr '_connectlist.mat'];
                load(fname);
                ids = primaryIds;

                switch(roitype)
                case 'hemiroi'
                    % primary, R/L, name order
                    load('data/hemiroi.mat');
                    labelNames = roiname(ids,1);
                otherwise
                    labelNames = {};
                    for k=1:length(primaryIds), labelNames{end+1} = num2str(k); end
                end
                cats=categorical(labelNames,labelNames);

                % plot just transmitter rate in each ROI.
                % because neuron number in each ROI is heavly different, so
                % absolute count is not so infomative.
                outNTypes = outNTypes ./ sum(outNTypes,2);
                figure; bar(cats,outNTypes(ids,:)','stacked'); ylim([0 1]);
                title([idstr ' out-neuron transmitters in each ROI']); legend(llabels);

                inNTypes = inNTypes ./ sum(inNTypes,2);
                figure; bar(cats,inNTypes(ids,:)','stacked'); ylim([0 1]);
                title([idstr ' in-neuron transmitters in each ROI']); legend(llabels);

                inSyTypes = inSyTypes ./ sum(inSyTypes,2);
                figure; bar(cats,inSyTypes(ids,:)','stacked'); ylim([0 1]);
                title([idstr ' in-synapse transmitters in each ROI']); legend(llabels);

                % compare with FC-SC correlation in each ROI.
                pftype = [smooth{n} nuisance{n} preproc idstr];
                aucmat = ['results/auc/' pftype '-fcauc.mat'];
                load(aucmat);
                FcSc = roiR(:,7); % already primaryId. (neuron count vs. m-FC(z))
                figure; sgtitle([pftype ' FC-SC vs. out-neuron ']);
                for i=2:7
                    [r,p] = corr(FcSc,outNTypes(ids,i)); 
                    disp([pftype ' FC-SC vs. out-neuron ' llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                    subplot(3,2,i-1);
                    scatter(FcSc,outNTypes(ids,i)); title([llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                end
                figure; sgtitle([pftype ' FC-SC vs. in-neuron ']);
                for i=2:7
                    [r,p] = corr(FcSc,inNTypes(ids,i)); 
                    disp([pftype ' FC-SC vs. in-neuron ' llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                    subplot(3,2,i-1);
                    scatter(FcSc,inNTypes(ids,i)); title([llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                end
                figure; sgtitle([pftype ' FC-SC vs. in-synapse ']);
                for i=2:7
                    [r,p] = corr(FcSc,inSyTypes(ids,i)); 
                    disp([pftype ' FC-SC vs. in-synapse ' llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                    subplot(3,2,i-1);
                    scatter(FcSc,inSyTypes(ids,i)); title([llabels{i} ' r=' num2str(r) ' p=' num2str(p)]);
                end
            end
        end
    end
end

function showSCdiffConnectionCount()
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {'hemiBranson7065km','hemiCmkm','hemiDistKm'};

    NSims = nan(length(roitypes),length(roinums),'single');
    SSims = NSims; xlabels = {}; labels = {}; axlabels = {}; 
    mN1 = NSims; mN2 = mN1; mS1 = mN1; mS2 = mN1;
    mN1in = mN1; mN2in = mN1; mS1in = mN1; mS2in = mN1;
    mN1oth = mN1; mN2oth = mN1; mS1oth = mN1; mS2oth = mN1;
    Vs = [];

    % load SC & atlas
    for i = 1:length(roitypes)
        labels{i} = roitypes{i}; labels{i+length(roitypes)} = [roitypes{i} 'Fw'];

        for j=1:length(roinums)
            % check atlas voxel size
            ii = j+(i-1)*length(roinums);
            axlabels{ii} = [roitypes{i} num2str(roinums(j))];
            info = niftiinfo(['atlas/' axlabels{ii} 'atlasCal.nii.gz']);
            V = niftiread(info);
            Vs(ii) = length(find(V>0));

            % check connection matrix
            xlabels{j} = ['roi' num2str(roinums(j))];
            roitype = [roitypes{i} num2str(roinums(j))];
            confile = ['results/sc/' lower(roitype) '_connectlist.mat'];
            confilefw = ['results/sc/' lower(roitype) '_fw0sr140_connectlist.mat'];
            if exist(confile,'file') && exist(confilefw,'file')
                t = load(confile);
                ids = t.primaryIds;
                N1 = t.ncountMat(ids,ids,2); S1 = t.sycountMat(ids,ids,2);
                E = logical(eye(roinums(j)));
                mN1(i,j) = mean(N1(:)); mS1(i,j) = mean(S1(:));
                mN1in(i,j) = mean(N1(E),'all'); mS1in(i,j) = mean(S1(E),'all');
                mN1oth(i,j) = mean(N1(~E),'all'); mS1oth(i,j) = mean(S1(~E),'all');
                clear t;
    
                t = load(confilefw);
                ids = t.primaryIds;
                N2 = t.ncountMat(ids,ids,2); S2 = t.sycountMat(ids,ids,2);
                mN2(i,j) = mean(N2(:)); mS2(i,j) = mean(S2(:));
                mN2in(i,j) = mean(N2(E),'all'); mS2in(i,j) = mean(S2(E),'all');
                mN2oth(i,j) = mean(N2(~E),'all'); mS2oth(i,j) = mean(S2(~E),'all');

                NSims(i,j) = getCosSimilarity(N1, N2);
                SSims(i,j) = getCosSimilarity(S1, S2);
%{
                % scatter plot of connected neuron count
                m = max([N1(:); N2(:)]);
                figure; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); ylim([0 m]); xlim([0 m]); daspect([1 1 1]);
                hold on; scatter(N1(:),N2(:)); hold off; xlabel('connected neuron count (FlyEM)'); ylabel('connected neuron count (FlyWire)');
                title(roitype);

                % scatter plot of connected post-synapse count
                m = max([S1(:); S2(:)]);
                figure; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); ylim([0 m]); xlim([0 m]); daspect([1 1 1]);
                hold on; scatter(S1(:),S2(:)); hold off; xlabel('connected synapse count (FlyEM)'); ylabel('connected synapse count (FlyWire)');
                title(roitype);
%}
            end
        end
    end

    % plot atlas total ROI voxels
%    cats=categorical(axlabels);
%    figure; bar(cats, Vs); title('atlas total ROI voxels');

    % plot cosine similarity result between FlyEM and FlyWire matrices
    figure; imagescLabel2(NSims,xlabels,roitypes,[0 1]); colorbar; colormap(hot); title(['Similarity of neuron count matrices between FlyEM and FlyWire.'])
    figure; imagescLabel2(SSims,xlabels,roitypes,[0 1]); colorbar; colormap(hot); title(['Similarity of synapse count matrices between FlyEM and FlyWire.'])

    % plot mean connected neuron & post-synapse num
    figure; plot([mN1' mN2']); legend(labels); setlineColors(3); title(['mean connected neuron count between FlyEM and FlyWire.'])
    figure; plot([mS1' mS2']); legend(labels); setlineColors(3); title(['mean connected synapse count between FlyEM and FlyWire.'])

    figure; plot([mN1in' mN2in']); legend(labels); setlineColors(3); title(['mean inside connected neuron count between FlyEM and FlyWire.'])
    figure; plot([mS1in' mS2in']); legend(labels); setlineColors(3); title(['mean inside connected synapse count between FlyEM and FlyWire.'])

    figure; plot([mN1oth' mN2oth']); legend(labels); setlineColors(3); title(['mean other connected neuron count between FlyEM and FlyWire.'])
    figure; plot([mS1oth' mS2oth']); legend(labels); setlineColors(3); title(['mean other connected synapse count between FlyEM and FlyWire.'])
end

function showAPLneuron()
    % load mushroom body ROIs.
%    aV = niftiread('atlas/hemiRoi68-59-87-106-50-27-54atlasCal.nii.gz');
%    idx = find(aV>0);

    % load connected post-synapse counts from APL pre-synapses (output)
    hV = niftiread('results/nifti/hemiAplOutputSynapses.nii.gz');
    wV = niftiread('results/nifti/wireAplOutputSynapses.nii.gz');

    % calculate Sørensen–Dice index
    bV = hV + wV;
    idx = find(bV>0);
    cV = hV .* wV;
    idx2 = find(cV>0);
    Didx = 2*length(idx2) / (length(find(hV>0))+length(find(wV>0)));
%    Didx = dice(logical(hV),logical(wV)); % same value

    smooth = [0 20 40 80];
    for i=1:length(smooth)
        sz = smooth(i);
        if sz > 0
           FWHM = [(sz/10)/(2.45/2.28) sz/10 (sz/10)/(3.715/2.28)]; % voxel size;
            sigma = FWHM / sqrt(8*log(2));
            filterSize = 2*ceil(2*sigma)+1;
                
            disp(['sz=' num2str(sz) ', sigma=' num2str(sigma(1)) ', flSz=' num2str(filterSize(1))]);
            hVt = imgaussfilt3(hV, sigma, 'FilterSize', filterSize);
            wVt = imgaussfilt3(wV, sigma, 'FilterSize', filterSize);
        else
            hVt = hV; wVt = wV;
        end

        hscnt = hVt(idx);
        wscnt = wVt(idx);
        mcnt = max([hscnt; wscnt]);
        r = corr(hscnt,wscnt);    % corr between FlyEM vs. FlyWire post-synapses to connect APL-R neuron
        figure; scatter(hscnt,wscnt,12,'x'); xlabel('FlyEM synapse count'); ylabel('FlyWire synapse count');
        hold on; plot([0 mcnt], [0 mcnt],':','Color',[0.5 0.5 0.5]); hold off; daspect([1 1 1]);
        title(['s' num2str(sz) ' corr between synapses FlyEM vs. FlyWire r=' num2str(r)]);
    end
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
    figure; imagescLabel(CM2logi.*E, labelNames, [], 'CM2 logi'); % ignore diag
    figure; imagescLabel(CM3.*E, labelNames, [], 'other count'); % ignore diag
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

function imagescLabel(mat, labelNames, range, titlestr)
    if isempty(range)
        imagesc(mat);
    else
        imagesc(mat, range);
    end
    colorbar; daspect([1 1 1]); title(titlestr);
    set(gca,'XTick',1:size(mat,2));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end

