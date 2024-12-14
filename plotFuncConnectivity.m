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

    % check smoothing result around 50 ROIs (s0 to s80)
    % roitype: FlyEM,FlyEMFw,Branson,BransonFw,Cm,CmFw,CmR1w1,Dist,Rand,Vrand
%    checkSmoothingResult50(vslabels);

    % check correlation result in each ROI num (s0 to s80, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
%    checkSmoothingByRoinum(vslabels);

    % check nuisance result round 50 ROIs (all nuisance)
    % roitype: FlyEM,FlyEmFw,Branson,Cm,CmR1w1,Dist,Rand,Vrand
%    checkNuisanceResult50(vslabels);

    % check correlation result in each ROI num (all nuisance, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
%    checkNuisanceByRoinum(vslabels);

    % check correlation result of large smoothing size (s0 to 300, roi 50 to 500, '' & poltcomp)
    % roitype: Cm,Dist
%    checkLargeSmoothingPoltcompByRoinum(vslabels);

    % check correlation result in each ROI num (roi 100 to 20000)
    % roitype: Cm,CmR1w1,Dist,Vand
%    checkNeuronVsSynapseByRoinum(vslabels);

    % check correlation result in each ROI num (s30 & s80, roi 100 to 10000, '' & poltcomp)
    % roitype: Cm,DistKm
%    checkSmoothingNuisanceByRoinum(vslabels);

    % check smoothing result of FlyEM vs. FlyWire around 50 ROIs (s0 to s80)
    % roitype: FlyEM,FlyEMFw,DistKm50,DistKm50Fw
%    checkSmoothNuisanceFlyWireResult50(vslabels);

    % check correlation result of FlyEM vs. FlyWire (s0 to s80, roi 20 to 1000)
    % roitype: Cm,CmFw,Branson,BransonFw,Dist,DistFw
    checkSmoothingFlyWireByRoinum(vslabels);

    % check SC matrix difference FlyEM vs. FlyWire (roi 20 to 1000)
    % roitype: Branson,Cm,DistKm
    checkSCdiffFlyWire(vslabels);

    % hemibrain ROI check other piece (Orphan) body type check.
%{
    load('data/flyemroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    checkOtherPieceSynapse('data/neuprint_connectlist.mat', ids, labelNames); % hemibrain ROI
%}
end

function checkSmoothingByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {{'hemiBranson7065km',''},{'hemiCmkm',''},{'hemiCmkm','r1w1'},{'hemiDistKm',''},{'hemiRand',''},{'hemiVrand',''}};
    roitypelabels = {'Branson','Cm','CmR1w1','Dist','Rand','Vand'};

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
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    roitypes = {'flyemroi','flyemroi_fw','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
    roitypelabels = {'FlyEM','FlyEMFw','Branson','Cm','CmR1w1','Dist','Rand','Vrand'};

    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = [];
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144 168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),smooth,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120]);% 144 168 192 216]);
    figure; plot([0:8]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);
    
    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144 168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),smooth,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96 120]);% 144 168 192 216]);
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
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120 144 168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),smooth,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120]);% 144 168 192 216]);
    figure; plot([0:8]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);
end

function checkSmoothNuisanceFlyWireResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roitypes = {'flyemroi','flyemroi_fw','hemiDistKm50','hemiDistKm50_fw'};
    roitypelabels = {'FlyEM','FlyEMFw','DistKm50','DistKm50Fw'};

    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = [];
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),smooth,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72]);
    figure; plot([0:8]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);
    
    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),smooth,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72]);
    figure; plot([0:8]',A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2);

    for i=[7 9]
        X = []; slabels = {};
        for j=1:length(roitypelabels)
            [m,idx] = max(A3(i+(j-1)*24,:),[],2);
            Y = squeeze(AA3(i+(j-1)*24,:,idx)); Y(Y==0)=nan; % ignore zero result
            X = [X, Y'];
            slabels = [slabels(:); [roitypelabels{j} ' ' smooth{idx}]];
        end
        figure; plot(X); legend(slabels); title(['FC-SC detection results by threshold in (' num2str(i) ') ' vslabels{i}]);
    end

    % both FC-SC correlation & detection
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),smooth,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0 24 48 72]);
    figure; plot([0:8]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);

    % checkSCdiffFlyWire
    roitypes = {'flyemroi'};

    NSims = nan(length(roitypes),1,'single');
    SSims = NSims; xlabels = {}; labels = {}; axlabels = {}; 
    mN1 = NSims; mN2 = mN1; mS1 = mN1; mS2 = mN1;
    mN1in = mN1; mN2in = mN1; mS1in = mN1; mS2in = mN1;
    mN1oth = mN1; mN2oth = mN1; mS1oth = mN1; mS2oth = mN1;
    Vs = [];

    rgbs = [107 41 147; 55 41 185; 0 0 0; 192 0 0; 254 254 41];
    gradmap = colormapGen(rgbs,[0,0.25,0.5,0.75,1],256);

    load('data/flyemroi.mat');

    % load SC & atlas
    for i = 1:length(roitypes)
        labels{i} = roitypes{i}; labels{i+length(roitypes)} = [roitypes{i} 'Fw'];

        % check connection matrix
        roitype = roitypes{i};
        confile = ['data/' lower(roitype) '_connectlist.mat'];
        confilefw = ['data/' lower(roitype) '_fw_connectlist.mat'];
        if exist(confile,'file') && exist(confilefw,'file')
            t = load(confile);
            ids = t.primaryIds;
            N1 = t.ncountMat(ids,ids,2); S1 = t.sycountMat(ids,ids,2);
            E = logical(eye(size(N1,1)));
            mN1(i) = mean(N1(:)); mS1(i) = mean(S1(:));
            mN1in(i) = mean(N1(E),'all'); mS1in(i) = mean(S1(E),'all');
            mN1oth(i) = mean(N1(~E),'all'); mS1oth(i) = mean(S1(~E),'all');
            clear t;

            t = load(confilefw);
            ids = t.primaryIds;
            N2 = t.ncountMat(ids,ids,2); S2 = t.sycountMat(ids,ids,2);
            mN2(i) = mean(N2(:)); mS2(i) = mean(S2(:));
            mN2in(i) = mean(N2(E),'all'); mS2in(i) = mean(S2(E),'all');
            mN2oth(i) = mean(N2(~E),'all'); mS2oth(i) = mean(S2(~E),'all');

            NSims(i) = getCosSimilarity(N1, N2);
            SSims(i) = getCosSimilarity(S1, S2);

            % plot SC matrix
            labelNames = roiname(ids,1);
            lN1 = log10(N1); lN1(isinf(lN1)) = 0; lN2 = log10(N2); lN2(isinf(lN2)) = 0; 
%            figure; imagesc(lN1); colorbar; daspect([1 1 1]); title([roitype ' neuron']);
%            figure; imagesc(lN2); colorbar; daspect([1 1 1]); title([roitype '\_fw neuron']);
            figure; imagescLabel(N1-N2,labelNames,[-1000 1000], [roitype ' Hemi-Wire neuron diff']); colormap(gradmap);

            figure; imagescLabel(S1-S2,labelNames,[-100000 100000], [roitype ' Hemi-Wire synapse diff']); colormap(gradmap);

            % scatter plot of connected neuron count
            m = max([N1(:); N2(:)]);
            figure; scatter(N1(:),N2(:)); ylim([0 m]); xlim([0 m]); daspect([1 1 1]);
            hold on; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected neuron count (FlyEM)'); ylabel('connected neuron count (FlyWire)');
            title([roitype ' similarity=' num2str(NSims(i))]);

            % scatter plot of connected post-synapse count
            m = max([S1(:); S2(:)]);
            figure; scatter(S1(:),S2(:)); ylim([0 m]); xlim([0 m]); daspect([1 1 1]);
            hold on; plot([0 m], [0 m],':','Color',[0.5 0.5 0.5]); hold off; xlabel('connected synapse count (FlyEM)'); ylabel('connected synapse count (FlyWire)');
            title([roitype ' similarity=' num2str(SSims(i))]);
        end
    end
end

function checkSmoothingFlyWireByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {{'hemiCmkm',''},{'hemiCmkm','_fw'},{'hemiBranson7065km',''},{'hemiBranson7065km','_fw'},{'hemiDistKm',''},{'hemiDistKm','_fw'}};
    roitypelabels = {'Cm','CmFw','Branson','BransonFw','Dist','DistFw'};

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
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    roitypes = {'flyemroi','flyemroi_fw','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
    roitypelabels = {'FlyEM','FlyEmFw','Branson','Cm','CmR1w1','Dist','Rand','Vrand',};

    ylabels = {}; R3 = []; A3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
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

function checkNuisanceByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', '_', '_', '_'}; %40, three dummies
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {{'hemiBranson7065km',''},{'hemiCmkm',''},{'hemiCmkm','r1w1'},{'hemiDistKm',''},{'hemiRand',''},{'hemiVrand',''}};
    roitypelabels = {'Branson','Cm','CmR1w1','Dist','Rand','Vand'};

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
                        xlabels{ii} = [nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
 %   figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24 48 72 96 120]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

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

function checkLargeSmoothingPoltcompByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80', 's90', 's100', ...
        's110', 's120', 's130', 's140', 's150', 's160', 's170', 's180', 's190', 's200', ...
        's210', 's220', 's230', 's240', 's250', 's260', 's270', 's280', 's290', 's300'};
    nuisance = {'','poltcomp'};
    roinums = [50 100 500];
    roitypes = {{'hemiCmkm',''},{'hemiDistKm',''}};
    roitypelabels = {'Cm','Dist'};

    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = []; ii=1; xlabels = {};
        for rr=1:length(roinums)
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                for n=1:length(nuisance)
                    for k=1:length(smooth)
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(roinums(rr)) roitypes{r}{2}];
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
        end

        R3 = [R3;Rm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA);
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24]);
    figure; plot(B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});

    % find best smooth size
    i = 2; step = length(smooth); BMidx = [];
    for rr=1:length(roinums)
        Bt = B(:,(i-1)*step+1:i*step); i=i+2;
        [m,idx] = max(Bt,[],2); BMidx = [BMidx, idx];
        disp([roitypelabels{1} '   neuron ' nuisance{2} 'roi' num2str(roinums(rr)) ' max : ' smooth{idx(7)} '=' num2str(m(7))]);
        disp([roitypelabels{2} ' neuron ' nuisance{2} 'roi' num2str(roinums(rr)) ' max : ' smooth{idx(7+24)} '=' num2str(m(7+24))]);
    end

    % show thresholded AUCs
    X = []; slabels = {}; % neuron only
    for i=[7 7+24]
        for rr=1:length(roinums)
            Y = squeeze(AA3(i,:,BMidx(i,rr))); Y(Y==0)=nan; % ignore zero result
            X = [X, Y'];
            slabels = [slabels(:); ['roi' num2str(roinums(rr)) ' ' smooth{BMidx(i,rr)}]];
        end
    end
    figure; plot(X); legend(slabels); title(['FC-SC detection results by threshold in ' ylabels{i}]);
end

function checkNeuronVsSynapseByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {''};
    roinums = [100 500 1000 5000 10000 15000 20000];
    roitypes = {{'hemiCmkm',''},{'hemiCmkm','r1w1'},{'hemiDistKm',''},{'hemiVrand',''}};
    roitypelabels = {'Cm','CmR1w1','Dist','Vand'};
    
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
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24 48 72]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48 72]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48 72]);
    figure; plot(B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSmoothingNuisanceByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's30', 's80'};
    nuisance = {'', 'poltcomp'};
    roinums = [100 500 1000 5000 10000];
    roitypes = {{'hemiCmkm',''},{'hemiDistKm',''}};
    roitypelabels = {'Cm','Dist'};
    
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
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
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
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24]);
    figure; plot(B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSCdiffFlyWire(vslabels)
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
            info = niftiinfo(['atlas\' axlabels{ii} 'atlasCal.nii.gz']);
            V = niftiread(info);
            Vs(ii) = length(find(V>0));

            % check connection matrix
            xlabels{j} = ['roi' num2str(roinums(j))];
            roitype = [roitypes{i} num2str(roinums(j))];
            confile = ['data/' lower(roitype) '_connectlist.mat'];
            confilefw = ['data/' lower(roitype) '_fw_connectlist.mat'];
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
    set(gca,'XTick',1:size(mat,1));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end

