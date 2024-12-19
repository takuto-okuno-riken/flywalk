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
    checkSmoothingResult50(vslabels);

    % check correlation result in each ROI num (s0 to s80, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
    checkSmoothingByRoinum(vslabels);

    % check nuisance result round 50 ROIs (all nuisance)
    % roitype: FlyEM,FlyEmFw,Branson,Cm,CmR1w1,Dist,Rand,Vrand
    checkNuisanceResult50(vslabels);

    % check correlation result in each ROI num (all nuisance, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
    checkNuisanceByRoinum(vslabels);

    % check correlation result of large smoothing size (s0 to 300, roi 50 to 500, '' & poltcomp)
    % roitype: Cm,Dist
    checkLargeSmoothingPoltcompByRoinum(vslabels);

    % check correlation result in each ROI num (roi 100 to 20000)
    % roitype: Cm,CmR1w1,Dist,Vand
    checkNeuronVsSynapseByRoinum(vslabels);

    % check correlation result in each ROI num (s30,80,150,230,300, roi 50 to 10000, '' & poltcomp)
    % because LargeSmoothing showed better result around s230, this function extended range.
    % roitype: Cm,DistKm
    checkSmoothingNuisanceByRoinum(vslabels);

    % check smoothing result of FlyEM vs. FlyWire around 50 ROIs (s0 to s80)
    % roitype: FlyEM,FlyEMFw,DistKm50,DistKm50Fw,DistKm50Avg
    checkSmoothNuisanceFlyWireResult50(vslabels);

    % check correlation result of FlyEM vs. FlyWire (s0 to s80, roi 20 to 1000)
    % roitype: Cm,CmFw,Branson,BransonFw,Dist,DistFw
    checkSmoothingFlyWireByRoinum(vslabels);

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

        R3 = [R3;Rm]; A3 = [A3;Am]; ii=ii-1;
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
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48 72 96 120]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSmoothingResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
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
    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiDistKm50','hemiDistKm50_fw0sr50','hemiDistKm50_avg'};
    roitypelabels = {'FlyEM','FlyEMFw','DistKm50','DistKm50Fw','DistKm50Avg'};

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
    I = getR3idx([7 9 19 21],[0 24 48 72 96]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),smooth,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96]);
    figure; plot([0:8]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);
    
    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),smooth,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0 24 48 72 96]);
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
    I = getR3idx([7 9 19 21],[0 24 48 72 96]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),smooth,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96]);
    figure; plot([0:8]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);
end

function checkSmoothingFlyWireByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roinums = [20 30 50 100 200 300 500 1000];
    roitypes = {{'hemiCmkm',''},{'hemiCmkm','_fw0sr50'},{'hemiBranson7065km',''},{'hemiBranson7065km','_fw0sr50'},{'hemiDistKm',''},{'hemiDistKm','_fw0sr50'}};
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
    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
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

        R3 = [R3;Rm]; A3 = [A3;Am]; ii=ii-1;
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
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48 72 96 120]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72 96 120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48 72 96 120]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
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

        R3 = [R3;Rm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA); ii=ii-1;
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
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});

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
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48 72]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24 48 72]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48 72]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSmoothingNuisanceByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's30', 's80', 's150', 's230', 's300'};
    nuisance = {'', 'poltcomp'};
    roinums = [50 100 500 1000 5000 10000];
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
        R3 = [R3;Rm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA); ii=ii-1;
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
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    for i=[7 9]
        X = []; slabels = {};
        for j=1:length(roitypelabels)
            idx = [];
            for k=1:length(roinums)
                [m,kk] = max(A3(i+(j-1)*24,(k-1)*12+1:k*12),[],2);
                idx = (k-1)*12+kk;
                slabels = [slabels(:); [roitypelabels{j} ' ' xlabels{idx}]];
                Y = squeeze(AA3(i+(j-1)*24,:,idx)); Y(Y==0)=nan; % ignore zero result
                X = [X, Y'];
            end
        end
        figure; plot(X); legend(slabels); title(['FC-SC detection results by threshold in (' num2str(i) ') ' vslabels{i}]); setlineColors(6);
    end

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0 24]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
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

