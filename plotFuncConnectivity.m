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
    % FlyEM,Cm,Dist matrix are used in figure.1
    % roitype: FlyEM,FlyEMFw,Branson,BransonFw,Cm,CmFw,CmR1w1,Dist,Rand,Vrand
%    checkSmoothingResult50(vslabels);

    % check smoothing result in several ROI nums (s0 to s80, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
%    checkSmoothingByRoinum(vslabels);

    % check nuisance result around 50 ROIs (all nuisance)
    % FlyEM,Cm,Dist are used in figure.1
    % roitype: FlyEM,FlyEmFw,Branson,Cm,CmR1w1,Dist,Rand,Vrand
%    checkNuisanceResult50(vslabels);

    % check nuisance result in several ROI nums (all nuisance, roi 20 to 1000)
    % roitype: Branson,Cm,CmR1w1,Dist,Rand,Vand
%    checkNuisanceByRoinum(vslabels);

    % check nuisance result in each hemi ROI (all nuisance)
    % to check inside neuropil FC-SC relation at 1 voxel resolution
    % roitype: hemiRoi1 to 113
%    checkNuisanceResultHemiROIs(vslabels);

    % check nuisance and smoothing in mushroom body (s30,80,150, '',6hm,tcomp,pol,poltcomp)
    % to check inside neuropil FC-SC relation at 1 voxel resolution based
    % on checkNuisanceResultHemiROIs result (poltcomp may not be the best).
    % roitype: hemiroi68-59-87-106-50-27-54 (mushroom body)
%    checkSmoothingNuisanceMushroomBody(vslabels);

    % check large smoothing size in several ROI nums and poltcomp (s0 to 300, roi 20 to 1000, '' & poltcomp)
    % Cm,Dist are used in figure.2
    % roitype: Cm,CmR1w1,Dist,
    checkLargeSmoothingPoltcompByRoinum(vslabels);

    % check high-pass filter in several ROI nums and poltcomp (s230, roi 20 to 1000, poltcomp)
    % Dist are used in ext figure.2
    % roitype: Dist,
%    checkHighpassFilterPoltcompByRoinum(vslabels);

    % check large ROI num result (roi 100 to 20000)
    % roitype: Cm,CmR1w1,Dist,Vand
%    checkLargeRoinumResult(vslabels);

    % check large smoothing size in several ROI nums and poltcomp (s30,80,150,230,300, roi 50 to 20000, '' & poltcomp)
    % because LargeSmoothing showed better result around s230, this function extended range.
    % Cm,Dist are used in figure.2    
    % roitype: Cm,DistKm,Cube4
%    checkSmoothingPoltcompByLargeRoinum(vslabels);

    % check DistKm1000 in each voxel num and poltcomp (s30,80,150,230,300, roi 1000, '' & poltcomp)
    % to check how smoothing and voxel size related, and inside or inter neuropil relation (1 voxel resolution).
    % In 1 voxel resolution, DistKm1000vox1 (very sparce) with s230poltcomp shows still
    % high correlation & detection result. Thus, high density neuropil needs different solution.
    % roitype: DistKm1000vox128,64,32,16,8,4,2,1
%    checkSmoothingNuisanceByDistKm1000vox(vslabels);

    % check smoothing result of FlyEM vs. FlyWire around 50 ROIs (s0 to s80, no nuisance)
    % roitype: FlyEM,FlyEMFw,DistKm50,DistKm50Fw,DistKm50Avg
%    checkSmoothingFlyWireResult50(vslabels);

    % check smoothing and several ROI nums of FlyEM vs. FlyWire (s0 to s80, roi 20 to 1000, no nuisance)
    % roitype: Cm,CmFw,Branson,BransonFw,Dist,DistFw
%    checkSmoothingFlyWireByRoinum(vslabels);

    % check reciprocal synapse distance thresholds of FlyEM vs. FlyWire (s0, poltcomp)
    % roitype: hemiroi
%    checkReciprocalDistanceRandByHemiroi(vslabels);

    % check synapse separation index thresholds of FlyEM vs. FlyWire (s0, poltcomp)
    % roitype: hemiroi
%    checkSeparationRandByHemiroi(vslabels);

%    checkRandSubsampleByHemiroi(vslabels); % find good rand param of permutation test of separation index & reciprocal

    checkRandSubsampleRankTestByHemiroi(vslabels); % show permutation test results (figure.4)

    % check reciprocal synapse distance thresholds of FlyEM vs. FlyWire (s230, poltcomp)
    % roitype: DistKm500
%    checkReciprocalDistanceRandByDistKm(vslabels);

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
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0:24:120]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0:24:120]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0:24:120]);
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
    I = getR3idx([7 9 19 21],[0:24:168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),smooth,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0:24:120]);% 144 168 192 216]);
    figure; plot([0:8]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);
    
    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0:24:168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),smooth,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    I = getR3idx([7 9],[0:24:120]);% 144 168 192 216]);
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
    I = getR3idx([7 9 19 21],[0:24:168]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),smooth,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    I = getR3idx([7 9],[0:24:120]);% 144 168 192 216]);
    figure; plot([0:8]',B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(2);
end

function checkSmoothingFlyWireResult50(vslabels)
    % around 50 clusters
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    nuisance = {''};
    roitypes = {'hemiroi','hemiroi_fw0sr140','hemiDistKm50','hemiDistKm50_fw0sr140','hemiDistKm50_avg'};
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
    roitypes = {{'hemiCmkm',''},{'hemiCmkm','_fw0sr140'},{'hemiBranson7065km',''},{'hemiBranson7065km','_fw0sr140'},{'hemiDistKm',''},{'hemiDistKm','_fw0sr140'}};
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
    K = [0:24:120];
    I = getR3idx([7 9 19 21],K);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);

    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],K);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],K);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9],K);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],K);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],K);
    figure; plot(B(I,:)'); legend(ylabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkReciprocalDistanceRandByHemiroi(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'poltcomp'};
    % !caution! hemiroi has 63 ROIs. but hemiroi_hb0sr80 has 52 ROIs. 52 ROIs should be used.
    roitypes = {{'hemiroi','_hb0sr80'}, ...
        {'hemiroi','_hb0sr80_rc20_only1'},{'hemiroi','_hb0sr80_rc40_only1'},{'hemiroi','_hb0sr80_rc100_only1'},{'hemiroi','_hb0sr80_rc500_only1'},{'hemiroi','_hb0sr80_rc1000_only1'},{'hemiroi','_hb0sr80_rc10000_only1'}, ...
        {'hemiroi','_fw0sr140'}, ...
        {'hemiroi','_fw0sr140_rc20_only1'},{'hemiroi','_fw0sr140_rc40_only1'},{'hemiroi','_fw0sr140_rc100_only1'},{'hemiroi','_fw0sr140_rc500_only1'},{'hemiroi','_fw0sr140_rc1000_only1'},{'hemiroi','_fw0sr140_rc10000_only1'}, ...
%        {'hemiroi','_hb0sr80_rc10000'},{'hemiroi','_hb0sr80_rc10000_rand1'},{'hemiroi','_hb0sr80_rc10000_rand2'},{'hemiroi','_hb0sr80_rc10000_rand3'}, ...
%        {'hemiroi','_hb0sr80_rc10000_xrand1'},{'hemiroi','_hb0sr80_rc10000_xrand2'},{'hemiroi','_hb0sr80_rc10000_xrand3'}, ...
%        {'hemiroi','_hb0sr80fw_rc20_xorand1'},{'hemiroi','_hb0sr80fw_rc20_xorand2'},{'hemiroi','_hb0sr80fw_rc20_xorand3'}, ...
%        {'hemiroi','_hb0sr80_rc20_xorand1'},{'hemiroi','_hb0sr80_rc20_xorand2'},{'hemiroi','_hb0sr80_rc20_xorand3'}, ...
%        {'hemiroi','_hb0sr80_rc10000_xorand1'},{'hemiroi','_hb0sr80_rc10000_xorand2'},{'hemiroi','_hb0sr80_rc10000_xorand3'}, ...
%        {'hemiroi','_fw0sr140_rc10000'},{'hemiroi','_fw0sr140_rc10000_rand1'},{'hemiroi','_fw0sr140_rc10000_rand2'},{'hemiroi','_fw0sr140_rc10000_rand3'}, ...
%        {'hemiroi','_fw0sr140_rc10000_xrand1'},{'hemiroi','_fw0sr140_rc10000_xrand2'},{'hemiroi','_fw0sr140_rc10000_xrand3'}, ...
%        {'hemiroi','_fw0sr140_rc20_xorand1'},{'hemiroi','_fw0sr140_rc20_xorand2'},{'hemiroi','_fw0sr140_rc20_xorand3'}, ...
%        {'hemiroi','_fw0sr140_rc10000_xorand1'},{'hemiroi','_fw0sr140_rc10000_xorand2'},{'hemiroi','_fw0sr140_rc10000_xorand3'}, ...
        };
    roitypelabels = {'HbRc0','o20','o40','o100','o500','o1000','o10000','FwRc0','',''};

    R3 = []; A3 = []; ii=1; rlabels = {}; ncounts = []; sycounts = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(A3,1),100);
                        R = nan(size(R3,1),1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end
        rlabels{ii} = [roitypes{r}{1} roitypes{r}{2} ]; ii=ii+1;
        conmat = ['results/sc/' roitypes{r}{1} roitypes{r}{2} '_connectlist.mat'];
        if exist(conmat,'file')
            load(conmat);
            nn = sum(ncountMat(primaryIds,primaryIds,2),'all');
            syn = sum(sycountMat(primaryIds,primaryIds,2),'all');
            ncounts = [ncounts; nn];
            sycounts = [sycounts; syn];
        else
            ncounts = [ncounts; NaN];
            sycounts = [sycounts; NaN];
        end

        R3 = [R3,Rm]; A3 = [A3,Am];
    end

    % show synapse count in each ROI type.
    cats=categorical(rlabels,rlabels);
    figure; bar(cats,ncounts); title('neuron count matrix total')
    figure; bar(cats,sycounts); title('synapse count matrix total')

    % FC-SC correlation (all)
    I = [7 9 19 21];
    figure; imagescLabel2(R3(I,:),roitypelabels,vslabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    figure; plot(R3(I,:)'); legend(vslabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    figure; imagescLabel2(A3(I,:),roitypelabels,vslabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    figure; plot(A3(I,:)'); legend(vslabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    figure; imagescLabel2(B(I,:),roitypelabels,vslabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    figure; plot(B(I,:)'); legend(vslabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkRandSubsampleByHemiroi(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'poltcomp'};
    roitypes = {{'hemiroi','_hb0sr80'}, ...
        {'hemiroi','_hb0sr80_rn50_orand1'},{'hemiroi','_hb0sr80_rn150_orand1'},{'hemiroi','_hb0sr80_rn500_orand1'},{'hemiroi','_hb0sr80_rn1000_orand1'}, ...
        {'hemiroi','_hb0sr80fw_rn50_orand1'},{'hemiroi','_hb0sr80fw_rn60_orand1'},{'hemiroi','_hb0sr80fw_rn70_orand1'}, ...
        {'hemiroi','_hb0sr80fw_rn130_orand1'},{'hemiroi','_hb0sr80fw_rn150_orand1'}, ...
        {'hemiroi','_hb0sr80fw_rn500_orand1'},{'hemiroi','_hb0sr80fw_rn600_orand1'},{'hemiroi','_hb0sr80fw_rn700_orand1'},{'hemiroi','_hb0sr80fw_rn710_orand1'},{'hemiroi','_hb0sr80fw_rn1000_orand1'}, ...
        {'hemiroi','_fw0sr140'}, ...
        {'hemiroi','_fw0sr140_rn50_orand1'},{'hemiroi','_fw0sr140_rn120_orand1'},{'hemiroi','_fw0sr140_rn130_orand1'},{'hemiroi','_fw0sr140_rn140_orand1'},{'hemiroi','_fw0sr140_rn150_orand1'}, ...
        {'hemiroi','_fw0sr140_rn200_orand1'},{'hemiroi','_fw0sr140_rn250_orand1'},{'hemiroi','_fw0sr140_rn500_orand1'},{'hemiroi','_fw0sr140_rn1000_orand1'},{'hemiroi','_fw0sr140_rn1530_orand1'},{'hemiroi','_fw0sr140_rn1550_orand1'},{'hemiroi','_fw0sr140_rn2000_orand1'}, ...
        };
    roitypelabels = {'HbRc0','r50','r150','r500','r1000','r50','r150','r500','r1000','FwRc0','FwRc10000','r50','r150','r500','r1000'};

    R3 = []; A3 = []; ii=1; rlabels = {}; ncounts = []; sycounts = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(A3,1),100);
                        R = nan(size(R3,1),1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end
        rlabels{ii} = [roitypes{r}{1} roitypes{r}{2} ]; ii=ii+1;
        conmat = ['results/sc/' roitypes{r}{1} roitypes{r}{2} '_connectlist.mat'];
        if exist(conmat,'file')
            load(conmat);
            nn = sum(ncountMat(primaryIds,primaryIds,2),'all');
            syn = sum(sycountMat(primaryIds,primaryIds,2),'all');
            ncounts = [ncounts; nn];
            sycounts = [sycounts; syn];
        else
            ncounts = [ncounts; NaN];
            sycounts = [sycounts; NaN];
        end

        R3 = [R3,Rm]; A3 = [A3,Am];
    end

    % show synapse count in each ROI type.
    cats=categorical(rlabels,rlabels);
    figure; bar(cats,ncounts); title('neuron count matrix total')
    figure; bar(cats,sycounts); title('synapse count matrix total')

    % FC-SC correlation (all)
    I = [7 9 19 21];
    figure; imagescLabel2(R3(I,:),roitypelabels,vslabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    figure; plot(R3(I,:)'); legend(vslabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    figure; imagescLabel2(A3(I,:),roitypelabels,vslabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    figure; plot(A3(I,:)'); legend(vslabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    figure; imagescLabel2(B(I,:),roitypelabels,vslabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    figure; plot(B(I,:)'); legend(vslabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkRandSubsampleRankTestByHemiroi(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'poltcomp'};
%    roitypes = {{'hemiroi','_hb0sr80fw_rd140-15'},{'hemiroi','_hb0sr80fw_rd705-40'},{'hemiroi','_hb0sr80fw_rd67-12'}, ...
%       {'hemiroi','_fw0sr140_rd185-20'},{'hemiroi','_fw0sr140_rd1560-80'},{'hemiroi','_fw0sr140_rd140-15'},};
    roitypes = {{'hemiroi','_fw0sr140_rd185-20'},{'hemiroi','_fw0sr140_rd1560-80'},{'hemiroi','_fw0sr140_rd140-15'},};
%    targets = {{'hemiroi','_hb0sr80_sp10db3000mi1_only1'},{'hemiroi','_hb0sr80_sp90db3000mi1'},{'hemiroi','_hb0sr80_rc20_only1'}, ...
%        {'hemiroi','_fw0sr140_sp10db3000mi1_only1'},{'hemiroi','_fw0sr140_sp90db3000mi1'},{'hemiroi','_fw0sr140_rc20_only1'},};
    targets = {{'hemiroi','_fw0sr140_sp10db3000mi1_only1'},{'hemiroi','_fw0sr140_sp90db3000mi1'},{'hemiroi','_fw0sr140_rc20_only1'},};
    rNums = 1:99; % for random subsampling number
    pval = 0.05 / 6; % bonferroni correction

    for r = 1:length(roitypes)
        Am = []; Rm = []; tAm = []; tRm = []; ncounts = []; sycounts = [];
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    % target ROI type
                    pftype = [smooth{k} hpfstr nuisance{n} preproc targets{r}{1} targets{r}{2}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(Am,1),100);
                        R = nan(size(Rm,1),1);
                    end
                    tAm = [tAm,nanmean(A,2)];
                    tRm = [tRm,R(:)];

                    % random subsampling ROI type
                    for ss = 1:length(rNums)
                        rstr = '';
                        if rNums(ss) > 0, rstr = ['-' num2str(rNums(ss))]; end
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2} rstr];
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
            % target ROI type
            conmat = ['results/sc/' targets{r}{1} targets{r}{2} '_connectlist.mat'];
            load(conmat);
            tncount = sum(ncountMat(primaryIds,primaryIds,2),'all');
            tsycount = sum(sycountMat(primaryIds,primaryIds,2),'all');
            % random subsampling ROI type
            for ss = 1:length(rNums)
                rstr = '';
                if rNums(ss) > 0, rstr = ['-' num2str(rNums(ss))]; end
                conmat = ['results/sc/' roitypes{r}{1} roitypes{r}{2} rstr '_connectlist.mat'];
                if exist(conmat,'file')
                    load(conmat);
                    nn = sum(ncountMat(primaryIds,primaryIds,2),'all');
                    syn = sum(sycountMat(primaryIds,primaryIds,2),'all');
                    ncounts = [ncounts; nn];
                    sycounts = [sycounts; syn];
                end
            end
        end
        if isempty(ncounts), continue; end

        % normality test
        [h,p] = lillietest(ncounts);
        disp(['neuron count, Lilliefors test p=' num2str(p)]);
        [h,p] = lillietest(sycounts);
        disp(['synapse count, Lilliefors test p=' num2str(p)]);
        pd1 = fitdist(ncounts,'Normal');
        p1 = 2*normcdf(-abs(tncount-pd1.mu)/pd1.sigma); % two tailed normal distribution based test
        pd2 = fitdist(sycounts,'Normal');
        p2 = 2*normcdf(-abs(tsycount-pd2.mu)/pd2.sigma); % two tailed normal distribution based test
        
        % show total synapse count histogram
        figure; h=histogram(ncounts, 10); title(['ncounts count matrix total : ' targets{r}{2} ' p=' num2str(p1)]);
        ymax = 5*ceil(max(h.Values)/5);
        [y, x]=ecdf([tncount; ncounts]); hold on; plot(x, y*ymax,'LineWidth',2); hold off;
        x=linspace(min(x),max(x)); y=normcdf(x,pd1.mu,pd1.sigma); hold on; plot(x, y*ymax,':k','LineWidth',0.5); hold off;
        xline(tncount,'r'); bx=norminv([pval,1-pval],pd1.mu,pd1.sigma); xline(bx,':r');

        figure; h=histogram(sycounts, 10); title(['synapse count matrix total : ' targets{r}{2} ' p=' num2str(p2)]);
        ymax = 5*ceil(max(h.Values)/5);
        [y, x]=ecdf([tsycount; sycounts]); hold on; plot(x, y*ymax,'LineWidth',2); hold off;
        x=linspace(min(x),max(x)); y=normcdf(x,pd2.mu,pd2.sigma); hold on; plot(x, y*ymax,':k','LineWidth',0.5); hold off;
        xline(tsycount,'r'); bx=norminv([pval,1-pval],pd2.mu,pd2.sigma); xline(bx,':r');
%        continue;

        % show R and AUC histogram
        B = abs(Rm) + abs(Am-0.5)*2;
        tB = abs(tRm) + abs(tAm-0.5)*2;
        I = [7 9];
        for i=1:length(I)
            ii = I(i);
            % normality test
            [h,p] = lillietest(Rm(ii,:));
            disp([targets{r}{2} ' corr, Lilliefors test p=' num2str(p)]);
            [h,p] = lillietest(Am(ii,:));
            disp([targets{r}{2} ' AUC, Lilliefors test p=' num2str(p)]);
            [h,p] = lillietest(B(ii,:));
            disp([targets{r}{2} ' all, Lilliefors test p=' num2str(p)]);
            pd1 = fitdist(Rm(ii,:)','Normal');
            p1 = 2*normcdf(-abs(tRm(ii)-pd1.mu)/pd1.sigma); % two tailed normal distribution based test
            pd2 = fitdist(Am(ii,:)','Normal');
            p2 = 2*normcdf(-abs(tAm(ii)-pd2.mu)/pd2.sigma); % two tailed normal distribution based test
            pd3 = fitdist(B(ii,:)','Normal');
            p3 = 2*normcdf(-abs(tB(ii)-pd3.mu)/pd3.sigma);  % two tailed normal distribution based test
            
            % show histogram
            figure; h=histogram(Rm(ii,:), 10); title(['FC-SC corr : ' targets{r}{2} ' : ' vslabels{ii} ' p=' num2str(p1)]);
            ymax = 5*ceil(max(h.Values)/5);
            [y, x]=ecdf([tRm(ii), Rm(ii,:)]); hold on; plot(x, y*ymax,'LineWidth',2); hold off;
            x=linspace(min(x),max(x)); y=normcdf(x,pd1.mu,pd1.sigma); hold on; plot(x, y*ymax,':k','LineWidth',0.5); hold off;
            xline(tRm(ii),'r'); bx=norminv([pval,1-pval],pd1.mu,pd1.sigma); xline(bx,':r');

            figure; h=histogram(Am(ii,:), 10); title(['FC-SC AUC : ' targets{r}{2} ' : ' vslabels{ii} ' p=' num2str(p2)]);
            ymax = 5*ceil(max(h.Values)/5);
            [y, x]=ecdf([tAm(ii), Am(ii,:)]); hold on; plot(x, y*ymax,'LineWidth',2); hold off;
            x=linspace(min(x),max(x)); y=normcdf(x,pd2.mu,pd2.sigma); hold on; plot(x, y*ymax,':k','LineWidth',0.5); hold off;
            xline(tAm(ii),'r'); bx=norminv([pval,1-pval],pd2.mu,pd2.sigma); xline(bx,':r');
%            figure; histogram(B(ii,:), 10); title(['FC-SC all : ' targets{r}{2} ' : ' vslabels{ii} ' p=' num2str(p3)]);
%            [y, x]=ecdf([tB(ii), B(ii,:)]); hold on; plot(x, y*25,'LineWidth',2); hold off;
%            x=linspace(min(x),max(x)); y=normcdf(x,pd3.mu,pd3.sigma); hold on; plot(x, y*25,':k','LineWidth',0.5); hold off;
%            xline(tB(ii),'--r');
        end
    end
end

function checkReciprocalDistanceRandByDistKm(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};
    roitypes = {
        {'hemidistkm500',''}, ...
        {'hemidistkm500','_hb0sr80_rc20_only1'},{'hemidistkm500','_hb0sr80_rc40_only1'},{'hemidistkm500','_hb0sr80_rc100_only1'},{'hemidistkm500','_hb0sr80_rc500_only1'},{'hemidistkm500','_hb0sr80_rc1000_only1'},{'hemidistkm500','_hb0sr80_rc10000_only1'}, ...
        {'hemidistkm500','_hb0sr80_rc20_xorand1'},{'hemidistkm500','_hb0sr80_rc20_xorand2'},{'hemidistkm500','_hb0sr80_rc20_xorand3'}, ...
        {'hemidistkm500','_fw0sr140'}, ...
        {'hemidistkm500','_fw0sr140_rc20_only1'},{'hemidistkm500','_fw0sr140_rc40_only1'},{'hemidistkm500','_fw0sr140_rc100_only1'},{'hemidistkm500','_fw0sr140_rc500_only1'},{'hemidistkm500','_fw0sr140_rc1000_only1'},{'hemidistkm500','_fw0sr140_rc10000_only1'}, ...
        {'hemidistkm500','_fw0sr140_rc20_xorand1'},{'hemidistkm500','_fw0sr140_rc20_xorand2'},{'hemidistkm500','_fw0sr140_rc20_xorand3'}, ...
        };
    roitypelabels = {'HbRc0','o20','o40','o100','o500','o1000','o10000','xor1','xor2','xor3','FwRc0','o20','o40','o100','o500','o1000','o10000','xor1','xor2','xor3'};

    R3 = []; A3 = []; ii=1; rlabels = {}; ncounts = []; sycounts = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(A3,1),100);
                        R = nan(size(R3,1),1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end
        rlabels{ii} = [roitypes{r}{1} roitypes{r}{2} ]; ii=ii+1;
        conmat = ['results/sc/' roitypes{r}{1} roitypes{r}{2} '_connectlist.mat'];
        if exist(conmat,'file')
            load(conmat);
            nn = sum(ncountMat(primaryIds,primaryIds,2),'all');
            syn = sum(sycountMat(primaryIds,primaryIds,2),'all');
            ncounts = [ncounts; nn];
            sycounts = [sycounts; syn];
        else
            ncounts = [ncounts; NaN];
            sycounts = [sycounts; NaN];
        end

        R3 = [R3,Rm]; A3 = [A3,Am];
    end

    % show synapse count in each ROI type.
    cats=categorical(rlabels,rlabels);
    figure; bar(cats,ncounts); title('neuron count matrix total')
    figure; bar(cats,sycounts); title('synapse count matrix total')

    % FC-SC correlation (all)
    I = [7 9 19 21];
    figure; imagescLabel2(R3(I,:),roitypelabels,vslabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    figure; plot(R3(I,:)'); legend(vslabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    figure; imagescLabel2(A3(I,:),roitypelabels,vslabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    figure; plot(A3(I,:)'); legend(vslabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    figure; imagescLabel2(B(I,:),roitypelabels,vslabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    figure; plot(B(I,:)'); legend(vslabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSeparationRandByHemiroi(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'poltcomp'};
    % !caution! hemiroi has 63 ROIs. but hemiroi_hb0sr80 has 52 ROIs. 52 ROIs should be used.
    roitypes = {
        {'hemiroi','_hb0sr80'}, ...
        {'hemiroi','_hb0sr80_sp10db3000mi1'},{'hemiroi','_hb0sr80_sp20db3000mi1'},{'hemiroi','_hb0sr80_sp40db3000mi1'},{'hemiroi','_hb0sr80_sp60db3000mi1'},{'hemiroi','_hb0sr80_sp80db3000mi1'},{'hemiroi','_hb0sr80_sp90db3000mi1'}, ...
        {'hemiroi','_hb0sr80_sp10db3000mi1_only1'},{'hemiroi','_hb0sr80_sp20db3000mi1_only1'},{'hemiroi','_hb0sr80_sp40db3000mi1_only1'},{'hemiroi','_hb0sr80_sp60db3000mi1_only1'},{'hemiroi','_hb0sr80_sp80db3000mi1_only1'},{'hemiroi','_hb0sr80_sp90db3000mi1_only1'}, ...
        {'hemiroi','_fw0sr140'}, ...
        {'hemiroi','_fw0sr140_sp10db3000mi1'},{'hemiroi','_fw0sr140_sp20db3000mi1'},{'hemiroi','_fw0sr140_sp40db3000mi1'},{'hemiroi','_fw0sr140_sp60db3000mi1'},{'hemiroi','_fw0sr140_sp80db3000mi1'},{'hemiroi','_fw0sr140_sp90db3000mi1'}, ...
        {'hemiroi','_fw0sr140_sp10db3000mi1_only1'},{'hemiroi','_fw0sr140_sp20db3000mi1_only1'},{'hemiroi','_fw0sr140_sp40db3000mi1_only1'},{'hemiroi','_fw0sr140_sp60db3000mi1_only1'},{'hemiroi','_fw0sr140_sp80db3000mi1_only1'},{'hemiroi','_fw0sr140_sp90db3000mi1_only1'}, ...
        };
    roitypelabels = {'Hb','10','20','40','60','80','90','o10','o20','o40','o60','o80','o90','Fw','10','20','40','60','80','90'};

    R3 = []; A3 = []; ii=1; rlabels = {}; ncounts = []; sycounts = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; 
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2}];
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(size(A3,1),100);
                        R = nan(size(R3,1),1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end
        rlabels{ii} = [roitypes{r}{1} roitypes{r}{2} ]; ii=ii+1;
        conmat = ['results/sc/' roitypes{r}{1} roitypes{r}{2} '_connectlist.mat'];
        if exist(conmat,'file')
            load(conmat);
            nn = sum(ncountMat(primaryIds,primaryIds,2),'all');
            syn = sum(sycountMat(primaryIds,primaryIds,2),'all');
            ncounts = [ncounts; nn];
            sycounts = [sycounts; syn];
        else
            ncounts = [ncounts; NaN];
            sycounts = [sycounts; NaN];
        end

        R3 = [R3,Rm]; A3 = [A3,Am];
    end

    % show synapse count in each ROI type.
    cats=categorical(rlabels,rlabels);
    figure; bar(cats,ncounts); title('neuron count matrix total')
    figure; bar(cats,sycounts); title('synapse count matrix total')

    % FC-SC correlation (all)
    I = [7 9 19 21];
    figure; imagescLabel2(R3(I,:),roitypelabels,vslabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    figure; plot(R3(I,:)'); legend(vslabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    figure; imagescLabel2(A3(I,:),roitypelabels,vslabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    figure; plot(A3(I,:)'); legend(vslabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    figure; imagescLabel2(B(I,:),roitypelabels,vslabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    figure; plot(B(I,:)'); legend(vslabels(I)); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
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
%    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
%    roitypelabels = {'FlyEM','FlyEmFw','Branson','Cm','CmR1w1','Dist','Rand','Vrand',};
    roitypes = {'hemiroi','hemiCmkm50','hemiDistKm50'};
    roitypelabels = {'FlyEM','Cm','Dist'};

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
%    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    I = getR3idx([7 9],[0:24:48]);
    figure; imagescLabel2(R3(I,:),nuisance,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation around 50 ROIs']); colormap(hot);
%    figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation around 50 ROIs']); setlineColors(24);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
%    I = getR3idx([7 9],[0:24:144]);
    figure; plot([1:37]',R3(I,:)'); legend(ylabels(I)); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2);

    % FC-SC detection (all)
%    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),nuisance,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection around 50 ROIs']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection around 50 ROIs']); setlineColors(24);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
%    I = getR3idx([7 9],[0:24:144]);
    figure; plot([1:37]',A3(I,:)'); legend(ylabels(I)); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2);

    % both FC-SC correlation & detection
    B = abs(R3) + abs(A3-0.5)*2;
%    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),nuisance,ylabels(I),[0 1.6]); colorbar; title('FC-SC correlation & detection around 50 ROIs'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection around 50 ROIs'); setlineColors(24);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
%    I = getR3idx([7 9],[0:24:144]);
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
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0:24:120]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0:24:120]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0:24:120]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0:24:120]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkNuisanceResultHemiROIs(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
    roitypes = {'hemiRoi1','hemiRoi5','hemiRoi7','hemiRoi27','hemiRoi30','hemiRoi32','hemiRoi43','hemiRoi52', ...
        'hemiRoi54','hemiRoi57','hemiRoi59','hemiRoi63','hemiRoi65','hemiRoi67','hemiRoi78','hemiRoi82', ...
        'hemiRoi89','hemiRoi93','hemiRoi95','hemiRoi100','hemiRoi101','hemiRoi106','hemiRoi113'};
    roitypelabels = {'Roi1','Roi5','Roi7','Roi27','Roi30','Roi32','Roi43','Roi52', ...
        'Roi54','Roi57','Roi59','Roi63','Roi65','Roi67','Roi78','Roi82', ...
        'Roi89','Roi93','Roi95','Roi100','Roi101','Roi106','Roi113'};

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
                        A = nan(24,100);
                        R = nan(24,1);
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
    nI = getR3idx([7],[0:24:24*(length(roitypes)-1)]);  % show only Traced neuron
    sI = getR3idx([9],[0:24:24*(length(roitypes)-1)]);  % show only Traced neuron, synapse
%    figure; imagescLabel2(R3(nI,:),nuisance,ylabels(nI),[0 0.5]); colorbar; title(['FC-SC correlation neuron hemi ROIs']); colormap(hot);
%    figure; imagescLabel2(R3(sI,:),nuisance,ylabels(sI),[0 0.5]); colorbar; title(['FC-SC correlation synapse hemi ROIs']); colormap(hot);

    % FC-SC correlation m-FCz Traced neuron vs synapse (neuron count shows better result)
    figure; plot([1:37]',R3(nI,:)'); legend(ylabels(nI)); title('FC-SC correlation Traced neuron'); setlineColors(5);
    hold on; boxplot(R3(nI,:)); hold off;
    figure; plot([1:37]',R3(sI,:)'); legend(ylabels(sI)); title('FC-SC correlation Traced synapse'); setlineColors(5);
    hold on; boxplot(R3(sI,:)); hold off;

    % FC-SC detection (all)
%    figure; imagescLabel2(A3(nI,:),nuisance,ylabels(nI),[0.5 1]); colorbar; title(['FC-SC detection hemi ROIs']); colormap(hot);

    % FC-SC detection m-FCz neuron vs synapse (neuron count shows better result)
    figure; plot([1:37]',A3(nI,:)'); legend(ylabels(nI)); title('FC-SC detection Traced neuron'); setlineColors(5);
    hold on; boxplot(A3(nI,:)); hold off;
    figure; plot([1:37]',A3(sI,:)'); legend(ylabels(sI)); title('FC-SC detection Traced synapse'); setlineColors(5);
    hold on; boxplot(A3(sI,:)); hold off;

    % both FC-SC correlation & detection
    B = abs(R3) + abs(A3-0.5)*2;
%    figure; imagescLabel2(B(nI,:),nuisance,ylabels(nI),[0 1.5]); colorbar; title('FC-SC correlation & detection hemi ROIs'); colormap(hot);

    % FC-SC correlation & detection Full vs. Traced (which is best?)
    figure; plot([1:37]',B(nI,:)'); legend(ylabels(nI)); title('FC-SC correlation & detection neuron hemi ROIs'); setlineColors(5);
    hold on; boxplot(B(nI,:)); hold off;
    figure; plot([1:37]',B(sI,:)'); legend(ylabels(nI)); title('FC-SC correlation & detection synapse hemi ROIs'); setlineColors(5);
    hold on; boxplot(B(sI,:)); hold off;
end

function checkSmoothingNuisanceMushroomBody(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'','s30','s40','s60','s80','s100','s150'};
    nuisance = {'','6hm','tcomp','pol','poltcomp' };
    roitypes = {'hemiroi68-59-87-106-50-27-54'};
    roitypelabels = {'MB'};

    ylabels = {}; R3 = []; A3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; ii=1; xlabels = {};
        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
            for n=1:length(nuisance)
                for k=1:length(smooth)
                    pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}];
                    xlabels{ii} = [smooth{k} nuisance{n}]; ii=ii+1;
                    aucmat = ['results/auc/' pftype '-fcauc.mat'];
                    if exist(aucmat,'file')
                        load(aucmat);
                    else
                        A = nan(24,100);
                        R = nan(24,1);
                    end
                    Am = [Am,nanmean(A,2)];
                    Rm = [Rm,R(:)];
                end
            end
        end

        R3 = [R3;Rm]; A3 = [A3;Am]; ii=ii+1;
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9 19 21],[0]);
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9 19 21], [0]);
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9 19 21],[0]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkLargeSmoothingPoltcompByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80', 's90', 's100', ...
        's110', 's120', 's130', 's140', 's150', 's160', 's170', 's180', 's190', 's200', ...
        's210', 's220', 's230', 's240', 's250', 's260', 's270', 's280', 's290', 's300'};
%    nuisance = {'','poltcomp'};
    nuisance = {'poltcomp'};
    roinums = [20 50 100 200 500 1000];
%    roitypes = {{'hemiCmkm',''},{'hemiCmkm','r1w1'},{'hemiDistKm',''}};
%    roitypelabels = {'Cm','CmR1w1','Dist'};
    roitypes = {{'hemiCmkm',''},{'hemiDistKm',''}};
    roitypelabels = {'Cm','Dist'};

    ylabels = {}; R3 = []; SR3 = []; A3 = []; AA3 = []; R3r = []; SR3r = []; A3r = []; I=[7 9];
    for r = 1:length(roitypes)
        Am = []; Rm = []; SRm = []; AA = []; ii=1; xlabels = {};
        for rr=1:length(roinums)
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                for n=1:length(nuisance)
                    Am2 = []; Rm2 = []; SRm2 = [];
                    for k=1:length(smooth)
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(roinums(rr)) roitypes{r}{2}];
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
                        if exist(aucmat,'file')
                            load(aucmat);
                        else
                            A = nan(size(Am,1),100);
                            R = nan(1,size(Rm,1));
                            SR = nan(1,size(SRm,1));
                        end
                        AA = cat(3,AA,A);
                        Am = [Am,nanmean(A,2)];
                        Am2 = [Am2,nanmean(A(I,:),2)];
                        Rm = [Rm,R(:)];
                        Rm2 = [Rm2,R(I)'];
                        SRm = [SRm,SR(:)];
                        SRm2 = [SRm2,SR(I)'];
                    end
                    R3r = [R3r;Rm2]; SR3r = [SR3r;SRm2]; A3r = [A3r;Am2];
                end
            end
        end

        R3 = [R3;Rm]; SR3 = [SR3;SRm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA); ii=ii-1;
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    L = (length(roitypes)-1) * 24;
    I = getR3idx([7 9 19 21],[0:24:L]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0:24:L]);
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC correlation (spearman) Traced neuron vs synapse
    I = getR3idx([7 9],[0:24:L]);
    figure; plot(SR3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation (spearman) Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0:24:L]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0:24:L]);
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0:24:L]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0:24:L]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
%{
    % comment out. These lines need '' & poltcomp, r1w1.
    % used for poltcomp-roi20-1000_smooth0-300ns.xlsx (Figure.2) 
    x = [3 4 7 8 11 12 15 16 19 20 23 24 27 28];
    y = [1 2 5 6 9 10 13 14 17 18 21 22 25 26];
    X = [x, x+56, y, y+56]; % excel order (copy & paste R3r(:,X), A3r(:,X), Br(:,X))
    R3r = R3r'; figure; plot(R3r); title('FC-SC correlation by smoothing size'); xlabel('smoothing size [voxel]');
    SR3r = SR3r'; figure; plot(SR3r); title('FC-SC correlation (Spearman) by smoothing size'); xlabel('smoothing size [voxel]');
    A3r = A3r'; figure; plot(A3r); title('FC-SC detection by smoothing size'); xlabel('smoothing size [voxel]');
    Br = abs(R3r) + abs(A3r-0.5)*2;
    figure; boxplot(Br(:,[1 3 13 15 25 27 [1 3 13 15 25 27]+56])); title('FC-SC detection & correlation, raw vs. poltcomp');

    p(1) = ranksum(Br(:,1),Br(:,3));
    p(2) = ranksum(Br(:,13),Br(:,15));
    p(3) = ranksum(Br(:,25),Br(:,27));
    p(4) = ranksum(Br(:,1+56),Br(:,3+56));
    p(5) = ranksum(Br(:,13+56),Br(:,25+56));
    p(6) = ranksum(Br(:,25+56),Br(:,27+56));
    disp(['ranksums : p=' num2str(p)]);

    % find best smooth size
    i = 1; step = length(smooth); BMidx = []; N = length(nuisance);
    for rr=1:length(roinums)
        Bt = B(:,(i-1)*step+1:i*step); i=i+N;
        [m,idx] = max(Bt,[],2); BMidx = [BMidx, idx];
        disp([roitypelabels{1} '   neuron ' nuisance{N} 'roi' num2str(roinums(rr)) ' max : ' smooth{idx(7)} '=' num2str(m(7))]);
        disp([roitypelabels{2} ' neuron ' nuisance{N} 'roi' num2str(roinums(rr)) ' max : ' smooth{idx(7+24)} '=' num2str(m(7+24))]);
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
%}
end

function checkHighpassFilterPoltcompByRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0, 0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.001]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};
    roinums = [20 50 100 200 500 1000];
    roitypes = {{'hemiDistKm',''}};
    roitypelabels = {'Dist'};

    ylabels = {}; R3 = []; SR3 = []; A3 = []; AA3 = []; R3r = []; SR3r = []; A3r = []; I=[7 9];
    for r = 1:length(roitypes)
        Am = []; Rm = []; SRm = []; AA = []; ii=1; xlabels = {};
        for rr=1:length(roinums)
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    Am2 = []; Rm2 = []; SRm2 = [];
                    for h=1:length(hpfTh)
                        hpfstr = '';
                        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(roinums(rr)) roitypes{r}{2}];
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
                        if exist(aucmat,'file')
                            load(aucmat);
                        else
                            A = nan(size(Am,1),100);
                            R = nan(1,size(Rm,1));
                            SR = nan(1,size(SRm,1));
                        end
                        AA = cat(3,AA,A);
                        Am = [Am,nanmean(A,2)];
                        Am2 = [Am2,nanmean(A(I,:),2)];
                        Rm = [Rm,R(:)];
                        Rm2 = [Rm2,R(I)'];
                        SRm = [SRm,SR(:)];
                        SRm2 = [SRm2,SR(I)'];
                    end
                    R3r = [R3r;Rm2]; SR3r = [SR3r;SRm2]; A3r = [A3r;Am2];
                end
            end
        end

        R3 = [R3;Rm]; SR3 = [R3;SRm]; A3 = [A3;Am]; AA3 = cat(1,AA3,AA); ii=ii-1;
        C = cell(24,1); C(1:24) = {[roitypelabels{r} ' ']};
        ylabels = [ylabels(:); strcat(C(:),vslabels(:))];
    end

    % FC-SC correlation (all)
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3(I,:)'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3(I,:)'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % both FC-SC correlation & detection (all)
    B = abs(R3) + abs(A3-0.5)*2;
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});

    % used for poltcomp-roi20-1000_highpass.xlsx (ext Figure.2) 
    %(copy & paste R3r(:,:), A3r(:,:), Br(:,:))
    R3r = R3r'; figure; plot(R3r); title('FC-SC correlation by high-pass filter'); xlabel('high-pass filter');
    SR3r = SR3r'; figure; plot(SR3r); title('FC-SC correlation (Spearman) by high-pass filter'); xlabel('high-pass filter');
    A3r = A3r'; figure; plot(A3r); title('FC-SC detection by high-pass filter'); xlabel('high-pass filter');
    Br = abs(R3r) + abs(A3r-0.5)*2;
end

function checkLargeRoinumResult(vslabels)
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

function checkSmoothingPoltcompByLargeRoinum(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's30', 's80', 's150', 's230', 's300'};
    nuisance = {'', 'poltcomp'};
    roinums = [20 50 100 500 1000 5000 10000 20000];
    roitypes = {{'hemiCmkm',''},{'hemiDistKm',''},{'hemiCube4',''}};
    roitypelabels = {'Cm','Dist','Cube4'};
    
    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = []; ii=1; xlabels = {};
        for rr=1:length(roinums)
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                for n=1:length(nuisance)
                    for k=1:length(smooth)
                        if r==3 && roinums(rr) == 5000 % Cube4 case
                            pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} roitypes{r}{2}];
                        else
                            pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(roinums(rr)) roitypes{r}{2}];
                        end
                        xlabels{ii} = [smooth{k} nuisance{n} 'roi' num2str(roinums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
                        if exist(aucmat,'file')
                            load(aucmat);
                        else
                            A = nan(24,100);
                            R = nan(24,1);
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
    I = getR3idx([7 9 19 21],[0 24 48]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9],[0 24 48]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % this is for figure.2 SC-FC corr & detection from 20 to 10000 ROIs.
    for i=[7]
        X = []; K = []; slabels = {};
        for j=1:length(roitypelabels)
            idx = [];
            for k=1:length(roinums)
                [m,k1] = max(R3(i+(j-1)*24,(k-1)*12+1:k*12),[],2);
                [n,k2] = max(A3(i+(j-1)*24,(k-1)*12+1:k*12),[],2);
                idx = (k-1)*12+k1;
                slabels = [slabels(:); [roitypelabels{j} ' ' xlabels{idx}]];
                X = [X; [m,n]];
                K = [K; [k1,k2]];
            end
        end
        figure; plot(X); xticks(1:length(X)); xticklabels(slabels); title(['FC-SC corr & detection results by threshold in (' num2str(i) ') ' vslabels{i}]); setlineColors(6);
    end

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0 24 48]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9], [0 24 48]);
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
    I = getR3idx([7 9 19 21],[0 24 48]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9],[0 24 48]);
    figure; plot(B(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation & detection'); setlineColors(2); setlineStyles({'-','--'});
end

function checkSmoothingNuisanceByDistKm1000vox(vslabels)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'', 's30', 's80', 's150', 's230', 's300'};
    nuisance = {'', 'poltcomp'};
    voxnums = [128 64 32 16 8 4 2 1];
    roitypes = {{'hemiDistKm1000vox',''}};
    roitypelabels = {'Dist'};

    ylabels = {}; R3 = []; A3 = []; AA3 = [];
    for r = 1:length(roitypes)
        Am = []; Rm = []; AA = []; ii=1; xlabels = {};
        for rr=1:length(voxnums)
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
                for n=1:length(nuisance)
                    for k=1:length(smooth)
                        pftype = [smooth{k} hpfstr nuisance{n} preproc roitypes{r}{1} num2str(voxnums(rr)) roitypes{r}{2}];
                        xlabels{ii} = [smooth{k} nuisance{n} 'vox' num2str(voxnums(rr))]; ii=ii+1;
                        aucmat = ['results/auc/' pftype '-fcauc.mat'];
                        if exist(aucmat,'file')
                            load(aucmat);
                        else
                            A = nan(24,100);
                            R = nan(24,1);
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
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(R3(I,:),xlabels,ylabels(I),[0.2 0.9]); colorbar; title(['FC-SC correlation (All) ']); colormap(hot);
 %   figure; plot(R3'); legend(ylabels); title(['FC-SC correlation (All)']); setlineColors(24);
    
    % FC-SC correlation Traced neuron vs synapse
    I = getR3idx([7 9 19 21],[0]);
 %   figure; imagescLabel2(R3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC correlation Full vs. Traced');
    figure; plot(R3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC correlation Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    % FC-SC detection (all)
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I),[0.5 1]); colorbar; title(['FC-SC detection (All) ']); colormap(hot);
%    figure; plot(A3'); legend(ylabels); title(['FC-SC detection (All)']); setlineColors(24);

    % FC-SC detection Traced neuron vs synapse
    I = getR3idx([7 9 19 21], [0]);
%    figure; imagescLabel2(A3(I,:),xlabels,ylabels(I)); colorbar; title('FC-SC detection Full vs. Traced');
    figure; plot(A3(I,:)'); legend(ylabels(I)); xticks(1:ii); xticklabels(xlabels); title('FC-SC detection Traced neuron vs synapse'); setlineColors(2); setlineStyles({'-','--'});

    for i=[7 9]
        X = []; slabels = {};
        for j=1:length(roitypelabels)
            idx = [];
            for k=1:length(voxnums)
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
    I = getR3idx([7 9 19 21],[0]);  % show only Traced neuron, synapse
    figure; imagescLabel2(B(I,:),xlabels,ylabels(I),[0 1.5]); colorbar; title('FC-SC correlation & detection'); colormap(hot);
%    figure; plot(B'); legend(ylabels); title('FC-SC correlation & detection'); setlineColors(24);

    % FC-SC correlation & detection Traced neuron vs synapse (which is best?)
    I = getR3idx([7 9 19 21],[0]);
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

