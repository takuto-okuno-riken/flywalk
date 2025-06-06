% analyze functional connectivity matrix (1st and 2nd level analysis).
% this script should run after makeVoxelROIatlas.m, makeStructConnectivity.m, extractROItimeseries.m first.

function analyzeFuncConnectivity
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = [0]; % high-pass filter threshold
%    hpfTh = [0, 0.1, 0.05, 0.025, 0.02, 0.01, 0.009, 0.008, 0.005, 0.001]; % high-pass filter threshold
%    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80', 's90', 's100', 's110', 's120', 's130', 's140', 's150', 's160', 's170', 's180', 's190', 's200', 's210', 's220', 's230', 's240', 's250', 's260', 's270', 's280', 's290', 's300'};
%    smooth = {'', 's30', 's80'}; %, 's150','s230','s300'};
%    smooth = {'', 's30', 's40', 's60', 's80', 's100', 's150'};
%    smooth = {'s230'}; % DistKm
%    smooth = {''}; % hemiroi
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp'}; % , ...
%        'polacompn','poltcompn','poltacompn'}; % normalized for checking diff with polacomp,etc.
%    nuisance = {'6hmtacomp'}; % good for bransonhemi, branson7065km50
%    nuisance = {'','poltcomp'}; % good for DistKm
%    nuisance = {'','6hm','tcomp','pol','poltcomp'}; % good for hemiRoiXX
%    nuisance = {''};
    nuisance = {'poltcomp'}; % good for hemiroi, DistKm

    % file number setting for random subsampling
    rNums = [0]; % no number setting
%    rNums = 1:99; % for random subsampling number

    % using subjects (flys). sbj 7 shows NaN row in FC matrix
    sbjids = [1 2 3 4 5 6 8 9];

    % ROI name
%    roitypes = {'hemiroif'};  % flyem ROI full
%    roitypes = {'hemiroi','bransonhemi'}; % flyem ROI (Turner compatible)
%    roitypes = {'hemiroi_hb0sr80','hemiroi_fw0','hemiroi_avg0'}; % flyem ROI (Primary, FlyEM vs. FlyWire vs. Average)
%    roitypes = {'hemiBranson7065'};
%    roitypes = {'hemiBranson7065km20','hemiBranson7065km30','hemiBranson7065km50','hemiBranson7065km100','hemiBranson7065km200', ...
%        'hemiBranson7065km300','hemiBranson7065km500','hemiBranson7065km300'};
%    roitypes = {'hemiCmkm10','hemiDistKm10','hemiCmkm10r1w1'};
%    roitypes = {'hemiCmkm20','hemiCmkm30','hemiCmkm50','hemiCmkm100','hemiCmkm200', ...
%        'hemiCmkm300','hemiCmkm500','hemiCmkm1000'};
%    roitypes = {'hemiCmkm20r1w1','hemiCmkm30r1w1','hemiCmkm50r1w1','hemiCmkm100r1w1','hemiCmkm200r1w1', ...
%        'hemiCmkm300r1w1','hemiCmkm500r1w1','hemiCmkm1000r1w1'};
%    roitypes = {'hemiCmkm20r1w1','hemiCmkm30r1w1','hemiCmkm50r1w1','hemiCmkm100r1w1','hemiCmkm200r1w1', ...
%        'hemiCmkm300r1w1','hemiCmkm500r1w1','hemiCmkm1000r1w1'};
    roitypes = {'hemiDistKm20','hemiDistKm30','hemiDistKm50','hemiDistKm100','hemiDistKm200', ...
        'hemiDistKm300','hemiDistKm500','hemiDistKm1000'};
%    roitypes = {'hemiRand20','hemiRand30','hemiRand50','hemiRand100','hemiRand200', ...
%        'hemiRand300','hemiRand500','hemiRand1000'};
%    roitypes = {'hemiVrand20','hemiVrand30','hemiVrand50','hemiVrand100','hemiVrand200', ...
%        'hemiVrand300','hemiVrand500','hemiVrand1000'};
%    roitypes = {'hemiCube4'};
%    roitypes = {'hemiPiece12','hemiPiece8','hemiPiece4'};
%    roitypes = {'hemiPiece3','hemiPiece2'};
    % neuropil FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
%    roitypes = {'hemiRoi101','hemiRoi57','hemiRoi57-51','hemiRoi51-62-20-111-100'};
%    roitypes = {'hemiRoi1','hemiRoi5','hemiRoi7','hemiRoi27','hemiRoi30','hemiRoi32','hemiRoi43','hemiRoi52', ...
%        'hemiRoi54','hemiRoi57','hemiRoi59','hemiRoi63','hemiRoi65','hemiRoi67','hemiRoi78','hemiRoi82', ...
%        'hemiRoi89','hemiRoi93','hemiRoi95','hemiRoi100','hemiRoi101','hemiRoi106','hemiRoi113'};
%    roitypes = {'hemiroi','hemiroi_fw0sr140','hemiBranson7065km50','hemiBranson7065km50_fw0sr140','hemiCmkm50','hemiCmkm50_fw0sr140', ...
%        'hemiCmkm50r1w1','hemiDistKm50','hemiDistKm50_fw0sr140','hemiRand50','hemiVrand50'};
%    roitypes = {'hemiBranson7065km30','hemiCmkm30','hemiCmkm30r1w1','hemiDistKm30','hemiRand30','hemiVrand30'};
%    roitypes = {'hemiCmkm20000','hemiCmkm20000r1w1','hemiDistKm20000','hemiVrand20000'};
%    roitypes = {'hemiCmkm100','hemiDistKm100','hemiCmkm500','hemiDistKm500','hemiCmkm1000','hemiDistKm1000'};
%    roitypes = {'hemiCmkm5000','hemiDistKm5000','hemiCmkm10000','hemiDistKm10000'};
    % FlyEM vs. FlyWire in each roi num by smoothing 0 to 80 (no nuisanse)
%    roitypes = {'hemiBranson7065km20_fw0sr140','hemiBranson7065km30_fw0sr140','hemiBranson7065km50_fw0sr140','hemiBranson7065km100_fw0sr140','hemiBranson7065km200_fw0sr140','hemiBranson7065km300_fw0sr140','hemiBranson7065km500_fw0sr140', 'hemiBranson7065km1000_fw0sr140', ...
%        'hemiCmkm20_fw0sr140','hemiCmkm30_fw0sr140','hemiCmkm50_fw0sr140','hemiCmkm100_fw0sr140','hemiCmkm200_fw0sr140','hemiCmkm300_fw0sr140', 'hemiCmkm500_fw0sr140', 'hemiCmkm1000_fw0sr140',...
%        'hemiDistKm20_fw0sr140','hemiDistKm30_fw0sr140','hemiDistKm50_fw0sr140','hemiDistKm100_fw0sr140','hemiDistKm200_fw0sr140','hemiDistKm300_fw0sr140','hemiDistKm500_fw0sr140','hemiDistKm1000_fw0sr140'};
%    roitypes = {'hemiCmkm50','hemiDistKm50','hemiCmkm100','hemiDistKm100','hemiCmkm500','hemiDistKm500'};  % for large smoothing size & no nuisanse, poltcomp
%    roitypes = {'hemiroi','hemiroi_fw0sr140','hemiDistKm50','hemiDistKm50_fw0sr140','hemiDistKm50_avg'};  % for s0 to s80 (no nuisanse) % for all nuisanse & s30, s80
%    roitypes = {'hemiroi_hb0sr50','hemiroi_hb0sr60','hemiroi_hb0sr70','hemiroi_hb0sr80','hemiroi_hb0sr90', ... 
%            'hemiroi_fw0sr50','hemiroi_fw0sr70','hemiroi_fw0sr100','hemiroi_fw0sr130','hemiroi_fw0sr140','hemiroi_fw0sr150', ...
%            'hemiroi_hb5sr60','hemiroi_fw5sr50','hemiroi_fw5sr140'};  % for s30 & s80, '' & poltcomp
%    roitypes = {'hemiDistKm500_hb0sr50','hemiDistKm500_hb0sr60','hemiDistKm500_hb0sr70','hemiDistKm500_hb0sr80','hemiDistKm500_hb0sr90', ... 
%            'hemiDistKm500_fw0sr50','hemiDistKm500_fw0sr70','hemiDistKm500_fw0sr100','hemiDistKm500_fw0sr130','hemiDistKm500_fw0sr140','hemiDistKm500_fw0sr150', ...
%            'hemiDistKm500_fw5sr50','hemiDistKm500_fw5sr140'};  % for s30 & s80, '' & poltcomp
%    roitypes = {'hemiRoi68-59-87-106-50-27-54'};  % s30,80,150, 6hm,pol,tcomp,poltcomp
%    roitypes = {'hemiRoi68-59-87-106-50-27-54DistKm200'};
%    roitypes = {'hemiDistKm1000vox128','hemiDistKm1000vox64','hemiDistKm1000vox32','hemiDistKm1000vox16','hemiDistKm1000vox8','hemiDistKm1000vox4','hemiDistKm1000vox2','hemiDistKm1000vox1'};
%    roitypes = {'hemiroi_hb0sr80_sp10db3000mi1','hemiroi_hb0sr80_sp20db3000mi1','hemiroi_hb0sr80_sp40db3000mi1','hemiroi_hb0sr80_sp60db3000mi1','hemiroi_hb0sr80_sp80db3000mi1','hemiroi_hb0sr80_sp90db3000mi1', ...
%            'hemiroi_fw0sr140_sp10db3000mi1','hemiroi_fw0sr140_sp20db3000mi1','hemiroi_fw0sr140_sp40db3000mi1','hemiroi_fw0sr140_sp60db3000mi1','hemiroi_fw0sr140_sp80db3000mi1','hemiroi_fw0sr140_sp90db3000mi1'};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80_sp10db3000mi1_only1','hemiroi_hb0sr80_sp20db3000mi1_only1','hemiroi_hb0sr80_sp40db3000mi1_only1','hemiroi_hb0sr80_sp60db3000mi1_only1','hemiroi_hb0sr80_sp80db3000mi1_only1','hemiroi_hb0sr80_sp90db3000mi1_only1', ...
%            'hemiroi_fw0sr140_sp10db3000mi1_only1','hemiroi_fw0sr140_sp20db3000mi1_only1','hemiroi_fw0sr140_sp40db3000mi1_only1','hemiroi_fw0sr140_sp60db3000mi1_only1','hemiroi_fw0sr140_sp80db3000mi1_only1','hemiroi_fw0sr140_sp90db3000mi1_only1'};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80_rc20','hemiroi_hb0sr80_rc40','hemiroi_hb0sr80_rc100','hemiroi_hb0sr80_rc500', 'hemiroi_hb0sr80_rc1000','hemiroi_hb0sr80_rc10000', ...
%            'hemiroi_fw0sr140_rc20','hemiroi_fw0sr140_rc40','hemiroi_fw0sr140_rc100','hemiroi_fw0sr140_rc500','hemiroi_fw0sr140_rc1000','hemiroi_fw0sr140_rc10000'};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80_sp5db3000mi1_rc40','hemiroi_hb0sr80_sp5db3000mi1_rc10000', ...
%            'hemiroi_fw0sr140_sp5db3000mi1_rc40','hemiroi_fw0sr140_sp5db3000mi1_rc10000'};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80_rc10000_rand1','hemiroi_hb0sr80_rc10000_rand2','hemiroi_hb0sr80_rc10000_rand3', ...
%            'hemiroi_hb0sr80_rc10000_xrand1','hemiroi_hb0sr80_rc10000_xrand2','hemiroi_hb0sr80_rc10000_xrand3',...
%            'hemiroi_hb0sr80_rc10000_xorand1','hemiroi_hb0sr80_rc10000_xorand2','hemiroi_hb0sr80_rc10000_xorand3', ...
%            'hemiroi_fw0sr140_rc10000_rand1','hemiroi_fw0sr140_rc10000_rand2','hemiroi_fw0sr140_rc10000_rand3', ...
%            'hemiroi_fw0sr140_rc10000_xrand1','hemiroi_fw0sr140_rc10000_xrand2','hemiroi_fw0sr140_rc10000_xrand3', ...
%            'hemiroi_fw0sr140_rc10000_xorand1','hemiroi_fw0sr140_rc10000_xorand2','hemiroi_fw0sr140_rc10000_xorand3'};  % for s0, poltcomp
%    roitypes = {'hemidistkm500_hb0sr80_rc20_xorand1','hemidistkm500_hb0sr80_rc20_xorand2','hemidistkm500_hb0sr80_rc20_xorand3', ...
%            'hemidistkm500_fw0sr140_rc20_xorand1','hemidistkm500_fw0sr140_rc20_xorand2','hemidistkm500_fw0sr140_rc20_xorand3'};  % for s230, poltcomp
%    roitypes = {'hemiroi_hb0sr80_rc20_only1','hemiroi_hb0sr80_rc40_only1','hemiroi_hb0sr80_rc100_only1','hemiroi_hb0sr80_rc500_only1','hemiroi_hb0sr80_rc1000_only1','hemiroi_hb0sr80_rc10000_only1', ...
%            'hemiroi_fw0sr140_rc20_only1','hemiroi_fw0sr140_rc40_only1','hemiroi_fw0sr140_rc100_only1','hemiroi_fw0sr140_rc500_only1','hemiroi_fw0sr140_rc1000_only1','hemiroi_fw0sr140_rc10000_only1'};  % for s0, poltcomp
%    roitypes = {'hemidistkm500_hb0sr80_rc20_only1','hemidistkm500_hb0sr80_rc40_only1','hemidistkm500_hb0sr80_rc100_only1','hemidistkm500_hb0sr80_rc500_only1','hemidistkm500_hb0sr80_rc1000_only1','hemidistkm500_hb0sr80_rc10000_only1', ...
%            'hemidistkm500_fw0sr140_rc20_only1','hemidistkm500_fw0sr140_rc40_only1','hemidistkm500_fw0sr140_rc100_only1','hemidistkm500_fw0sr140_rc500_only1','hemidistkm500_fw0sr140_rc1000_only1','hemidistkm500_fw0sr140_rc10000_only1'};  % for s230, poltcomp
%    roitypes = {'hemiroi_hb0sr80_rn50_orand1','hemiroi_hb0sr80_rn150_orand1','hemiroi_hb0sr80_rn500_orand1','hemiroi_hb0sr80_rn1000_orand1', ...
%            'hemiroi_hb0sr80fw_rn50_orand1','hemiroi_hb0sr80fw_rn60_orand1','hemiroi_hb0sr80fw_rn70_orand1','hemiroi_hb0sr80fw_rn130_orand1','hemiroi_hb0sr80fw_rn150_orand1', ...
%            'hemiroi_hb0sr80fw_rn500_orand1','hemiroi_hb0sr80fw_rn600_orand1','hemiroi_hb0sr80fw_rn700_orand1','hemiroi_hb0sr80fw_rn710_orand1','hemiroi_hb0sr80fw_rn1000_orand1', ...
%            'hemiroi_fw0sr140_rn50_orand1','hemiroi_fw0sr140_rn120_orand1','hemiroi_fw0sr140_rn130_orand1','hemiroi_fw0sr140_rn140_orand1','hemiroi_fw0sr140_rn150_orand1','hemiroi_fw0sr140_rn200_orand1','hemiroi_fw0sr140_rn250_orand1','hemiroi_fw0sr140_rn500_orand1', ...
%            'hemiroi_fw0sr140_rn1000_orand1','hemiroi_fw0sr140_rn1530_orand1','hemiroi_fw0sr140_rn1550_orand1','hemiroi_fw0sr140_rn2000_orand1',};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80fw_rc20_xorand1','hemiroi_hb0sr80fw_rc20_xorand2','hemiroi_hb0sr80fw_rc20_xorand3', ...
%            'hemiroi_hb0sr80_rc20_xorand1','hemiroi_hb0sr80_rc20_xorand2','hemiroi_hb0sr80_rc20_xorand3', ...
%            'hemiroi_fw0sr140_rc20_xorand1','hemiroi_fw0sr140_rc20_xorand2','hemiroi_fw0sr140_rc20_xorand3'};  % for s0, poltcomp
%    roitypes = {'hemiroi_hb0sr80fw_rd67-12','hemiroi_hb0sr80fw_rd140-15','hemiroi_hb0sr80fw_rd705-40', ...
%            'hemiroi_fw0sr140_rd140-15','hemiroi_fw0sr140_rd185-20','hemiroi_fw0sr140_rd1560-80'};  % for s0, poltcomp, rNums

    %to check scatter and AUC graph (need to comment out plot lines)
%    roitypes = {'hemiroi_hb0sr80','hemiroi_hb0sr80_rc20_only1','hemiroi_hb0sr80_rn150_orand1','hemiroi_hb0sr80_rn10_orand1', ...
%            'hemiroi_fw0sr140','hemiroi_fw0sr140_rc20_only1','hemiroi_fw0sr140_rn10_orand1'};  % for s0, poltcomp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:length(roitypes)
        for r = 1:length(rNums)
            rstr = '';
            if rNums(r) > 0, rstr = ['-' num2str(rNums(r))]; end
            analyzeFcROItype([roitypes{n} rstr], preproc, hpfTh, smooth, nuisance, sbjids, false)
        end
    end
end

function analyzeFcROItype(roitype, preproc, hpfTh, smooth, nuisance, sbjids, isplot)
    AUCVER = 3;

    % load structural connectivity matrix (from makeStructConnectivity.m)
    switch(roitype)
    case 'hemiroif'
        load('results/sc/hemiroi_connectlist.mat');
        ids = 1:roiNum;
    otherwise
        roitype = lower(roitype);
        load(['results/sc/' roitype '_connectlist.mat']);
        ids = primaryIds;
    end

    nweightMat(isnan(nweightMat)) = 0;
    isw2 = ~isempty(nweightMat);
    issw = ~isempty(syweightMat);
    C2 = ncountMat(ids,ids,1); S = sycountMat(ids,ids,1); W2 = []; Wo = []; Sw = []; W3 = [];
    if isw2
        outweightMat(isnan(outweightMat)) = 0;
        W3 = nweightMat(ids,ids,1); Wo = outweightMat(ids,ids,1);
        W2 = W3 .* S; % pure ROI-input neuron connection weight
    end
    if issw
        Sw = syweightMat(ids,ids,1);
    end
    C2b = ncountMat(ids,ids,2); Sb = sycountMat(ids,ids,2); W2b = []; Wob = []; Swb = []; W3b = [];
    if isw2
        W3b = nweightMat(ids,ids,2); Wob = outweightMat(ids,ids,2);
        W2b = W3b .* Sb; % pure ROI-input neuron connection weight
    end
    if issw
        Swb = syweightMat(ids,ids,2);
    end

    % show corr between neurons v. synapse weight
    if isplot
        if isw2
            r = corr(W2b(:),C2b(:));    % corr between neurons v. synapse weight
            disp(['corr between synapse weight vs. neurons. r=' num2str(r)]);
            figure; scatter(W2b(:),C2b(:)); xlabel('synapse weight2'); ylabel('neurons2');
        end
        r = corr(Sb(:),C2b(:));    % corr between neurons v. synapse weight
        disp(['corr between synapses vs. neurons. r=' num2str(r)]);
        figure; scatter(Sb(:),C2b(:)); xlabel('synapses'); ylabel('neurons2');
    end

    n = length(ids);
    E = eye(n); E = logical(1-E);

    lC2 = log10(C2); lC2(lC2<0) = 0;
    lS = log10(S); lS(lS<0) = 0;
    lW2 = log10(W2); lW2(lW2<0) = 0;
    lC2b = log10(C2b); lC2b(lC2b<0) = 0;
    lSb = log10(Sb); lSb(lSb<0) = 0;
    lW2b = log10(W2b); lW2b(lW2b<0) = 0;

%    figure; imagesc(lC2); colorbar; title('lC2'); % to check ordered matrix

    if ~exist('results/auc','dir'), mkdir('results/auc'); end
    if ~exist('results/fc','dir'), mkdir('results/fc'); end

    sbjR = [];
    Rm = []; SRm = []; rlabel = {}; ii=1;
    AUC = [];
    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                pftype = [smooth{k} hpfstr nuisance{n} preproc roitype];
                rlabel{ii} = [smooth{k} hpfstr nuisance{n}]; ii=ii+1;

                % load ROI time-series (from extractROItimeseries.m)
                str = split(pftype,'_');
                load(['results/ts/' str{1} '-ts.mat']); % FlyWire SC calculation.
    
                CM = {};
                for i=1:length(sbjids)
                    CX{sbjids(i)} = CX{sbjids(i)}';
                    CM{i} = single(corr(CX{sbjids(i)})); %calcPartialCorrelation_(CX{i}',[],[],[],false,1e-3);
                end
                cmlen = length(CM);
                    
                % transport first
                outfile = ['results/fc/' pftype '-func.mat'];
                if exist(outfile,'file')
                    load(outfile);
                else
                    [B2, RSS2, T2, df] = calcSeedCorrMixed(CX(sbjids));
                    
                    % output beta matrix
                    save(outfile,'T2','-v7.3');
                end

                F3 = []; F3z = [];
                for i=1:cmlen
                    F = CM{i}(ids,ids);
                    F3 = cat(3,F3,F);
                end
                F3z = atanh(F3); % z transformed (better FC-SC corr).
                clear CM;

                % mean group data
                mFz = nanmean(F3z,3);
                mFz(isinf(mFz)) = max(mFz(~isinf(mFz)));
%                figure; imagesc(abs(mFz)); colorbar; daspect([1 1 1]); title([pftype ' m-FCz matrix']);

                T3 = T2(ids,ids);
                T3(isinf(T3)) = max(T3(~isinf(T3)));
%                figure; imagesc(abs(T3)); colorbar; daspect([1 1 1]); title([pftype ' T-value matrix']);
                lT3 = log(T3); lT3(lT3<0) = 0;

                % each flys [log10(neurons) vs. m-FCz]
                for i=1:cmlen
                    Dz = F3z(:,:,i);
                    Dz(isinf(Dz)) = max(Dz(~isinf(Dz)));
                    sbjR(k,i) = corr(lC2(:),abs(Dz(:)));
                end

                % each ROIs (Traced vs. m-FCz and only)
                roiR = nan(size(lC2,1),24,'single');
                for i=1:size(lC2,1)
                    roimFz = abs([mFz(i,:)';mFz(:,i)]);
                    roiR(i,7) = corr([lC2b(i,:)';lC2b(:,i)],roimFz);
                    roiR(i,9) = corr([lSb(i,:)';lSb(:,i)],roimFz);
                    if isw2
                        roiR(i,8) = corr([lW2b(i,:)';lW2b(:,i)],roimFz);
                        roiR(i,10) = corr([W3b(i,:)';W3b(:,i)],roimFz);
                        roiR(i,12) = corr([Wob(i,:)';Wob(:,i)],roimFz);
                    end
                    if issw
                        roiR(i,11) = corr([Swb(i,:)';Swb(:,i)],roimFz);
                    end
                end

                % full ROIs (vs. mean group data)
                R = nan(1,24,'single');
                SR = nan(1,24,'single');
                R(1) = corr(lC2(:),abs(mFz(:))); % whole vs. m-FCz
                disp(['prefix=' pftype ' : log10(neurons2) vs. m-FCz = ' num2str(R(1))]);
                R(3) = corr(lS(:),abs(mFz(:)));
                disp(['prefix=' pftype ' : log10(synapses) vs. m-FCz = ' num2str(R(3))]);
                if isw2
                    R(2) = corr(lW2(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : log10(synapse weight2) vs. m-FCz = ' num2str(R(2))]);
                    R(4) = corr(W3(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : ROI in-neuron weight vs. m-FCz = ' num2str(R(4))]);
                    if issw 
                        R(5) = corr(Sw(:),abs(mFz(:)));
                        disp(['prefix=' pftype ' : ROI in-synapse weight vs. m-FCz = ' num2str(R(5))]);
                    end
                    R(6) = corr(Wo(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : ROI out-neuron weight vs. m-FCz = ' num2str(R(6))]);
                end

                R(7) = corr(lC2b(:),abs(mFz(:))); % Traced vs. m-FCz
                disp(['prefix=' pftype ' : log10(neurons2b) vs. m-FCz = ' num2str(R(7))]);
%                figure; scatter(lC2b(:),abs(mFz(:)),18,'x'); xlabel('log10(neurons2b)'); ylabel('m-FCz'); title([pftype ' r='num2str(R(7))]); ylim([0 3]);
                R(9) = corr(lSb(:),abs(mFz(:)));
                disp(['prefix=' pftype ' : log10(synapses b) vs. m-FCz = ' num2str(R(9))]);
%                figure; scatter(lSb(:),abs(mFz(:)),18,'x'); xlabel('log10(synapses b)'); ylabel('m-FCz'); title([pftype ' r=' num2str(R(9))]);
                SR(7) = corr(lC2b(:),abs(mFz(:)),'Type','Spearman'); % Traced vs. m-FCz
                SR(9) = corr(lSb(:),abs(mFz(:)),'Type','Spearman');
                if isw2
                    R(8) = corr(lW2b(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : log10(synapse weight2b) vs. m-FCz = ' num2str(R(8))]);
                    R(10) = corr(W3b(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : ROI in-neuron weight b vs. m-FCz = ' num2str(R(10))]);
                    if issw 
                        R(11) = corr(Swb(:),abs(mFz(:)));
                        disp(['prefix=' pftype ' : ROI in-synapse weight b vs. m-FCz = ' num2str(R(11))]);
                    end
                    R(12) = corr(Wob(:),abs(mFz(:)));
                    disp(['prefix=' pftype ' : ROI out-neuron weight b vs. m-FCz = ' num2str(R(12))]);
                end

                R(13) = corr(lC2(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(neurons2) vs. FC-Tval = ' num2str(R(13))]);
                R(15) = corr(lS(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapses) vs. FC-Tval = ' num2str(R(15))]);
                if isw2
                    R(14) = corr(lW2(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : log10(synapse weight2) vs. FC-Tval = ' num2str(R(14))]);
                    R(16) = corr(W3(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : ROI in-neuron weight vs. FC-Tval = ' num2str(R(16))]);
%                    figure; scatter(W3(:),abs(T3(:))); xlabel('ROI in-neuron weight'); ylabel('FC-Tval'); title([pftype ' r=' num2str(R(16))]);
                    if issw 
                        R(17) = corr(Sw(:),abs(T3(:)));
                        disp(['prefix=' pftype ' : ROI in-synapse weight vs. FC-Tval = ' num2str(R(17))]);
                    end
                    R(18) = corr(Wo(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : ROI out-neuron weight vs. FC-Tval = ' num2str(R(18))]);
                end

                R(19) = corr(lC2b(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(neurons2b) vs. FC-Tval = ' num2str(R(19))]);
                R(21) = corr(lSb(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapses b) vs. FC-Tval = ' num2str(R(21))]);
                SR(19) = corr(lC2b(:),abs(T3(:)),'Type','Spearman');
                SR(21) = corr(lSb(:),abs(T3(:)),'Type','Spearman');
                if isw2
                    R(20) = corr(lW2b(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : log10(synapse weight2b) vs. FC-Tval = ' num2str(R(20))]);
                    R(22) = corr(W3b(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : ROI in-neuron weight b vs. FC-Tval = ' num2str(R(22))]);
                    if issw 
                        R(23) = corr(Swb(:),abs(T3(:)));
                        disp(['prefix=' pftype ' : ROI in-synapse weight b vs. FC-Tval = ' num2str(R(23))]);
                    end
                    R(24) = corr(Wob(:),abs(T3(:)));
                    disp(['prefix=' pftype ' : ROI out-neuron weight b vs. FC-Tval = ' num2str(R(24))]);
                end
                Rm = [Rm, R']; SRm = [Rm, SR']; aucver = 3;
                
                disp(['prefix=' pftype ' : log10(neurons2b) vs. m-FCz (sp)= ' num2str(SR(7))]);
                disp(['prefix=' pftype ' : log10(synapses b) vs. m-FCz (sp)= ' num2str(SR(9))]);
                disp(['prefix=' pftype ' : log10(neurons2b) vs. FC-Tval (sp)= ' num2str(SR(19))]);
                disp(['prefix=' pftype ' : log10(synapses b) vs. FC-Tval (sp)= ' num2str(SR(21))]);

                % calculate AUC
                aucmat = ['results/auc/' pftype '-fcauc.mat'];
                if exist(aucmat,'file')
                    % load beta volumes
                    load(aucmat);
                else
                    thN = 100;
                    aths = cell(thN,1);
                    XY = cell(thN,1);
                    sbths = [];
%                    for th = 1:thN
                    parfor th = 1:thN
                        % include injection voxel in ground truth
                        c2th = prctile(C2(C2>0),th-1);
                        ct2 = C2; ct2(ct2<c2th) = 0; ct2(ct2>0) = 1;
                        sth = prctile(S(S>0),th-1);
                        st = S; st(st<sth) = 0; st(st>0) = 1;
                        w2th = prctile(W2(W2>0),th-1);
                        wt2 = W2; wt2(wt2<w2th) = 0; wt2(wt2>0) = 1;
                        w3th = prctile(W3(W3>0),th-1);
                        wt3 = W3; wt3(wt3<w3th) = 0; wt3(wt3>0) = 1;
                        swth = prctile(Sw(Sw>0),th-1);
                        swt = Sw; swt(swt<swth) = 0; swt(swt>0) = 1;
                        woth = prctile(Wo(Wo>0),th-1);
                        wot = Wo; wot(wot<woth) = 0; wot(wot>0) = 1;

                        c2bth = prctile(C2b(C2b>0),th-1);
                        ct2b = C2b; ct2b(ct2b<c2bth) = 0; ct2b(ct2b>0) = 1;
                        sbth = prctile(Sb(Sb>0),th-1);
                        stb = Sb; stb(stb<sbth) = 0; stb(stb>0) = 1; sbths(th) = sbth;
                        w2bth = prctile(W2b(W2b>0),th-1);
                        wt2b = W2b; wt2b(wt2b<w2bth) = 0; wt2b(wt2b>0) = 1;
                        w3bth = prctile(W3b(W3b>0),th-1);
                        wt3b = W3b; wt3b(wt3b<w3bth) = 0; wt3b(wt3b>0) = 1;
                        swbth = prctile(Swb(Swb>0),th-1);
                        swtb = Swb; swtb(swtb<swbth) = 0; swtb(swtb>0) = 1;
                        wobth = prctile(Wob(Wob>0),th-1);
                        wotb = Wob; wotb(wotb<wobth) = 0; wotb(wotb>0) = 1;

                        aucs = cell(24,1);
                        [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(mFz(:)'), ['m-FCz vs. neurons2 th=' num2str(th-1)], false);
%                        figure; imagesc(ct2); colorbar; figure; imagesc(abs(mFz)); colorbar; % to check matrix
                        aucs{1} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(mFz(:)'), ['m-FCz vs. synapses th=' num2str(th-1)], false);
                        aucs{3} = single(auc);
                        if isw2
                            [~, ~, auc] = calcShowGroupROCcurve(wt2(:)', abs(mFz(:)'), ['m-FCz vs. synapse weight2 th=' num2str(th-1)], false);
                            aucs{2} = single(auc);
                            [~, ~, auc] = calcShowGroupROCcurve(wt3(:)', abs(mFz(:)'), ['m-FCz vs. ROI in-neuron weight th=' num2str(th-1)], false);
                            aucs{4} = single(auc);
                            if issw 
                                [~, ~, auc] = calcShowGroupROCcurve(swt(:)', abs(mFz(:)'), ['m-FCz vs. ROI in-synapse weight th=' num2str(th-1)], false);
                                aucs{5} = single(auc);
                            end
                            [~, ~, auc] = calcShowGroupROCcurve(wot(:)', abs(mFz(:)'), ['m-FCz vs. ROI out-neuron weight th=' num2str(th-1)], false);
                            aucs{6} = single(auc);
                        end

                        [X, Y, auc] = calcShowGroupROCcurve(ct2b(:)', abs(mFz(:)'), ['m-FCz vs. neurons2b th=' num2str(th-1)], false);
                        aucs{7} = single(auc); % input node is one, so it is not vector.
                        [X, Y, auc] = calcShowGroupROCcurve(stb(:)', abs(mFz(:)'), ['m-FCz vs. synapses b th=' num2str(th-1)], false);
                        aucs{9} = single(auc);
                        if isw2
                            [~, ~, auc] = calcShowGroupROCcurve(wt2b(:)', abs(mFz(:)'), ['m-FCz vs. synapse weight2b th=' num2str(th-1)], false);
                            aucs{8} = single(auc);
                            [~, ~, auc] = calcShowGroupROCcurve(wt3b(:)', abs(mFz(:)'), ['m-FCz vs. ROI in-neuron weight b th=' num2str(th-1)], false);
                            aucs{10} = single(auc);
                            if issw 
                                [~, ~, auc] = calcShowGroupROCcurve(swtb(:)', abs(mFz(:)'), ['m-FCz vs. ROI in-synapse weight b th=' num2str(th-1)], false);
                                aucs{11} = single(auc);
                            end
                            [~, ~, auc] = calcShowGroupROCcurve(wotb(:)', abs(mFz(:)'), ['m-FCz vs. ROI out-neuron weight b th=' num2str(th-1)], false);
                            aucs{12} = single(auc);
                        end

                        [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(T3(:)'), ['FC-Tval vs. neurons2 th=' num2str(th-1)], false);
                        aucs{13} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(T3(:)'), ['FC-Tval vs. synapses th=' num2str(th-1)], false);
                        aucs{15} = single(auc);
                        if isw2
                            [~, ~, auc] = calcShowGroupROCcurve(wt2(:)', abs(T3(:)'), ['FC-Tval vs. synapse weight2 th=' num2str(th-1)], false);
                            aucs{14} = single(auc);
                            [~, ~, auc] = calcShowGroupROCcurve(wt3(:)', abs(T3(:)'), ['FC-Tval vs. ROI in-neuron weight th=' num2str(th-1)], false);
                            aucs{16} = single(auc);
                            if issw 
                                [~, ~, auc] = calcShowGroupROCcurve(swt(:)', abs(T3(:)'), ['FC-Tval vs. ROI in-synapse weight th=' num2str(th-1)], false);
                                aucs{17} = single(auc);
                            end
                            [~, ~, auc] = calcShowGroupROCcurve(wot(:)', abs(T3(:)'), ['FC-Tval vs. ROI out-neuron weight th=' num2str(th-1)], false);
                            aucs{18} = single(auc);
                        end

                        [~, ~, auc] = calcShowGroupROCcurve(ct2b(:)', abs(T3(:)'), ['FC-Tval vs. neurons2b th=' num2str(th-1)], false);
                        aucs{19} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(stb(:)', abs(T3(:)'), ['FC-Tval vs. synapses b th=' num2str(th-1)], false);
                        aucs{21} = single(auc);
                        if isw2
                            [~, ~, auc] = calcShowGroupROCcurve(wt2b(:)', abs(T3(:)'), ['FC-Tval vs. synapse weight2b th=' num2str(th-1)], false);
                            aucs{20} = single(auc);
                            [~, ~, auc] = calcShowGroupROCcurve(wt3b(:)', abs(T3(:)'), ['FC-Tval vs. ROI in-neuron weight b th=' num2str(th-1)], false);
                            aucs{22} = single(auc);
                            if issw 
                                [~, ~, auc] = calcShowGroupROCcurve(swtb(:)', abs(T3(:)'), ['FC-Tval vs. ROI in-synapse weight b th=' num2str(th-1)], false);
                                aucs{23} = single(auc);
                            end
                            [~, ~, auc] = calcShowGroupROCcurve(wotb(:)', abs(T3(:)'), ['FC-Tval vs. ROI out-neuron weight b th=' num2str(th-1)], false);
                            aucs{24} = single(auc);
                        end
                        aths{th} = aucs;
%                        XY{th} = [X; Y]; % keep X,Y for ROC curve plot
                    end
%{
                    figure; plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]); ylim([0 1]); xlim([0 1]); daspect([1 1 1]);
                    for th = 1:thN
                        c = 0.95 - th*0.8/100;
                        hold on; plot(XY{th}(1,:),XY{th}(2,:),'Color',[c c c 0.7]); hold off;
                    end
                    title([pftype 'log10(synapses b) vs. m-FCz ROC curve']); xlabel('False Positive Rate'); ylabel('True Positive Rate');
%}
%{
                    figure; plot(sbths); xlabel('percentile'); ylabel('threshold number');
                    title([pftype ' (synapses b) thresholds']);

                    figure; histogram(Sb(Sb>0)); ylabel('element count');
                    title([pftype ' (synapses b) histogram']);
%}
                    A = nan(24,thN,'single');
                    for th = 1:thN
                        for j=1:24
                            if ~isempty(aths{th}{j})
                                A(j,th) = aths{th}{j}; % input node is one, so it is not vector.
                            end
                        end
                    end
                end
                if aucver <= AUCVER
                    aucver = aucver + 0.1;
                    save(aucmat,'A','R','SR','roiR','aucver','-v7.3');
                end
                AUC = cat(3,AUC,A);

%                I=[7 9]; 
%                figure; plot(A(I,:)'); title([pftype ' : FC AUC result vs. ground truth SC']);
%                xlabel('percentile'); legend({'log10(neurons) vs. m-FCz','log10(synapses) vs. m-FCz'});
            end
        end
    end
    if ~isplot, return; end

    % FC-SC correlation (6-type mixed box plot)
    dlabels = {'neuron vs. m-FC(z)','synapse vs. m-FC(z)'};
    figure; boxplot(Rm([7 9],:),'Labels',rlabel); title([roitype ' FC-SC correlation (neuron vs. synapse)']);
    hold on; plot(Rm([7 9],:)'); hold off; legend(dlabels);

    % FC-SC detection (m-FCz vs. neurons)
    % neurons: only ROI was transformed. neurons2: synapse points were transformed and re-counted in all ROI. 
%{
    A1 = squeeze(AUC(1,:,:)); % for internal check
    figure; boxplot(A1,'Labels',rlabel); title([roitype ' FC-SC detection (m-FCz vs. neurons)']);
    A5 = squeeze(AUC(5,:,:));
    figure; boxplot(A5,'Labels',rlabel); title([roitype ' FC-SC detection (m-FCz vs. neurons2)']);

    A2 = squeeze(AUC(2,:,:)); % for internal check
    figure; boxplot(A2,'Labels',rlabel); title([roitype ' FC-SC detection (m-FCz vs. synapse weight)']);
    A6 = squeeze(AUC(6,:,:));
    figure; boxplot(A6,'Labels',rlabel); title([roitype ' FC-SC detection (m-FCz vs. synapse weight2)']);
%}
    AA = squeeze(nanmean(AUC,2));
    figure; boxplot(AA([7 9],:),'Labels',rlabel); title([roitype ' FC-SC detection  (neuron vs. synapse)']);
    hold on; plot(AA([7 9],:)'); hold off; legend(dlabels);
end
