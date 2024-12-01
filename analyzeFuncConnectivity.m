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
    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
%    nuisance = {'6hmtacomp'}; % good for bransonhemi, branson7065km50
    nuisance = {''};

    % using subjects (flys). sbj 7 shows NaN row in FC matrix
    sbjids = [1 2 3 4 5 6 8 9];

    % ROI name
%    roitypes = {'flyemroif'};  % flyem ROI full
%    roitypes = {'flyemroi','bransonhemi'}; % flyem ROI (Turner compatible)
%    roitypes = {'hemiBranson7065'};
%    roitypes = {'hemiBranson7065km20','hemiBranson7065km30','hemiBranson7065km50','hemiBranson7065km100','hemiBranson7065km200', ...
%        'hemiBranson7065km300','hemiBranson7065km500','hemiBranson7065km1000'};
%    roitypes = {'hemiCmkm20','hemiCmkm30','hemiCmkm50','hemiCmkm100','hemiCmkm200', ...
%        'hemiCmkm300','hemiCmkm500','hemiCmkm1000'};
%    roitypes = {'hemiCmkm20r1w1','hemiCmkm30r1w1','hemiCmkm50r1w1','hemiCmkm100r1w1','hemiCmkm200r1w1', ...
%        'hemiCmkm300r1w1','hemiCmkm500r1w1','hemiCmkm1000r1w1'};
%    roitypes = {'hemiDistKm20','hemiDistKm30','hemiDistKm50','hemiDistKm100','hemiDistKm200', ...
%        'hemiDistKm300','hemiDistKm500','hemiDistKm1000'};
%    roitypes = {'hemiRand20','hemiRand30','hemiRand50','hemiRand100','hemiRand200', ...
%        'hemiRand300','hemiRand500','hemiRand1000'};
%    roitypes = {'hemiCube12','hemiCube8','hemiCube4'};
%    roitypes = {'hemiPiece12','hemiPiece8','hemiPiece4'};
%    roitypes = {'hemiPiece3','hemiPiece2'};
    % neuropil FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
%    roitypes = {'hemiRoi101','hemiRoi57','hemiRoi57-51','hemiRoi51-62-20-111-100'};
%    roitypes = {'hemiRoi1','hemiRoi5','hemiRoi7','hemiRoi27','hemiRoi30','hemiRoi32','hemiRoi43','hemiRoi52', ...
%        'hemiRoi54','hemiRoi57','hemiRoi59','hemiRoi63','hemiRoi65','hemiRoi67','hemiRoi78','hemiRoi82', ...
%        'hemiRoi89','hemiRoi93','hemiRoi95','hemiRoi100','hemiRoi101','hemiRoi106','hemiRoi113'};
    roitypes = {'flyemroi','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:length(roitypes)
        analyzeFcROItype(roitypes{n}, preproc, hpfTh, smooth, nuisance, sbjids)
    end
end

function analyzeFcROItype(roitype, preproc, hpfTh, smooth, nuisance, sbjids)

    % load structural connectivity matrix (from makeStructConnectivity.m)
    switch(roitype)
    case 'flyemroi'
        load('data/neuprint_connectlist.mat');
        ids = primaryIds;
    case 'flyemroif'
        load('data/neuprint_connectlist.mat');
        ids = 1:roiNum;
    case 'bransonhemi'
        load('data/branson_connectlist.mat');
        ids = primaryIds;
    otherwise
        roitype = lower(roitype);
        load(['data/' roitype '_connectlist.mat']);
        ids = primaryIds;
        weightMat2(isnan(weightMat2)) = 0;
    end
    C2 = countMat2(ids,ids,1); S = sycountMat(ids,ids,1); W2 = weightMat2(ids,ids,1); Wo = outweightMat(ids,ids,1); Sw = syweightMat(ids,ids,1);
    W3 = W2 ./ S; W3(S==0) = 0; % pure ROI-input neuron connection weight
    C2b = countMat2(ids,ids,2); Sb = sycountMat(ids,ids,2); W2b = weightMat2(ids,ids,2); Wob = outweightMat(ids,ids,2); Swb = syweightMat(ids,ids,2);
    W3b = W2b ./ Sb; W3b(Sb==0) = 0; % pure ROI-input neuron connection weight

    % show corr between cell count v. synapse weight
    r = corr(W2(:),C2(:));    % corr between cell count v. synapse weight
    disp(['corr between synapse weight2 vs. cell count2. r=' num2str(r)]);
    figure; scatter(W2(:),C2(:)); xlabel('synapse weight2'); ylabel('cell count2');
%{
    r = corr(C(:),C2(:));    % corr between cell count v. synapse weight (for internal check)
    disp(['corr between cell count vs. cell count2. r=' num2str(r)]);
    figure; scatter(C2(:),C(:)); xlabel('cell count 2'); ylabel('cell count')
%}
    r = corr(S(:),C2(:));    % corr between cell count v. synapse weight
    disp(['corr between synapse count vs. cell count2. r=' num2str(r)]);
    figure; scatter(S(:),C2(:)); xlabel('synapse count'); ylabel('cell count2');

    n = length(ids);
    E = eye(n); E = logical(1-E);

    lC2 = log10(C2); lC2(lC2<0) = 0;
    lS = log10(S); lS(lS<0) = 0;
    lW2 = log10(W2); lW2(lW2<0) = 0;
    lC2b = log10(C2b); lC2b(lC2b<0) = 0;
    lSb = log10(Sb); lSb(lSb<0) = 0;
    lW2b = log10(W2b); lW2b(lW2b<0) = 0;

%    figure; imagesc(lC2); colorbar; title('lC2'); % to check ordered matrix

    sbjR = [];
    roiR = [];
    Rm = []; rlabel = {}; ii=1;
    AUC = [];
    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                pftype = [smooth{k} hpfstr nuisance{n} preproc roitype];
                rlabel{ii} = [smooth{k} hpfstr nuisance{n}]; ii=ii+1;

                % load ROI time-series (from extractROItimeseries.m)
                load(['results/' pftype '-ts.mat']);
    
                CM = {};
                for i=1:length(sbjids)
                    CX{sbjids(i)} = CX{sbjids(i)}';
                    CM{i} = single(corr(CX{sbjids(i)})); %calcPartialCorrelation_(CX{i}',[],[],[],false,1e-3);
                end
                    
                % transport first
                outfile = ['results/' pftype '-func.mat'];
                if exist(outfile,'file')
                    load(outfile);
                else
                    [B2, RSS2, T2, df] = calcSeedCorrMixed(CX(sbjids));
                    
                    % output beta matrix
                    save(outfile,'T2','-v7.3');
                end

                D3 = []; D3z = [];
                for i=1:length(CM)
                    D = CM{i}(ids,ids);
                    D3 = cat(3,D3,D);
                end
                D3z = atanh(D3); % z transformed (better FC-SC corr).

                % mean group data
                Dm = nanmean(D3,3);
                Dmz = nanmean(D3z,3);
                Dmz(isinf(Dmz)) = max(Dmz(~isinf(Dmz)));
%                figure; imagesc(abs(Dmz)); colorbar; title([pftype ' FC(z) matrix']);

                T3 = T2(ids,ids);
                T3(isinf(T3)) = max(T3(~isinf(T3)));
%                figure; imagesc(abs(T3)); colorbar; title([pftype ' T-value matrix']);
                lT3 = log(T3); lT3(lT3<0) = 0;

                % each flys [log10(cell count) vs. FC(z)]
                for i=1:length(CM)
                    Dz = D3z(:,:,i);
                    Dz(isinf(Dz)) = max(Dz(~isinf(Dz)));
                    sbjR(k,i) = corr(lC2(:),abs(Dz(:)));
                end

                % each ROIs [log10(cell count) vs. FC(z)]
                for i=1:size(lC2,1)
                    roiR(k,i) = corr([lC2(i,:)';lC2(:,i)],abs([Dmz(i,:)';Dmz(:,i)]));
                end

                % full ROIs (vs. mean group data)
                R(1) = corr(lC2(:),abs(Dmz(:))); % whole vs. FC(z)
                disp(['prefix=' pftype ' : log10(cell count2) vs. FC(z) = ' num2str(R(1))]);
                R(2) = corr(lW2(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : log10(synapse weight2) vs. FC(z) = ' num2str(R(2))]);
                R(3) = corr(lS(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : log10(synapse count) vs. FC(z) = ' num2str(R(3))]);
                R(4) = corr(W3(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI in-neuron weight vs. FC(z) = ' num2str(R(4))]);
                R(5) = corr(Sw(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI in-synapse weight vs. FC(z) = ' num2str(R(5))]);
                R(6) = corr(Wo(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI out-neuron weight vs. FC(z) = ' num2str(R(6))]);

                R(7) = corr(lC2b(:),abs(Dmz(:))); % Traced vs. FC(z)
                disp(['prefix=' pftype ' : log10(cell count2b) vs. FC(z) = ' num2str(R(7))]);
%                figure; scatter(lC2b(:),abs(Dmz(:))); xlabel('log10(cell count2b)'); ylabel('FC(z)'); title(pftype);
                R(8) = corr(lW2b(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : log10(synapse weight2b) vs. FC(z) = ' num2str(R(8))]);
                R(9) = corr(lSb(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : log10(synapse count b) vs. FC(z) = ' num2str(R(9))]);
                R(10) = corr(W3b(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI in-neuron weight b vs. FC(z) = ' num2str(R(10))]);
                R(11) = corr(Swb(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI in-synapse weight b vs. FC(z) = ' num2str(R(11))]);
                R(12) = corr(Wob(:),abs(Dmz(:)));
                disp(['prefix=' pftype ' : ROI out-neuron weight b vs. FC(z) = ' num2str(R(12))]);

                R(13) = corr(lC2(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(cell count2) vs. FC T-val = ' num2str(R(13))]);
                R(14) = corr(lW2(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapse weight2) vs. FC T-val = ' num2str(R(14))]);
                R(15) = corr(lS(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapse count) vs. FC T-val = ' num2str(R(15))]);
                R(16) = corr(W3(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI in-neuron weight vs. FC T-val = ' num2str(R(16))]);
%                figure; scatter(W3(:),abs(T3(:))); xlabel('ROI in-neuron weight'); ylabel('FC T-val'); title(pftype);
                R(17) = corr(Sw(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI in-synapse weight vs. FC T-val = ' num2str(R(17))]);
                R(18) = corr(Wo(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI out-neuron weight vs. FC T-val = ' num2str(R(18))]);

                R(19) = corr(lC2b(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(cell count2b) vs. FC T-val = ' num2str(R(19))]);
                R(20) = corr(lW2b(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapse weight2b) vs. FC T-val = ' num2str(R(20))]);
                R(21) = corr(lSb(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapse count b) vs. FC T-val = ' num2str(R(21))]);
                R(22) = corr(W3b(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI in-neuron weight b vs. FC T-val = ' num2str(R(22))]);
                R(23) = corr(Swb(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI in-synapse weight b vs. FC T-val = ' num2str(R(23))]);
                R(24) = corr(Wob(:),abs(T3(:)));
                disp(['prefix=' pftype ' : ROI out-neuron weight b vs. FC T-val = ' num2str(R(24))]);

                Rm = [Rm, R'];

                % calculate AUC
                aucmat = ['results/' pftype '-fcauc.mat'];
                if exist(aucmat,'file')
                    % load beta volumes
                    load(aucmat);
                else
                    thN = 100;
                    aths = cell(thN,1);
%                    for th = [1 50 99] %1:thN
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
                        stb = Sb; stb(stb<sbth) = 0; stb(stb>0) = 1;
                        w2bth = prctile(W2b(W2b>0),th-1);
                        wt2b = W2b; wt2b(wt2b<w2bth) = 0; wt2b(wt2b>0) = 1;
                        w3bth = prctile(W3b(W3b>0),th-1);
                        wt3b = W3b; wt3b(wt3b<w3bth) = 0; wt3b(wt3b>0) = 1;
                        swbth = prctile(Swb(Swb>0),th-1);
                        swtb = Swb; swtb(swtb<swbth) = 0; swtb(swtb>0) = 1;
                        wobth = prctile(Wob(Wob>0),th-1);
                        wotb = Wob; wotb(wotb<wobth) = 0; wotb(wotb>0) = 1;

                        aucs = cell(24,1);
                        [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(Dmz(:)'), ['FC(z) vs. cell count2 th=' num2str(th-1)], false);
%                        figure; imagesc(ct2); colorbar; figure; imagesc(abs(Dmz)); colorbar; % to check matrix
                        aucs{1} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt2(:)', abs(Dmz(:)'), ['FC(z) vs. synapse weight2 th=' num2str(th-1)], false);
                        aucs{2} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(Dmz(:)'), ['FC(z) vs. synapse count th=' num2str(th-1)], false);
                        aucs{3} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt3(:)', abs(Dmz(:)'), ['FC(z) vs. ROI in-neuron weight th=' num2str(th-1)], false);
                        aucs{4} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(swt(:)', abs(Dmz(:)'), ['FC(z) vs. ROI in-synapse weight th=' num2str(th-1)], false);
                        aucs{5} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wot(:)', abs(Dmz(:)'), ['FC(z) vs. ROI out-neuron weight th=' num2str(th-1)], false);
                        aucs{6} = single(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(ct2b(:)', abs(Dmz(:)'), ['FC(z) vs. cell count2b th=' num2str(th-1)], false);
                        aucs{7} = single(auc); % input node is one, so it is not vector.
                        [~, ~, auc] = calcShowGroupROCcurve(wt2b(:)', abs(Dmz(:)'), ['FC(z) vs. synapse weight2b th=' num2str(th-1)], false);
                        aucs{8} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(stb(:)', abs(Dmz(:)'), ['FC(z) vs. synapse count b th=' num2str(th-1)], false);
                        aucs{9} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt3b(:)', abs(Dmz(:)'), ['FC(z) vs. ROI in-neuron weight b th=' num2str(th-1)], false);
                        aucs{10} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(swtb(:)', abs(Dmz(:)'), ['FC(z) vs. ROI in-synapse weight b th=' num2str(th-1)], false);
                        aucs{11} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wotb(:)', abs(Dmz(:)'), ['FC(z) vs. ROI out-neuron weight b th=' num2str(th-1)], false);
                        aucs{12} = single(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(T3(:)'), ['FC T-val vs. cell count2 th=' num2str(th-1)], false);
                        aucs{13} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt2(:)', abs(T3(:)'), ['FC T-val vs. synapse weight2 th=' num2str(th-1)], false);
                        aucs{14} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(T3(:)'), ['FC T-val vs. synapse count th=' num2str(th-1)], false);
                        aucs{15} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt3(:)', abs(T3(:)'), ['FC T-val vs. ROI in-neuron weight th=' num2str(th-1)], false);
                        aucs{16} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(swt(:)', abs(T3(:)'), ['FC T-val vs. ROI in-synapse weight th=' num2str(th-1)], false);
                        aucs{17} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wot(:)', abs(T3(:)'), ['FC T-val vs. ROI out-neuron weight th=' num2str(th-1)], false);
                        aucs{18} = single(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(ct2b(:)', abs(T3(:)'), ['FC T-val vs. cell count2b th=' num2str(th-1)], false);
                        aucs{19} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt2b(:)', abs(T3(:)'), ['FC T-val vs. synapse weight2b th=' num2str(th-1)], false);
                        aucs{20} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(stb(:)', abs(T3(:)'), ['FC T-val vs. synapse count b th=' num2str(th-1)], false);
                        aucs{21} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wt3b(:)', abs(T3(:)'), ['FC T-val vs. ROI in-neuron weight b th=' num2str(th-1)], false);
                        aucs{22} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(swtb(:)', abs(T3(:)'), ['FC T-val vs. ROI in-synapse weight b th=' num2str(th-1)], false);
                        aucs{23} = single(auc);
                        [~, ~, auc] = calcShowGroupROCcurve(wotb(:)', abs(T3(:)'), ['FC T-val vs. ROI out-neuron weight b th=' num2str(th-1)], false);
                        aucs{24} = single(auc);
                        aths{th} = aucs;
                    end
                    A = zeros(24,thN,'single');
                    for th = 1:thN
                        for j=1:24
                            A(j,th) = aths{th}{j}; % input node is one, so it is not vector.
                        end
                    end
                    save(aucmat,'A','R');
                end
                AUC = cat(3,AUC,A);

%                figure; plot(A'); title(['prefix=' pftype ' : FC AUC result vs. ground truth SC']);
%                xlabel('percentile');
            end
        end
    end
    % FC-SC correlation (6-type mixed box plot)
    figure; boxplot(Rm(7:12,:),'Labels',rlabel); title([roitype ' FC-SC correlation (6-type mixed plot)']);
    hold on; plot(Rm(7:12,:)'); hold off; legend;

    % FC-SC detection (FC(z) vs. cell count)
    % cell count: only ROI was transformed. cell count2: synapse points were transformed and re-counted in all ROI. 
%{
    A1 = squeeze(AUC(1,:,:)); % for internal check
    figure; boxplot(A1,'Labels',rlabel); title([roitype ' FC-SC detection (FC(z) vs. cell count)']);
    A5 = squeeze(AUC(5,:,:));
    figure; boxplot(A5,'Labels',rlabel); title([roitype ' FC-SC detection (FC(z) vs. cell count2)']);

    A2 = squeeze(AUC(2,:,:)); % for internal check
    figure; boxplot(A2,'Labels',rlabel); title([roitype ' FC-SC detection (FC(z) vs. synapse weight)']);
    A6 = squeeze(AUC(6,:,:));
    figure; boxplot(A6,'Labels',rlabel); title([roitype ' FC-SC detection (FC(z) vs. synapse weight2)']);
%}
    AA = squeeze(nanmean(AUC,2));
    figure; boxplot(AA(7:12,:),'Labels',rlabel); title([roitype ' FC-SC detection (6-type mixed plot)']);
    hold on; plot(AA(7:12,:)'); hold off; legend;
end
