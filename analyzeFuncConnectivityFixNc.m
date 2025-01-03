% re-calculation FC-SC of post-synapse data.
% this script should run after analyzeFuncConnectivity.m.

function analyzeFuncConnectivityFixNc
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = [0]; % high-pass filter threshold
%    hpfTh = [0, 0.1, 0.05, 0.025, 0.02, 0.01, 0.009, 0.008, 0.005, 0.001]; % high-pass filter threshold
%    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
%    smooth = {'s90', 's100', 's110', 's120', 's130', 's140', 's150', 's160', 's170', 's180', 's190', 's200', 's210', 's220', 's230', 's240', 's250', 's260', 's270', 's280', 's290', 's300'};
%    smooth = {'', 's30', 's80'};
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
%    nuisance = {'6hmtacomp'}; % good for bransonhemi, branson7065km50
%    nuisance = {'','poltcomp'}; % good for DistKm(synapse)
    nuisance = {''};

    % using subjects (flys). sbj 7 shows NaN row in FC matrix
    sbjids = [1 2 3 4 5 6 8 9];

    % ROI name
%    roitypes = {'hemiroif'};  % flyem ROI full
%    roitypes = {'hemiroi','bransonhemi'}; % flyem ROI (Turner compatible)
%    roitypes = {'hemiroi_hb0sr80','hemiroi_fw0','hemiroi_avg0'}; % flyem ROI (Primary, FlyEM vs. FlyWire vs. Average)
%    roitypes = {'hemiBranson7065'};
%    roitypes = {'hemiBranson7065km20','hemiBranson7065km30','hemiBranson7065km50','hemiBranson7065km100','hemiBranson7065km200', ...
%        'hemiBranson7065km300','hemiBranson7065km500','hemiBranson7065km300'};
%    roitypes = {'hemiCmkm20','hemiCmkm30','hemiCmkm50','hemiCmkm100','hemiCmkm200', ...
%        'hemiCmkm300','hemiCmkm500','hemiCmkm1000'};
%    roitypes = {'hemiCmkm20r1w1','hemiCmkm30r1w1','hemiCmkm50r1w1','hemiCmkm100r1w1','hemiCmkm200r1w1', ...
%        'hemiCmkm300r1w1','hemiCmkm500r1w1','hemiCmkm1000r1w1'};
%    roitypes = {'hemiDistKm20','hemiDistKm30','hemiDistKm50','hemiDistKm100','hemiDistKm200', ...
%        'hemiDistKm300','hemiDistKm500','hemiDistKm1000'};
%    roitypes = {'hemiRand20','hemiRand30','hemiRand50','hemiRand100','hemiRand200', ...
%        'hemiRand300','hemiRand500','hemiRand1000'};
%    roitypes = {'hemiVrand20','hemiVrand30','hemiVrand50','hemiVrand100','hemiVrand200', ...
%        'hemiVrand300','hemiVrand500','hemiVrand1000'};
%    roitypes = {'hemiCube12','hemiCube8','hemiCube4'};
%    roitypes = {'hemiPiece12','hemiPiece8','hemiPiece4'};
%    roitypes = {'hemiPiece3','hemiPiece2'};
    % neuropil FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
%    roitypes = {'hemiRoi101','hemiRoi57','hemiRoi57-51','hemiRoi51-62-20-111-100'};
%    roitypes = {'hemiRoi1','hemiRoi5','hemiRoi7','hemiRoi27','hemiRoi30','hemiRoi32','hemiRoi43','hemiRoi52', ...
%        'hemiRoi54','hemiRoi57','hemiRoi59','hemiRoi63','hemiRoi65','hemiRoi67','hemiRoi78','hemiRoi82', ...
%        'hemiRoi89','hemiRoi93','hemiRoi95','hemiRoi100','hemiRoi101','hemiRoi106','hemiRoi113'};
    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiBranson7065km50','hemiBranson7065km50_fw0sr50','hemiCmkm50','hemiCmkm50_fw0sr50', ...
        'hemiCmkm50r1w1','hemiDistKm50','hemiDistKm50_fw0sr50','hemiRand50','hemiVrand50'};
%    roitypes = {'hemiBranson7065km30','hemiCmkm30','hemiCmkm30r1w1','hemiDistKm30','hemiRand30','hemiVrand30'};
%    roitypes = {'hemiCmkm20000','hemiCmkm20000r1w1','hemiDistKm20000','hemiVrand20000'};
%    roitypes = {'hemiCmkm100','hemiDistKm100','hemiCmkm500','hemiDistKm500','hemiCmkm1000','hemiDistKm1000'};
%    roitypes = {'hemiCmkm5000','hemiDistKm5000','hemiCmkm10000','hemiDistKm10000'};
    % FlyEM vs. FlyWire in each roi num by smoothing 0 to 80 (no nuisanse)
    roitypes = {'hemiBranson7065km20_fw0sr50','hemiBranson7065km30_fw0sr50','hemiBranson7065km50_fw0sr50','hemiBranson7065km100_fw0sr50','hemiBranson7065km200_fw0sr50','hemiBranson7065km300_fw0sr50','hemiBranson7065km500_fw0sr50','hemiBranson7065km1000_fw0sr50', ...
        'hemiCmkm20_fw0sr50','hemiCmkm30_fw0sr50','hemiCmkm50_fw0sr50','hemiCmkm100_fw0sr50','hemiCmkm200_fw0sr50','hemiCmkm300_fw0sr50','hemiCmkm500_fw0sr50', 'hemiCmkm1000_fw0sr50',...
        'hemiDistKm20_fw0sr50','hemiDistKm30_fw0sr50','hemiDistKm50_fw0sr50','hemiDistKm100_fw0sr50','hemiDistKm200_fw0sr50','hemiDistKm300_fw0sr50','hemiDistKm500_fw0sr50','hemiDistKm1000_fw0sr50'};
%    roitypes = {'hemiCmkm50','hemiDistKm50','hemiCmkm100','hemiDistKm100','hemiCmkm500','hemiDistKm500'}; % for large smoothing size & no nuisanse, poltcomp
%    roitypes = {'hemiroi','hemiroi_fw0sr50','hemiDistKm50','hemiDistKm50_fw0sr50','hemiDistKm50_avg'}; % for all nuisanse & s30, s80 % for s0 to s80 (no nuisanse)
%    roitypes = {'hemiroi_hb0sr50','hemiroi_hb0sr60','hemiroi_hb0sr70','hemiroi_hb0sr80','hemiroi_hb0sr90', ... % for s30 & s80, '' & poltcomp
%            'hemiroi_fw0sr50','hemiroi_fw0sr70','hemiroi_fw0sr100','hemiroi_fw0sr130','hemiroi_fw0sr140','hemiroi_fw0sr150'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:length(roitypes)
        analyzeFcROItype(roitypes{n}, preproc, hpfTh, smooth, nuisance, sbjids)
    end
end

function analyzeFcROItype(roitype, preproc, hpfTh, smooth, nuisance, sbjids)

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
    if ~exist('ncountMat','var'), disp(['cannot find ncountMat in ' roitype ' ... skip']); return; end

    C2 = ncountMat(ids,ids,1); S = sycountMat(ids,ids,1);
    C2b = ncountMat(ids,ids,2); Sb = sycountMat(ids,ids,2);

    lC2 = log10(C2); lC2(lC2<0) = 0;
    lS = log10(S); lS(lS<0) = 0;
    lC2b = log10(C2b); lC2b(lC2b<0) = 0;
    lSb = log10(Sb); lSb(lSb<0) = 0;

    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                pftype = [smooth{k} hpfstr nuisance{n} preproc roitype];

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
                load(outfile);

                F3 = [];
                for i=1:cmlen
                    F = CM{i}(ids,ids);
                    F3 = cat(3,F3,F);
                end
                F3z = atanh(F3); % z transformed (better FC-SC corr).
                clear CM;

                % mean group data
                mFz = nanmean(F3z,3);
                mFz(isinf(mFz)) = max(mFz(~isinf(mFz)));

                T3 = T2(ids,ids);
                T3(isinf(T3)) = max(T3(~isinf(T3)));

                aucmat = ['results/auc/' pftype '-fcauc.mat'];
                if ~exist(aucmat,'file'), disp(['cannot find ' aucmat]); end
                load(aucmat);

                % full ROIs (vs. mean group data)
                R(1) = corr(lC2(:),abs(mFz(:))); % whole vs. m-FCz
                disp(['prefix=' pftype ' : log10(neurons2) vs. m-FCz = ' num2str(R(1))]);
                R(3) = corr(lS(:),abs(mFz(:)));
                disp(['prefix=' pftype ' : log10(synapses) vs. m-FCz = ' num2str(R(3))]);
                R(7) = corr(lC2b(:),abs(mFz(:))); % Traced vs. m-FCz
                disp(['prefix=' pftype ' : log10(neurons2b) vs. m-FCz = ' num2str(R(7))]);
                R(9) = corr(lSb(:),abs(mFz(:)));
                disp(['prefix=' pftype ' : log10(synapses b) vs. m-FCz = ' num2str(R(9))]);
                R(13) = corr(lC2(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(neurons2) vs. FC-Tval = ' num2str(R(13))]);
                R(15) = corr(lS(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapses) vs. FC-Tval = ' num2str(R(15))]);
                R(19) = corr(lC2b(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(neurons2b) vs. FC-Tval = ' num2str(R(19))]);
                R(21) = corr(lSb(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapses b) vs. FC-Tval = ' num2str(R(21))]);

                thN = 100;
                aths = cell(thN,1);
                for th = 1:thN
                    for j=1:24
                        if ~isempty(A(j,th))
                            aths{th}{j} = A(j,th); % input node is one, so it is not vector.
                        end
                    end
                end
%                for th = 1:thN
                parfor th = 1:thN
                    % include injection voxel in ground truth
                    c2th = prctile(C2(C2>0),th-1);
                    ct2 = C2; ct2(ct2<c2th) = 0; ct2(ct2>0) = 1;
                    sth = prctile(S(S>0),th-1);
                    st = S; st(st<sth) = 0; st(st>0) = 1;

                    c2bth = prctile(C2b(C2b>0),th-1);
                    ct2b = C2b; ct2b(ct2b<c2bth) = 0; ct2b(ct2b>0) = 1;
                    sbth = prctile(Sb(Sb>0),th-1);
                    stb = Sb; stb(stb<sbth) = 0; stb(stb>0) = 1;

                    aucs = aths{th};
                    [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(mFz(:)'), ['m-FCz vs. neurons2 th=' num2str(th-1)], false);
                    aucs{1} = single(auc);
                    [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(mFz(:)'), ['m-FCz vs. synapses th=' num2str(th-1)], false);
                    aucs{3} = single(auc);
                    [X, Y, auc] = calcShowGroupROCcurve(ct2b(:)', abs(mFz(:)'), ['m-FCz vs. neurons2b th=' num2str(th-1)], false);
                    aucs{7} = single(auc); % input node is one, so it is not vector.
                    [~, ~, auc] = calcShowGroupROCcurve(stb(:)', abs(mFz(:)'), ['m-FCz vs. synapses b th=' num2str(th-1)], false);
                    aucs{9} = single(auc);
                    [~, ~, auc] = calcShowGroupROCcurve(ct2(:)', abs(T3(:)'), ['FC-Tval vs. neurons2 th=' num2str(th-1)], false);
                    aucs{13} = single(auc);
                    [~, ~, auc] = calcShowGroupROCcurve(st(:)', abs(T3(:)'), ['FC-Tval vs. synapses th=' num2str(th-1)], false);
                    aucs{15} = single(auc);
                    [~, ~, auc] = calcShowGroupROCcurve(ct2b(:)', abs(T3(:)'), ['FC-Tval vs. neurons2b th=' num2str(th-1)], false);
                    aucs{19} = single(auc);
                    [~, ~, auc] = calcShowGroupROCcurve(stb(:)', abs(T3(:)'), ['FC-Tval vs. synapses b th=' num2str(th-1)], false);
                    aucs{21} = single(auc);
                    aths{th} = aucs;
                end

                A = nan(24,thN,'single');
                for th = 1:thN
                    for j=1:24
                        if ~isempty(aths{th}{j})
                            A(j,th) = aths{th}{j}; % input node is one, so it is not vector.
                        end
                    end
                end
                save(aucmat,'A','R');
            end
        end
    end
end
