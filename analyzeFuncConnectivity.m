% analyze functional connectivity matrix (1st, 2nd level analysis).
% need to run makeStructConnectivity.m, extractROItimeseries.m first.

function analyzeFuncConnectivity
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
%    hpfTh = [0]; % high-pass filter threshold
    hpfTh = [0, 0.1, 0.05, 0.025, 0.02, 0.01, 0.009, 0.008, 0.005, 0.001]; % high-pass filter threshold
%    prefix = {'', 'hf', 's10', 's10hf', 's20', 's20hf', 's30', 's30hf', 's40', 's40hf', ...
%         's50', 's50hf', 's60', 's50hf'};
%    prefix = {'m10'};
    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
%    nuisance = {'', 'poltcomp'}; % good for flyemroi
    nuisance = {'', 'tacomp'}; % good for bransonhemi

    % using subjects (flys). sbj 7 shows NaN row in FC matrix
    sbjids = [1 2 3 4 5 6 8 9];

    % ROI name
%    roitype = 'flyemroi';
%    roitype = 'bransonhemi';
    roitype = 'hemicube4';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load structural connectivity matrix (from makeStructConnectivity.m)
    switch(roitype)
    case 'flyemroi'
        load('data/neuprint_connectlist.mat');
    
        ids = primaryIds;
%        ids = [20	111	59	68	51	62	106	87	24	27	75	50	54]; % MB only
        ids = [103	107	20	111	59	68	65	78  49	51	62	106	87	47 100 24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	80	1	102	63	95	56];
    case 'bransonhemi'
        load('data/branson_connectlist.mat');
        ids = primaryIds;
    case 'hemicube4'
    end
    C = countMat(ids,ids); W = weightMat(ids,ids);

    % show corr between cell count v. synapse weight
    corr(W(:),C(:))    % corr between cell count v. synapse weight
    figure; scatter(W(:),C(:)); xlabel('synapse weight'); ylabel('cell count')

    n = length(ids);
    E = eye(n); E = logical(1-E);

    lC = log10(C); lC(lC<0) = 0;
    lW = log10(W); lW(lW<0) = 0;
%    lC = C;

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

                % transport first
                outfile = ['results/' pftype '-func.mat'];
                if exist(outfile,'file')
                    load(outfile);
                else
                    % load ROI time-series (from extractROItimeseries.m)
                    load(['results/' pftype '-ts.mat']);
        
                    CM = {};
                    for i=1:length(sbjids)
                        CX{sbjids(i)} = CX{sbjids(i)}';
%                        CX{i} = highpass(CX{i},hpfTh,1/TR);
                        CM{i} = corr(CX{sbjids(i)}); %calcPartialCorrelation_(CX{i}',[],[],[],false,1e-3);
                    end
                    [B2, RSS2, T2, df] = calcSeedCorrMixed(CX(sbjids));
                    
                    % output beta matrix
                    save(outfile,'df','B2','RSS2','T2','CM','-v7.3');
                end

                D3 = []; D3z = [];
                for i=1:length(CM)
                    D = CM{i}(ids,ids);
                    D3 = cat(3,D3,D);
                end
                D3z = atanh(D3); % z transformed (better FC-SC corr).
    %           D3z = D3; % no transform

                % mean group data
                Dm = nanmean(D3,3);
                Dmz = nanmean(D3z,3);
    %            Dmz = atanh(Dmz);
                Dmz(isinf(Dmz)) = max(Dmz(~isinf(Dmz)));

                T3 = T2(ids,ids);
                T3(isinf(T3)) = max(T3(~isinf(T3)));
                figure; imagesc(abs(T3)); colorbar; title([pftype ' ROI FC matrix']);
                lT3 = log(T3); lT3(lT3<0) = 0;

                % each flys [log10(cell count) vs. FC(z)]
                for i=1:length(CM)
                    Dz = D3z(:,:,i);
                    Dz(isinf(Dz)) = max(Dz(~isinf(Dz)));
                    sbjR(k,i) = corr(lC(:),abs(Dz(:)));
                end

                % each ROIs [log10(cell count) vs. FC(z)]
                for i=1:size(lC,1)
                    roiR(k,i) = corr([lC(i,:)';lC(:,i)],abs([Dmz(i,:)';Dmz(:,i)]));
                end

                % full ROIs (vs. mean group data)
                R(1) = corr(lC(:),abs(Dmz(:)));
%                figure; scatter(lC(:),abs(Dmz(:))); xlabel('log10(cell count)'); ylabel('FC(z)'); title(pftype);
                disp(['prefix=' pftype ' : log10(cell count) vs. FC(z) = ' num2str(R(1))]);

                R(2) = corr(lW(:),abs(Dmz(:)));
%                figure; scatter(lW(:),abs(Dmz(:))); xlabel('log10(synapse weight'); ylabel('FC(z)'); title(pftype);
                disp(['prefix=' pftype ' : log10(synapse weight) vs. FC(z) = ' num2str(R(2))]);

                R(3) = corr(lC(:),abs(T3(:)));
%                figure; scatter(lC(:),abs(T3(:))); xlabel('log10(cell count)'); ylabel('FC T-val'); title(pftype);
                disp(['prefix=' pftype ' : log10(cell count) vs. FC T-val = ' num2str(R(3))]);

                R(4) = corr(lW(:),abs(T3(:)));
                disp(['prefix=' pftype ' : log10(synapse weight) vs. FC T-val = ' num2str(R(4))]);
                Rm = [Rm, R'];

                % calculate AUC
                aucmat = ['results/' pftype '-fcauc.mat'];
                if exist(aucmat,'file')
                    % load beta volumes
                    load(aucmat);
                else
                    thN = 100;
                    A = zeros(4,thN);
                    As = cell(4,thN);
                    for th = 1:thN
                        % include injection voxel in ground truth
                        cth = prctile(C(C>0),th-1);
                        ct = C; ct(ct<cth) = 0; ct(ct>0) = 1;
                        wth = prctile(W(W>0),th-1);
                        wt = W; wt(wt<wth) = 0; wt(wt>0) = 1;

                        [~, ~, auc] = calcShowGroupROCcurve(ct(:)', abs(Dmz(:)'), ['FC(z) vs. cell count th=' num2str(th-1)], false);
                        As{1,th} = auc;
                        A(1,th) = nanmean(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(wt(:)', abs(Dmz(:)'), ['FC(z) vs. synapse weight th=' num2str(th-1)], false);
                        As{2,th} = auc;
                        A(2,th) = nanmean(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(ct(:)', abs(T3(:)'), ['FC T-val vs. cell count th=' num2str(th-1)], false);
                        As{3,th} = auc;
                        A(3,th) = nanmean(auc);

                        [~, ~, auc] = calcShowGroupROCcurve(wt(:)', abs(T3(:)'), ['FC T-val vs. synapse weight th=' num2str(th-1)], false);
                        As{4,th} = auc;
                        A(4,th) = nanmean(auc);
                    end
                    save(aucmat,'A','As');
                end
                AUC = cat(3,AUC,A);

%                figure; plot(A'); title(['prefix=' pftype ' : FC AUC result vs. ground truth SC']);
%                xlabel('percentile');
            end
        end
    end
    % FC-SC correlation (4-type mixed box plot)
    figure; boxplot(Rm,'Labels',rlabel); title('FC-SC correlation (4-type mixed plot)');

    % FC-SC detection (FC(z) vs. cell count)
    A1 = squeeze(AUC(1,:,:));
    figure; boxplot(A1,'Labels',rlabel); title('FC-SC detection (FC(z) vs. cell count)');

    A2 = squeeze(AUC(2,:,:));
    figure; boxplot(A2,'Labels',rlabel); title('FC-SC detection (FC(z) vs. synapse weight)');

    AA = reshape(AUC,[4*size(AUC,2),size(AUC,3)]);
    figure; boxplot(A2,'Labels',rlabel); title('FC-SC detection (4-type mixed plot)');
end
