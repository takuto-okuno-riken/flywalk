% make structural connectivity matrix from neuprint connectivity list (csv).
% countMat: only ROI was transformed. countMat2: synapse points were transformed and re-counted in all ROI. 
% this script should run after makeVoxelROIatlas.m

function makeStructConnectivity
    SCVER = 5;
    hbSth = 80;     % FlyEM hemibrain synapse threshold
    fwSth = 140;    % FlyWire synapse score threshold
    synTh = 0;      % synapse count at one neuron threshold
    epsilon = 3000; % dbscan range (nanometer)
    minpts = 1;     % set 1, includes isolated synapse

    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs.
    % list data was acquired by c.fetch_roi_connectivity() of neuprint python api.
%%{
    % primary ROIs
    primaryIds = [103	107	20	111	59	68	65	78	34	4	49	51	62	106	87	47	100	24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	74	42	80	1	102	63	95	56];
    roiNum = 114;
    turnerIds = [103	107	20	111	59	68	65	78  49	51	62	106	87	47 100 24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	80	1	102	63	95	56];

    % load matfile
    idstr = 'hemiroi';
    fname = ['results/sc/' lower(idstr) '_connectlist.mat']; scver = 1;
    load(fname);

    %
    if ~exist('countMat','var')
        countMat = zeros(roiNum,roiNum,'single'); % neurons matrix
        weightMat = zeros(roiNum,roiNum,'single'); % synapse weight matrix
        for i=1:size(connectlist,1)
            countMat(connectlist(i,1),connectlist(i,2)) = connectlist(i,3);
            weightMat(connectlist(i,1),connectlist(i,2)) = connectlist(i,4);
        end
    end
    if ~exist('syweightMat','var')
        [roiIdxs, sz, aV] = getHemiroiRoiIdxs(turnerIds);
%        niftiwrite(aV,'atlas/hemiFlyem52atlasCal.nii','Compressed',true);

        [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));

        save([fname(1:end-4) '_cnids.mat'],'Cnids','-v7.3');
        scver = 3;
    end

    % check large sparse version
%{
        roiIdxs = {};
        listing = dir(['atlas/hemiroi/*.nii.gz']);
        for i=1:length(listing)
            V = niftiread(['atlas/hemiroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            roiIdxs{i} = find(V>0);
            sz = size(V);
        end
        [countMat3, sycountMat3] = makeSCcountMatrixLarge(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
        dcm = ncountMat(:,:,2) - full(countMat3);
        disp(['ncountMat-countMat3 : '  num2str(sum(abs(dcm),'all'))]);
        dsm = sycountMat(:,:,2) - full(sycountMat3);
        disp(['sycountMat2-sycountMat3 : '  num2str(sum(abs(dsm),'all'))]);

        [weightMat3] = makeSCweightMatrixLarge(sycountMat3, 'hemiroi');
        swm = single(full(weightMat3) ./ sycountMat(:,:,2)); swm(isnan(swm)) = 0;
        dwm = nweightMat(:,:,2) - sycountMat(:,:,2) .* swm;
        disp(['nweightMat-weightMat3 : '  num2str(sum(abs(dwm),'all'))]);
%}

    % find primary ROI ids from list 
    % !caution! neurite does not return full primary ROI connections
%    id1 = unique(connectlist(:,1));
%    id2 = unique(connectlist(:,2));
%    id3 = unique([id1, id2]);
%{
    ids = primaryIds; CM = countMat(ids,ids); WM = weightMat(ids,ids); 
    figure; imagesc(log(CM)); colorbar; title('hemibrain neurons matrix');
    figure; imagesc(log(WM)); colorbar; title('hemibrain all synapse weight matrix');
%}
    % primary, R/L, name order
    load('data/hemiroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagescLabel(log10(CM), labelNames, 'hemibrain neurons matrix');
    figure; imagescLabel(log10(SM), labelNames, 'hemibrain all synapses matrix');
    WMi = nweightMat(ids,ids,2);
    WMo = outweightMat(ids,ids,2);
    figure; imagescLabel(log(WMi'), labelNames, 'hemibrain all input-rate matrix');
    figure; imagescLabel(log(WMo), labelNames, 'hemibrain all output-rate matrix');
%{
    [roiIdxs, sz, aV] = getHemiroiRoiIdxs(ids);
    niftiwrite(aV,'atlas/hemiFlyemPrimaryatlasCal.nii','Compressed',true);
%}
    % reorder by neuron count matrix clustering
    % for Turner et al. (2021) compatible (around 50 ROIs)
    if scver <= SCVER
        % reorder by tree clustering
        ids = turnerIds;
        CM = ncountMat(ids,ids,2);
        eucD = pdist(CM,'euclidean');
        Z = linkage(eucD,'ward');
        [H,T,outperm] = dendrogram(Z,103);
        primaryIds = turnerIds(outperm); % use default leaf order
        scver = 5;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'connectlist','countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end
    ids = primaryIds;
    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log10(CM)); colorbar; daspect([1 1 1]); title('hemibrain neurons matrix');
    figure; imagesc(log10(SM)); colorbar; daspect([1 1 1]); title('hemibrain synapses matrix');
%    CM2b = ncountMat(ids,ids,2); SMb = sycountMat(ids,ids,2); WM2b = nweightMat(ids,ids,2); WMob = outweightMat(ids,ids,2);
%    WM3b = WM2b ./ SMb; WM3b(SMb==0) = 0; % pure ROI-input neuron connection weight
%}
    % ---------------------------------------------------------------------
    % calc ROI distance matrix based on distance based k-means atlas.
    %
%{
    fname = ['results/sc/' lower(idstr) '_dist.mat'];

    clfname = ['results/sc/' lower(idstr) '_connectlist.mat'];
    load(clfname);

    if exist(fname,'file')
        load(fname);
    else
        listing = dir(['atlas/hemiroi/*.nii.gz']);
        roimax = length(listing);
        P1 = zeros(roimax,1,3);
        P2 = zeros(1,roimax,3);
        for i=1:roimax
            V = niftiread(['atlas/hemiroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            sz = size(V);
            idx = find(V>0);
            X = zeros(length(idx),3);
            for j=1:length(idx)
                [x,y,z] = ind2sub(sz,idx(j));
                X(j,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
            end
            P1(i,1,:) = mean(X,1);
            P2(1,i,:) = mean(X,1);
        end
        P1 = repmat(P1,[1 roimax 1]);
        P2 = repmat(P2,[roimax 1 1]);
        distMat = single(sqrt(sum((P1 - P2).^2,3)));
        save(fname,'distMat','primaryIds','roiNum','scver','-v7.3');
    end

    ids = primaryIds;
    CM = ncountMat(ids,ids,2); DM = distMat(ids,ids);
    figure; imagesc(DM); colorbar; title([idstr ' distance matrix']);
    r = corr(CM(:),DM(:));
    figure; scatter(CM(:),DM(:)); title([idstr ' sc vs. dist r=' num2str(r)]); xlabel('neuron SC matrix'); ylabel('distance matrix');
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs by FlyWire EM data.
    %
%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth)];
    fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
    if exist(fname,'file')
        load(fname);
    else
        [roiIdxs, sz] = getHemiroiRoiIdxs();

        conf = getSCconfig('wire', synTh, fwSth);
        [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));

        countMat = []; weightMat = []; scver = 5;
    end
    checkScverAndSave(idstr, fname, ordertemp, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs
    % cPPSSI is applied for threshold
%{
    functype = 'fw'; %''; %
%{
    rnum = 0;
    rtype = 0; % for spiTh=90, logi inversion ~(0<=spiTh<=90)
%    rtype = 3; % for spiTh=10, logi as is (0<=spiTh<=10) '_only'
    mdstr = '_md'; % for morphological-based ''; % for straight-line distance
    ii = 1;
    for k=90%[10 20 40 60 80 90]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
%{
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
%}
%{
    rtype = 4;
    k = 0;
    ii = 1;
    for rnum = [5e5 6e5 7e5 13e5 15e5 5e6 6e6 7e6 7.1e6 10e6]
%}
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) '_sp' num2str(k) 'db' num2str(epsilon) 'mi' num2str(minpts) mdstr];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();

            if isempty(functype)
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), k*100, epsilon, minpts, 0, rtype, rnum, {2, 2});
            else
                conf = getSCconfig('hemi', synTh, hbSth);
%                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), k*100, epsilon, minpts, 0, rtype, rnum);
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), k*100, epsilon, minpts, 0, rtype, rnum, [], mdstr);
                nweightMat = []; outweightMat = []; syweightMat = [];
            end
            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_hb0sr80_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs by FlyWire EM data.
    % cPPSSI is applied for threshold
%{
%{
    rnum = 0;
    rtype = 0; % for spiTh=90, logi inversion ~(0<=spiTh<=90)
%    rtype = 3; % for spiTh=10, logi as is (0<=spiTh<=10) '_only'
    ii = 1;
    mdstr = '_md'; % for morphological-based ''; % for straight-line distance
    for k=90%[10 20 40 60 80 90]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
%{
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
%}
%{
    rtype = 4;
    k = 0;
    ii = 1;
    for rnum = [500000 1500000 5000000 10000000]
%}
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_sp' num2str(k) 'db' num2str(epsilon) 'mi' num2str(minpts) mdstr];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();

            conf = getSCconfig('wire', synTh, fwSth);
%            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), k*100, epsilon, minpts, 0, rtype, rnum, [], mdstr);
            [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), k*100, epsilon, minpts, 0, rtype, rnum, [], mdstr);
            nweightMat = []; outweightMat = []; syweightMat = [];

            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_hb0sr80_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs
    % full random subsampling is not good. Setting sparsity constraints.
%%{
    functype = 'fw'; %''; %
%    randrange = {[6.5e5, 1.5e5],[14e5, 1.5e5],[70.5e5, 4e5]}; % old. old.
%    randrange = {[32e5, 2e5, 0.4233, 0.015, 40],[42.5e5, 3.5e5, 0, 0, 150]}; % (new) morphological-based distance
    randrange = {[32e5, 0, 0.4233, 0, 40],[42.5e5, 0, 0, 0, 150],[14e5, 0, 0.3850, 0, 100]}; % (new) for check histogram. morphological-based distance
    for ii=2%1:length(randrange)
        param = randrange{ii};
        for k=1:99
            clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;

            idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rd' num2str(param(1)/10000) '-' num2str(param(2)/10000)  '-' num2str(param(5)) '-' num2str(k)];
            fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
            if exist(fname,'file')
                load(fname);
            else
                rnum = int32(param(1) + randn() * param(2)); % subsampling number
                spnum = param(3) + randn() * param(4); % sparcity rate
                [roiIdxs, sz] = getHemiroiRoiIdxs();

                if spnum > 0
                    spmat = getSparseRandMatix(['results/sc/hemiroi_hb0sr80' functype '_connectlist.mat'], spnum, length(roiIdxs));
                else
                    spmat = [];
                end

                if isempty(functype)
                    [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), 0, epsilon, minpts, 0, 4, rnum, {2, 2});
                else
                    conf = getSCconfig('hemi', synTh, hbSth);
                    spsubsamp = getSparseSubsampSynFw(conf, ['hemiroi_hb0sr80' functype], spmat, param(5));
%                    [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, 0, 4, rnum);
                    [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, 0, 6, rnum, spsubsamp);
                    nweightMat = []; outweightMat = []; syweightMat = [];
                end
                countMat = []; weightMat = []; scver = 5;
            end
            % not 52 ROIs, but this one is 63 ROIs
            checkScverAndSave(idstr, fname, ['results/sc/hemiroi_hb0sr80' functype '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
        end
    end
%}
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs by FlyWire EM data.
    % full random subsampling is not good. Setting sparsity constraints.
%{
%    randrange = {[14e5, 1.5e5, 0, 0, 0],[18.5e5, 2e5, 0, 0, 0],[155.5e5, 8e5, 0, 0, 0]}; % old. old.
%    randrange = {[19.5e5, 2e5, 0.4863, 0.015, 60],[32e5, 2e5, 0.4175, 0.015, 35]}; % old. straight-line distance
%    randrange = {[19.5e5, 0, 0.4863, 0, 60],[32e5, 0, 0.4175, 0.015, 35]}; % old. for check histogram. straight-line distance
    randrange = {[17e5, 1.5e5, 0.5203, 0.015, 35],[220e5, 8e5, 0, 0, 150]}; % (new) morphological-based distance
%    randrange = {[17e5, 0, 0.5203, 0, 35],[220e5, 0, 0, 0, 150]}; % (new) for check histogram. morphological-based distance
    for ii=1:length(randrange)
        param = randrange{ii};
        for k=1:99
            clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;

            idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rd' num2str(param(1)/10000) '-' num2str(param(2)/10000)  '-' num2str(param(5)) '-' num2str(k)];
            fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
            if exist(fname,'file')
                load(fname);
            else
                rnum = int32(param(1) + randn() * param(2)); % subsampling number
                spnum = param(3) + randn() * param(4); % sparcity rate
                [roiIdxs, sz] = getHemiroiRoiIdxs();

                if spnum > 0
                    spmat = getSparseRandMatix('results/sc/hemiroi_fw0sr140_connectlist.mat', spnum, length(roiIdxs));
                else
                    spmat = [];
                end
                conf = getSCconfig('wire', synTh, fwSth);
                spsubsamp = getSparseSubsampSynFw(conf, 'hemiroi_fw0sr140', spmat, param(5));

%                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, 0, 6, rnum, spsubsamp);
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, 0, 6, rnum, spsubsamp);
                nweightMat = []; outweightMat = []; syweightMat = [];

                countMat = []; weightMat = []; scver = 5;
            end
            % not 52 ROIs, but this one is 63 ROIs
            checkScverAndSave(idstr, fname, 'results/sc/hemiroi_fw0sr140_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
        end
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs
    % reciprocal synapse distance is applied for threshold
%{
    functype = ''; %'fw'; %
%{
    rnum = 0;
    rtype = 3;
    ii = 1;
    for k=[20 40 100 500 1000 10000]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rc' num2str(k)];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();

            if isempty(functype)
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum, {2, 2});
            else
                conf = getSCconfig('hemi', synTh, hbSth);
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum);
            end

            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_hb0sr80_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs by FlyWire EM data.
    % reciprocal synapse distance is applied for threshold
%{
%{
    rnum = 0;
    rtype = 3;
    ii = 1;
    for k=[20 40 100 500 1000 10000]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rc' num2str(k)];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum);
    
            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_fw0sr140_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs
    % separation index & reciprocal synapse distance are applied for threshold
%{
    SR = {[5 40],[5 10000]};
    for k=1:length(SR)
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        idstr = ['hemiroi_hb'  num2str(synTh) 'sr' num2str(hbSth) '_sp' num2str(SR{k}(1)) 'db' num2str(epsilon) 'mi' num2str(minpts) '_rc' num2str(SR{k}(2))];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();
    
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), SR{k}(1)*100, epsilon, minpts, SR{k}(2)*100);

            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_hb0sr80_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of flyem hemibrain neuropil ROIs by FlyWire EM data.
    % separation index & reciprocal synapse distance are applied for threshold
%{
    SR = {[5 40],[5 10000]};
    for k=1:length(SR)
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        idstr = ['hemiroi_fw'  num2str(synTh) 'sr' num2str(fwSth) '_sp' num2str(SR{k}(1)) 'db' num2str(epsilon) 'mi' num2str(minpts) '_rc' num2str(SR{k}(2))];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz] = getHemiroiRoiIdxs();

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), SR{k}(1)*100, epsilon, minpts, SR{k}(2)*100);
    
            countMat = []; weightMat = []; scver = 5;
        end
        % not 52 ROIs, but this one is 63 ROIs
        checkScverAndSave(idstr, fname, 'results/sc/hemiroi_hb0sr80_connectlist.mat', countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson matrix csv.
    % extract ROI ids from hemibrain mask
    % branson matrix csv was acquired from Turner et al. (2021).
%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    idstr = 'bransonhemi';
    fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
    if exist(fname,'file')
        load(fname);
    else
        % read hemibrain mask
        Vm = niftiread('template/jrc2018f_flyemhemibrainCal_invFDACal.nii.gz');
        Vm(Vm>0) = 1;
        Vm(Vm<1) = 0;
    
        % read branson atlas (FDA registered)
        Vb = niftiread('atlas/JRC2018_branson_atlasCal_invFDACal.nii.gz');
        % check full region is included
        primaryIds = [];
        Vbm = Vb .* Vm;
        for i=1:999
            len1 = length(find(Vb==i));
            len2 = length(find(Vbm==i));
            if len1==len2
                primaryIds = [primaryIds, i];
            end
        end
%        primaryIds = unique(Vb(Vm>0)); primaryIds(primaryIds==0) = [];
        roiNum = length(primaryIds);
    
        countMat = readmatrix('data/JRC2018_branson_cellcount_matrix.csv');
        countMat = countMat(2:end,2:end); % remove headers
        weightMat = readmatrix('data/JRC2018_branson_weighted_tbar_matrix.csv');
        weightMat = weightMat(2:end,2:end); % remove headers

        % primary id should have neurons
        c = sum(countMat,2);
        for i=1:length(primaryIds)
            if c(primaryIds(i)) == 0
                disp(['primaryId=' num2str(primaryIds(i)) ' does not have neurons']);
            end
        end
        save(fname,'countMat','weightMat','primaryIds','roiNum');
    end
    if ~exist('syweightMat','var')
        roiIdxs = {};
        % read branson atlas (FDA registered)
        Vb = niftiread('atlas/JRC2018_branson_atlasCal_invFDACal.nii.gz');
        sz = size(Vb);
        for i=1:999
            if ismember(i,primaryIds)
                roiIdxs{i} = find(Vb==i);
            else
                roiIdxs{i} = [];
            end
        end

        [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
        scver = 5;
    end
    checkScverAndSave(idstr, fname, [], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    %
%{
    idstr = 'hemibranson7065';
    fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    if exist(fname,'file')
        load(fname);
    else
        [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs('atlas/hemiBranson7065atlasCal.nii.gz');

        [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));

        countMat = []; weightMat = []; scver = 5;
    end
    checkScverAndSave(idstr, fname, [], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    %
%{
    idstr = 'wirebranson7065';
    fname = ['results/sc/' idstr '_connectlist.mat'];
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    if exist(fname,'file')
        load(fname);
    else
        [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs('atlas/hemiBranson7065atlasCal.nii.gz');

        conf = getSCconfig('wire', synTh, fwSth);
        [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));

        countMat = []; weightMat = []; scver = 5;
    end
    checkScverAndSave(idstr, fname, [], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas.
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas by FlyWire EM data
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k) '_fw'  num2str(synTh) 'sr' num2str(fwSth)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiBranson7065km' num2str(k) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas (flyWire based clustering).
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['wireBranson7065km' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas (flyWire based clustering) by FlyWire EM data
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['wireBranson7065km' num2str(k) '_fw'  num2str(synTh) 'sr' num2str(fwSth)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/wireBranson7065km' num2str(k) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for cube ROI.
    % extract ROI ids from hemicube4 mask
%{
    for k=4%[12 8 4]
        idstr = ['hemiCube' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' lower(idstr) 'atlasCal.nii.gz']);

            [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
    
            countMat = []; weightMat = []; syweightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for piece ROI.
    % extract ROI ids from hemipiece3 mask
%{
    for k=[12 8 4 3 2]
        idstr = ['hemiPiece' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' lower(idstr) 'atlasCal.nii.gz']);

            [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
    
            countMat = []; weightMat = []; scver = 5;
        end
        checkScverAndSave(idstr, fname, [], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix for neuropil specific voxels.
    % extract voxel ids from fanshape-body (FB), eliptic-body (EB), and other atlas.
%{
    roiids = {[68 59 87 106 50 27 54]}; % a'L(R)-aL(R)-b'L(R)-bL(R)-gL(R)-CA(R)-PED(R)
%    roiids = {[101],[57],[57,51],[51,62,20,111,100]}; % FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
%    roiids = {1	5	7	27	30	32	43	52	54	57	59	63	65	67	78	82	89	93	95	100	101	106	113};
    for k = 1:length(roiids)
        idstr = ['hemiRoi' num2str(roiids{k}(1))];
        for j=2:length(roiids{k}), idstr=[idstr '-' num2str(roiids{k}(j))]; end

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);
    
            [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
    
            countMat = []; weightMat = []; syweightMat = [];  scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix for neuropil specific voxels & distance k-means ROI atlas
    % extract voxel ids from mushroom body atlas.
%{
    roiids = {[68 59 87 106 50 27 54]}; % a'L(R)-aL(R)-b'L(R)-bL(R)-gL(R)-CA(R)-PED(R)
    for i = 1:length(roiids)
        roiname = ['hemiRoi' num2str(roiids{i}(1))];
        for j=2:length(roiids{i}), roiname=[roiname '-' num2str(roiids{i}(j))]; end

        for k=[200 400 800 1600 3200]
            idstr = [roiname 'DistKm' num2str(k)];
            fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

            clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
            if exist(fname,'file')
                load(fname);
            else
                [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);
        
                [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
        
                countMat = []; weightMat = []; syweightMat = [];  scver = 5;
            end
            checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
        end
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for all EM ROI voxels (except fibers).
%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    idstr = 'hemiRoiWhole';
    fname = ['results/sc/' idstr '_connectlist.mat'];
    if exist([fname(1:end-4) '_cm.mat'],'file')
        load([fname(1:end-4) '_cm.mat']);
        load([fname(1:end-4) '_sm.mat']);
    else
        [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

        [countMat, sycountMat] = makeSCcountMatrixLarge(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
        
        save([fname(1:end-4) '_cm.mat'],'countMat','-v7.3');
        save([fname(1:end-4) '_sm.mat'],'sycountMat','-v7.3');
    end
    if exist([fname(1:end-4) '_wm.mat'],'file')
        load([fname(1:end-4) '_wm.mat']);
    else
        roiNum = size(sycountMat,1);
        primaryIds = 1:roiNum;
        [weightMat] = makeSCweightMatrixLarge(sycountMat, lower(idstr));

        save([fname(1:end-4) '_wm.mat'],'weightMat','primaryIds','roiNum','scver','-v7.3');
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means atlas.
    %
%%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000]
        idstr = ['hemiCmkm' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means (smoothing) atlas.
    %
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000]
        idstr = ['hemiCmkm' num2str(k) 'r1w1'];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means atlas by FlyWire EM data
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiCmkm' num2str(k) '_fw'  num2str(synTh) 'sr' num2str(fwSth)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiCmkm' num2str(k) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas.
    %
%%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % calc ROI distance matrix based on distance based k-means atlas.
    %
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_dist.mat'];

        clfname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        load(clfname);

        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/' idstr 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);

            P1 = zeros(roimax,1,3);
            P2 = zeros(1,roimax,3);
            for i=1:roimax
                idx = find(atlV==i);
                X = zeros(length(idx),3);
                for j=1:length(idx)
                    [x,y,z] = ind2sub(sz,idx(j));
                    X(j,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
                end
                P1(i,1,:) = mean(X,1);
                P2(1,i,:) = mean(X,1);
            end
            P1 = repmat(P1,[1 roimax 1]);
            P2 = repmat(P2,[roimax 1 1]);
            distMat = single(sqrt(sum((P1 - P2).^2,3)));
            save(fname,'distMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; daspect([1 1 1]); title([idstr ' neurons matrix']);
        figure; imagesc(distMat(ids,ids)); colorbar; daspect([1 1 1]); title([idstr ' distance matrix']);
        r = corr(CM(:),distMat(:));
        figure; scatter(CM(:),distMat(:)); title([idstr ' sc vs. dist r=' num2str(r)]); xlabel('neuron SC matrix'); ylabel('distance matrix');
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means n voxels atlas.
    %
%{
    for k=[1000]
        for n=[1 2 4 8 16 32 64 128]
            idstr = ['hemiDistKm' num2str(k) 'vox' num2str(n)];
            fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
    
            clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
            if exist(fname,'file')
                load(fname);
            else
                [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);
    
                [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
    
                countMat = []; weightMat = []; syweightMat = []; scver = 5;
            end
            checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
        end
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas by FlyWire EM data.
    %
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000]
        idstr = ['hemiDistKm' num2str(k) '_fw' num2str(synTh) 'sr' num2str(fwSth)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiDistKm' num2str(k) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas
    % pre-post synapse separation index is applied for threshold
%{
    roisz = 500;
    for k=[5 10 15 20 30]
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        idstr = ['hemiDistKm' num2str(roisz) '_hb'  num2str(synTh) 'sr' num2str(hbSth) '_sp' num2str(k) 'db' num2str(epsilon) 'mi' num2str(minpts)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiDistKm' num2str(roisz) 'atlasCal.nii.gz']);
    
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), k*100, epsilon, minpts);

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas by FlyWire EM data.
    % pre-post synapse separation index is applied for threshold
%{
    roisz = 500;
    for k=[5 10 15 20 30]
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        idstr = ['hemiDistKm' num2str(roisz) '_fw'  num2str(synTh) 'sr' num2str(fwSth) '_sp' num2str(k) 'db' num2str(epsilon) 'mi' num2str(minpts)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiDistKm' num2str(roisz) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), k*100, epsilon, minpts);
    
            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas
    % reciprocal synapse distance is applied for threshold
%{
    roisz = 500;
    functype = ''; %'fw'; %
%{
    rnum = 0;
    rtype = 3;
    ii = 1;
    for k=[20 40 100 500 1000 10000]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
%{
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
%}
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiDistKm' num2str(roisz) '_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiDistKm' num2str(roisz) '_hb'  num2str(synTh) 'sr' num2str(hbSth) functype '_rc' num2str(k)];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiDistKm' num2str(roisz) 'atlasCal.nii.gz']);

            if isempty(functype)
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum, {2, 2});
            else
                conf = getSCconfig('hemi', synTh, hbSth);
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum);
            end

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas by FlyWire EM data.
    % reciprocal synapse distance is applied for threshold
%{
    roisz = 500;
%{
    rnum = 0;
    rtype = 3;
    ii = 1;
    for k=[20 40 100 500 1000 10000]
%}
%{
    rnum = 0;
    rtype = 5;
    k = 20;
    for ii=1:3
%}
%{
    rnum = 100000;
    rtype = 4;
    k = 0;
    for ii=1:3
%}
        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if rnum > 0
            idstr = ['hemiDistKm' num2str(roisz) '_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rn' num2str(rnum/10000)];
        else
            idstr = ['hemiDistKm' num2str(roisz) '_fw'  num2str(synTh) 'sr' num2str(fwSth) '_rc' num2str(k)];
        end
        switch(rtype)
        case 1,  idstr = [idstr '_rand' num2str(ii)];
        case 2,  idstr = [idstr '_xrand' num2str(ii)];
        case 3,  idstr = [idstr '_only' num2str(ii)];
        case 4,  idstr = [idstr '_orand' num2str(ii)];
        case 5,  idstr = [idstr '_xorand' num2str(ii)];
        end
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/hemiDistKm' num2str(roisz) 'atlasCal.nii.gz']);

            conf = getSCconfig('wire', synTh, fwSth);
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, conf, lower(idstr), 0, epsilon, minpts, k*100, rtype, rnum);

            countMat = []; weightMat = []; scver = 5;
        end
        str = split(idstr,'_');
        checkScverAndSave(idstr, fname, ['results/sc/' lower(str{1}) '_connectlist.mat'], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of hemibrain primary ROI atlas (flyem hemibrain) by taking mean from FlyEM and FlyWire.
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_avg_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            t1name = ['results/sc/' lower(idstr) '_connectlist.mat'];
            t1 = load(t1name);
            t2name = ['results/sc/' lower(idstr) '_fw'  num2str(synTh) 'sr' num2str(fwSth) '_connectlist.mat'];
            t2 = load(t2name);
            ncountMat = (t1.ncountMat + t2.ncountMat) / 2;
            nweightMat = (t1.nweightMat + t2.nweightMat) / 2;
            sycountMat = (t1.sycountMat + t2.sycountMat) / 2;
            syweightMat = (t1.syweightMat + t2.syweightMat) / 2;
            outweightMat = (t1.outweightMat + t2.outweightMat) / 2;
            primaryIds = t1.primaryIds;
            roiNum = t1.roiNum;

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverAndSave(idstr, fname, [], countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from fully random cluster atlas.
    %
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiRand' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from fully random voxel atlas.
    %
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000]
        idstr = ['hemiVrand' num2str(k)];
        fname = ['results/sc/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            [roiIdxs, sz, primaryIds, roiNum] = getRoiIdxs(['atlas/' idstr 'atlasCal.nii.gz']);

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, hbSth/100, synTh, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver);
    end
%}
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

function [roiIdxs, sz, aV] = getHemiroiRoiIdxs(ids)
    if nargin < 1, ids = []; end
    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    aV = niftiread(info); aV(:) = 0;
    cnt = 1;
    roiIdxs = {};
    listing = dir(['atlas/hemiroi/*.nii.gz']);
    for i=1:length(listing)
        V = niftiread(['atlas/hemiroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
        idx = find(V>0);
        roiIdxs{i} = idx;
        if any(ismember(ids,i))
            aV(idx) = cnt; cnt = cnt + 1;
        end
        sz = size(V);
    end
end

function spmat = getSparseRandMatix(ordertemp, spnum, roiNum)
    spmat = zeros(roiNum,roiNum,'logical');
    load(ordertemp);
    ids = primaryIds;
    spmat(ids,ids) = 1;
    idx = find(spmat);
    permidx = randperm(length(idx));
    spmat(idx(permidx(ceil(length(idx)*spnum):end))) = 0;
end

function spsubsamp = getSparseSubsampSynFw(conf, type, spmat, preTh)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;

    % read neuron info (id, type)
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)
    % read synapse info
    Sid = []; postNidx = []; preNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find
    if isempty(spmat)
        spsubsamp = valid & score;
        return;
    end

    load(['results/cache/' type '_Nin_Nout.mat']); % should exist
    roimax = size(spmat,1);
    % remove output neurons which belongs specific ROIs
    outnlogi = zeros(length(Nid),1,'int8');
    for i=1:roimax
        outnidx = Nout{i}{2}; % post synapse neuron
        for j=1:roimax
            if spmat(i,j) == 1, continue; end
            innidx = Nin{j}{2}; % pre synapse neuron
            % find input neuron from ROI(i)
            nlogi = ismember(innidx,outnidx);
            idx = innidx(nlogi);
            outnlogi(idx) = outnlogi(idx) + 1;
        end
    end
    nidx = find(outnlogi >= preTh);
    slogi1 = ismember(preNidx,nidx);
%    slogi1 = ismember(postNidx,nidx);
    spsubsamp = ~slogi1 & valid & score;
end

function checkScverAndSave(idstr, fname, ordertemp, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver)
    SCVER = 5;
    if scver <= SCVER && ~isempty(ordertemp)
        % set same order based on template file.
        cl = load(ordertemp);
        primaryIds = cl.primaryIds;
        clear cl;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end
%{
    ids = primaryIds; CM = ncountMat(ids,ids,1); WM = sycountMat(ids,ids,1);
    figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix full']);
    figure; imagesc(log(WM)); colorbar; title([idstr ' synapses matrix full']);
%}
    ids = primaryIds; CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
    figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
end

function checkScverTreeAndSave(idstr, fname, countMat, weightMat, ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, primaryIds, roiNum, scver)
    SCVER = 5;
    if scver <= SCVER
        % reorder by tree clustering
        cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
        eucD = pdist(cm2,'euclidean');
        Z = linkage(eucD,'ward');
        [H,T,outperm] = dendrogram(Z,roiNum);
        primaryIds = outperm; % use default leaf order
%        dendrogram(Z,200); % for figure.2
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end

    ids = primaryIds;
    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; daspect([1 1 1]); title([idstr ' neurons matrix']);
    figure; imagesc(log(SM)); colorbar; daspect([1 1 1]); title([idstr ' synapses matrix']);
end

function imagescLabel(mat, labelNames, titlestr)
    imagesc(mat); colorbar; daspect([1 1 1]); title(titlestr);
    set(gca,'XTick',1:size(mat,1));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end
