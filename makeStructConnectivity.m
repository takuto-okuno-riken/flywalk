% make structural connectivity matrix from neuprint connectivity list (csv).
% countMat: only ROI was transformed. countMat2: synapse points were transformed and re-counted in all ROI. 
% this script should run after makeVoxelROIatlas.m

function makeStructConnectivity
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from neuprint connectivity list.
    % list data was acquired by c.fetch_roi_connectivity() of neuprint python api.
%%{
    % primary ROIs
    primaryIds = [103	107	20	111	59	68	65	78	34	4	49	51	62	106	87	47	100	24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	74	42	80	1	102	63	95	56];
    roiNum = 114;
    turnerIds = [103	107	20	111	59	68	65	78  49	51	62	106	87	47 100 24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	80	1	102	63	95	56];

    % load matfile
    fname = 'data/neuprint_connectlist.mat';
    load(fname);

    %
    if ~exist('countMat','var')
        countMat = zeros(roiNum,roiNum,'single'); % cell count matrix
        weightMat = zeros(roiNum,roiNum,'single'); % synapse weight matrix
        for i=1:size(connectlist,1)
            countMat(connectlist(i,1),connectlist(i,2)) = connectlist(i,3);
            weightMat(connectlist(i,1),connectlist(i,2)) = connectlist(i,4);
        end
        save(fname,'connectlist','countMat','weightMat','primaryIds','roiNum');
    end

    if ~exist('syweightMat','var')
        info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
        aV = niftiread(info); aV(:) = 0;
        cnt = 1;
        roiIdxs = {};
        listing = dir(['atlas/flyemroi/*.nii.gz']);
        for i=1:length(listing)
            V = niftiread(['atlas/flyemroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            idx = find(V>0);
            roiIdxs{i} = idx;
            if any(ismember(turnerIds,i))
                aV(idx) = cnt; cnt = cnt + 1;
            end
            sz = size(V);
        end
%        niftiwrite(aV,'atlas/hemiFlyem52atlasCal.nii','Compressed',true);

        [countMat2, sycountMat, weightMat2, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'hemiroi');

        save(fname,'connectlist','countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum');
        save([fname(1:end-4) '_cnids.mat'],'Cnids','-v7.3');
    end

    % check large sparse version
%{
        roiIdxs = {};
        listing = dir(['atlas/flyemroi/*.nii.gz']);
        for i=1:length(listing)
            V = niftiread(['atlas/flyemroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            roiIdxs{i} = find(V>0);
            sz = size(V);
        end
        [countMat3, sycountMat3] = makeSCcountMatrixLarge(roiIdxs, sz, 0.8, 0, 'hemiroi');
        dcm = countMat2(:,:,2) - full(countMat3);
        disp(['countMat2-countMat3 : '  num2str(sum(abs(dcm),'all'))]);
        dsm = sycountMat(:,:,2) - full(sycountMat3);
        disp(['sycountMat2-sycountMat3 : '  num2str(sum(abs(dsm),'all'))]);

        [weightMat3] = makeSCweightMatrixLarge(sycountMat3, 'hemiroi');
        swm = single(full(weightMat3) ./ sycountMat(:,:,2)); swm(isnan(swm)) = 0;
        dwm = weightMat2(:,:,2) - sycountMat(:,:,2) .* swm;
        disp(['weightMat2-weightMat3 : '  num2str(sum(abs(dwm),'all'))]);
%}

    % find primary ROI ids from list 
    % !caution! neurite does not return full primary ROI connections
%    id1 = unique(connectlist(:,1));
%    id2 = unique(connectlist(:,2));
%    id3 = unique([id1, id2]);
%{
    ids = primaryIds; CM = countMat(ids,ids); WM = weightMat(ids,ids); 
    figure; imagesc(log(CM)); colorbar; title('hemibrain cell count matrix');
    figure; imagesc(log(WM)); colorbar; title('hemibrain all synapse weight matrix');
%}
    % primary, R/L, name order
    load('data/flyemroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    CM2 = countMat2(ids,ids,2); WM2 = weightMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagescLabel(log10(CM2), labelNames, 'hemibrain cell count matrix');
    figure; imagescLabel(log10(SM), labelNames, 'hemibrain all synapse count matrix');
    figure; imagescLabel(log10(WM2), labelNames, 'hemibrain all synapse weight matrix');
    WMi = WM2 ./ SM; WMi(SM==0) = 0; % pure ROI-input neuron connection weight
    WMo = outweightMat(ids,ids,2);
    figure; imagescLabel(log(WMi'), labelNames, 'hemibrain all input-rate matrix');
    figure; imagescLabel(log(WMo), labelNames, 'hemibrain all output-rate matrix');

    % reorder by neuron count matrix clustering
    % for Turner et al. (2021) compatible (around 50 ROIs)
    ids = turnerIds;
    CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
    eucD = pdist(CM2,'euclidean');
    Z = linkage(eucD,'ward');
    ids = optimalleaforder(Z,eucD);
    figure; imagesc(log10(CM2(ids,ids))); colorbar; title('hemibrain cell count matrix');
    figure; imagesc(log10(SM(ids,ids))); colorbar; title('hemibrain synapse count matrix');
    
%    CM2b = countMat2(ids,ids,2); SMb = sycountMat(ids,ids,2); WM2b = weightMat2(ids,ids,2); WMob = outweightMat(ids,ids,2);
%    WM3b = WM2b ./ SMb; WM3b(SMb==0) = 0; % pure ROI-input neuron connection weight
%}

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson matrix csv.
    % extract ROI ids from hemibrain mask
    % branson matrix csv was acquired from Turner et al. (2021).
%{
    clear countMat2; clear sycountMat; clear weightMat2;
    fname = 'data/branson_connectlist.mat';
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

        % primary id should have cell count
        c = sum(countMat,2);
        for i=1:length(primaryIds)
            if c(primaryIds(i)) == 0
                disp(['primaryId=' num2str(primaryIds(i)) ' does not have cell count']);
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

        [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson');
        save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum');
    end
%{
    ids = primaryIds; CM = countMat(ids,ids); WM = weightMat(ids,ids);
    figure; imagesc(log(CM)); colorbar; title('branson cell count matrix');
    figure; imagesc(log(WM)); colorbar; title('branson synapse weight matrix');
%}
    ids = primaryIds; CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM2)); colorbar; title('branson cell count 2 matrix');
    figure; imagesc(log(SM)); colorbar; title('branson synapse count matrix');

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    % extract ROI ids from hemibrain mask

    fname = 'data/hemibranson7065_connectlist.mat';
    clear countMat2; clear sycountMat; clear weightMat2;
    if exist(fname,'file')
        load(fname);
    else
        atlV = niftiread('atlas/hemiBranson7065atlasCal.nii.gz');
        roimax = max(atlV(:));
        sz = size(atlV);

        roiIdxs = {};
        for i=1:roimax
            roiIdxs{i} = find(atlV==i);
        end

        primaryIds = 1:roimax;
        roiNum = length(primaryIds);

        [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson7065');

        countMat = []; weightMat = [];
        save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum');
    end

    ids = primaryIds;
    CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM2)); colorbar; title(['branson7065 cell count 2 matrix']);
    figure; imagesc(log(SM)); colorbar; title(['branson7065 synapse count matrix']);

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas.
    % extract ROI ids from hemibrain mask
%}
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiBranson7065km' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        if primaryIds(1) == 1 && primaryIds(roiNum) == roiNum
            % reorder by tree clustering
            cm2 = countMat2(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            primaryIds = optimalleaforder(Z,eucD);
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end
%}
%{
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for cube ROI.
    % extract ROI ids from hemicube4 mask
    for k=[12 8 4]
        idstr = ['hemiCube' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/' idstr 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','primaryIds','roiNum');
        end
    
        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for piece ROI.
    % extract ROI ids from hemipiece3 mask
    for k=[12 8 4 3 2]
        idstr = ['hemiPiece' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/' idstr 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','primaryIds','roiNum','-v7.3');
        end
    
        ids = primaryIds;
        CM2 = countMat2(ids,ids); SM = sycountMat(ids,ids);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for fanshape body (FB) and others.
    % extract voxel ids from fanshape-body (FB), eliptic-body (EB), and other atlas.
%    roiids = {[101],[57],[57,51],[51,62,20,111,100]}; % FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
    roiids = {1	5	7	27	30	32	43	52	54	57	59	63	65	67	78	82	89	93	95	100	101	106	113};

    for k = 1:length(roiids)
        idstr = num2str(roiids{k}(1));
        for j=2:length(roiids{k}), idstr=[idstr '-' num2str(roiids{k}(j))]; end

        clear countMat2; clear sycountMat; clear weightMat2;
        fname = ['data/hemiroi' idstr '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiRoi' idstr 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, ['hemiroi' idstr]);
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','primaryIds','roiNum');
        end
    
        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title(['hemiroi ' idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title(['hemiroi ' idstr ' synapse count matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for all EM ROI voxels (except fibers).
%{
    clear countMat2; clear sycountMat; clear weightMat2;
    fname = ['data/hemiroiwhole_connectlist.mat'];
    if exist([fname(1:end-4) '_cm.mat'],'file')
        load([fname(1:end-4) '_cm.mat']);
        load([fname(1:end-4) '_sm.mat']);
    else
        atlV = niftiread(['atlas/hemiRoiWholeatlasCal.nii.gz']);
        roimax = max(atlV(:));
        sz = size(atlV);

        roiIdxs = {};
        for i=1:roimax
            roiIdxs{i} = find(atlV==i);
        end

        primaryIds = 1:roimax;
        roiNum = length(primaryIds);

        [countMat, sycountMat] = makeSCcountMatrixLarge(roiIdxs, sz, 0.8, 0, 'hemiRoiWhole');
        
        save([fname(1:end-4) '_cm.mat'],'countMat','-v7.3');
        save([fname(1:end-4) '_sm.mat'],'sycountMat','-v7.3');
    end
    if exist([fname(1:end-4) '_wm.mat'],'file')
        load([fname(1:end-4) '_wm.mat']);
    else
        roiNum = size(sycountMat,1);
        primaryIds = 1:roiNum;
        [weightMat] = makeSCweightMatrixLarge(sycountMat, 'hemiRoiWhole');

        save([fname(1:end-4) '_wm.mat'],'weightMat','primaryIds','roiNum','-v7.3');
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means atlas.
    % extract ROI ids from hemibrain mask

    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiCmkm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiCmkm' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        if primaryIds(1) == 1 && primaryIds(roiNum) == roiNum
            % reorder by tree clustering
            cm2 = countMat2(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            primaryIds = optimalleaforder(Z,eucD);
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means (smoothing) atlas.
    % extract ROI ids from hemibrain mask

    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiCmkm' num2str(k) 'r1w1'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiCmkm' num2str(k) 'r1w1atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        if primaryIds(1) == 1 && primaryIds(roiNum) == roiNum
            % reorder by tree clustering
            cm2 = countMat2(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            primaryIds = optimalleaforder(Z,eucD);
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas.
    % extract ROI ids from hemibrain mask

    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiDistKm' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        if primaryIds(1) == 1 && primaryIds(roiNum) == roiNum
            % reorder by tree clustering
            cm2 = countMat2(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            primaryIds = optimalleaforder(Z,eucD);
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapse count matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from fully random cluster atlas.
    % extract ROI ids from hemibrain mask

    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiRand' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear sycountMat; clear weightMat2;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiRand' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
    
            countMat = []; weightMat = [];
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        if primaryIds(1) == 1 && primaryIds(roiNum) == roiNum
            % reorder by tree clustering
            cm2 = countMat2(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            primaryIds = optimalleaforder(Z,eucD);
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','syweightMat','primaryIds','roiNum','-v7.3');
        end

        ids = primaryIds;
        CM2 = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM2),[8 10]); colorbar; title([idstr ' cell count 2 matrix']);
        figure; imagesc(log(SM),[8 12]); colorbar; title([idstr ' synapse count matrix']);
    end
end

function imagescLabel(mat, labelNames, titlestr)
    imagesc(mat); colorbar; daspect([1 1 1]); title(titlestr);
    set(gca,'XTick',1:size(mat,1));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end

function [countMat, sycountMat, weightMat, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrix(roiIdxs, sz, rateTh, synTh, type)

    % read neuron info (id, connection number, size)
    Nid = []; Nstatus = [];
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % read synapse info
    Sdir = []; StoN = []; Srate = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear StoS; clear Sloc;
    Sid = uint32(1:length(StoN));

    % read synapse location in FDA
    SlocFc = [];
    load('data/synapseloc_fdacal.mat');

    % make presynapse index
    cfile = 'results/hemibrain_v1_2_synapseCell.mat';
    if exist(cfile,'file')
        load(cfile);
    else
        C = cell(sz(1),sz(2),sz(3));
        for i=1:size(SlocFc,1)
            t = ceil(SlocFc(i,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                C{t(1),t(2),t(3)} = [C{t(1),t(2),t(3)},i];
            else
                disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
            end
        end
        save(cfile,'C','-v7.3');
    end

    % use only accurate synapse more than 'rate'
    idx = find(Srate < rateTh); % use only accurate synapse more than 'rate'
    Sdir(idx) = 0; 
    clear Srate;

    isoutsyweight = (nargout >= 5);
    roimax = length(roiIdxs);
    nfile = ['results/cache-' type '_Nin_Nout.mat'];
    if exist(nfile,'file')
        load(nfile);
    else
        Nin = cell(roimax,1); Nout = cell(roimax,1);
        Sin = cell(roimax,1); Sout = cell(roimax,1);
        for i=1:roimax
            if isempty(roiIdxs{i}), continue; end
            disp(['Nin_Nout ' num2str(i) ' / ' num2str(roimax)]);
    
            % find post synapses in that Cube ROI.
    %        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
    %        C{x,y,z}
            D = C(roiIdxs{i});
            sids = [];
            for j=1:length(D)
                sids = [sids, D{j}];
            end
            sids = unique(sids);
            idx = find(Sdir(sids)==2); % get post-synapse ids in this ROI
            postsids = sids(idx);
            outnids = unique(StoN(postsids)); % get (post-synapse) traced & orphan body-ids
            % get pre and post neurons in ROI(i)
            Nout{i}{1} = outnids; % all output cells (including orphan, etc)
            logis = ismember(Nid,outnids);
            Nout{i}{2} = Nid(logis & Nstatus==1); % output traced neuron in ROI(i)
            logis = ismember(outnids, Nout{i}{2});
            Nout{i}{3} = outnids(~logis); % output orphan bodys
    
            idx = find(Sdir(sids)==1); % get pre-synapse ids in this ROI
            presids = sids(idx);
            innids = unique(StoN(presids)); % get (pre-synapse) traced & orphan body-ids
            Nin{i}{1} = innids; % all input cells (including orphan, etc)
            logis = ismember(Nid,innids);
            Nin{i}{2} = Nid(logis & Nstatus==1); % input traced neuron in ROI(i)
            logis = ismember(innids, Nin{i}{2});
            Nin{i}{3} = innids(~logis); % input orphan bodys
%{
            % get pre and post synapses in ROI(i)
            Sout{i}{1} = postsids; % all output synapses (including orphan, etc)
            logis = ismember(StoN,Nout{i}{2});
            logis = ismember(postsids,Sid(logis)); % find may be slow, but Sid will consume memory.
            Sout{i}{2} = postsids(logis); % output traced synapses in ROI(i)
            Sout{i}{3} = postsids(~logis); % output orphan bodys' synapses
%}
            % syweightMat is heavy. calculate is if only it is required.
            if isoutsyweight
                Sin{i}{1} = presids; % all input synapses (including orphan, etc)
                logis = ismember(StoN,Nin{i}{2});
                logis = ismember(presids,Sid(logis)); % find may be slow, but Sid will consume memory.
                Sin{i}{2} = presids(logis); % output traced synapses in ROI(i)
                Sin{i}{3} = presids(~logis); % output orphan bodys' synapses
            end
        end
        if isoutsyweight
            save(nfile,'Nin','Nout','Sin','-v7.3');
        else
            save(nfile,'Nin','Nout','-v7.3');
        end
    end
    clear Nid; clear Nstatus;

    roimax = length(roiIdxs);
    sycountMat = nan(roimax,roimax,3,'single');
    countMat = nan(roimax,roimax,3,'single');
    weightMat = nan(roimax,roimax,3,'single');
    outweightMat = nan(roimax,roimax,3,'single');
    if isoutsyweight
        syweightMat = nan(roimax,roimax,3,'single');
    end

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(32);

    CC = cell(roimax,1);
%    for i=1:roimax
    parfor i=1:roimax
        if isempty(roiIdxs{i}), continue; end
        disp(['get synaptic connections of ROI ' num2str(i) ' / ' num2str(roimax)]);

        % find post synapses in that Cube ROI.
%        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
%        C{x,y,z}

        % three patterns, full (including orphan, etc), neurons, others (orphan, etc)
        CX = cell(roimax,3);
        for p=1:3
            outnids = Nout{i}{p};

            % ROI(i) output all cells to pre-synapses for other ROIs
            logi = ismember(StoN,outnids); % find synapses which belong to ROI(i) output neurons
            outsids = Sid(logi);
            idx = find(Sdir(outsids)==1); % get pre-synapse ids of output neurons
            presids = outsids(idx);
    
            % get connected synapse counts in each ROI (from ROI to connected ROI)
            Vs = zeros(sz(1),sz(2),sz(3),'int32');
            conSlocFc = SlocFc(presids,:); % get 3D location in FDA Cal template.
            N = cell(sz(1),sz(2),sz(3));
            for j=1:size(conSlocFc,1)
                t = ceil(conSlocFc(j,:));
                if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                    Vs(t(1),t(2),t(3)) = Vs(t(1),t(2),t(3)) + 1;
                    N{t(1),t(2),t(3)} = [N{t(1),t(2),t(3)},StoN(presids(j))];
                else
                    disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
                end
            end
            X = zeros(1,roimax,'int32');
            for j=1:roimax
                if isempty(roiIdxs{j}), continue; end
                B = Vs(roiIdxs{j});
                X(j) = nansum(B,1);
            end
            sycountMat(i,:,p) = X;

            % get connected neuron counts in each ROI (from ROI to connected ROI)
            X = zeros(1,roimax,'int32');
            for j=1:roimax
                D = N(roiIdxs{j});
                nids = [];
                for k=1:length(D)
                    nids = [nids, D{k}];
                end
                numsyn = groupcounts(nids'); % number of synapse in each neuron
                nids = unique(nids);
                CX{j,p} = nids(numsyn >= synTh);
                X(j) = length(CX{j,p}); % synapse number threshold for neurons in ROI(j)
            end
            countMat(i,:,p) = X;
        end
        CC{i} = CX;
    end
    clear SlocFc; clear Sdir;

    % calculate weight matrix (full, neurons, others)
    for p=1:3
%        for i=1:roimax
        parfor i=1:roimax
            if isempty(Nout{i}), continue; end
            disp(['ROI ' num2str(i) ' / ' num2str(roimax)]);
            outnids = Nout{i}{p};
            X = zeros(1,roimax,'single');
            Y = zeros(1,roimax,'single');
            SX = zeros(1,roimax,'single');
            for j=1:roimax
                if isempty(Nin{j}), continue; end
                innids = Nin{j}{p};
                insids = Sin{j}{p};
                % find input neuron rate from ROI(i)
                logi = ismember(innids,outnids);
                X(j) = single(sum(logi)) / length(innids); % in-weight (from i to j)
                % syweightMat is heavy. calculate is if only it is required.
                if isoutsyweight
                    logis = ismember(StoN,innids(logi));
                    logis = ismember(insids,Sid(logis));
                    SX(j) = single(sum(logis)) / length(insids); % in-synaptic-weight (from i to j)
                end
                % find output neuron rate from ROI(i)
                logi2 = ismember(outnids,innids); % actually, same as x(j)
                Y(j) = single(sum(logi2)) / length(outnids); % out-weight (from i to j)
                if sum(logi) ~= sum(logi2)
                    disp(['error logi: ' num2str(i) '-' num2str(j)]); % error check. this should not show.
                end
            end
            weightMat(i,:,p) = X;
            outweightMat(i,:,p) = Y;
            if isoutsyweight
                syweightMat(i,:,p) = SX;
            end
        end
    end
    weightMat = weightMat .* sycountMat;

    % pure input and output cells (full, neurons, others count)
    if nargout >= 6
        Ncount = zeros(roimax,2,3,'single');
        for p=1:3
            for i=1:roimax
                if ~isempty(Nin{i})
                    Ncount(i,1,p) = length(Nin{i}{p});
                end
                if ~isempty(Nout{i})
                    Ncount(i,2,p) = length(Nout{i}{p});
                end
            end
        end
    end

    % connected neuron ids (cell), only output neurons and others, but not full.
    if nargout >= 7
        Cnids = cell(roimax,roimax,3);
        for p=2:3
            for i=1:roimax
                for j=1:roimax
                    if ~isempty(CC{i})
                        Cnids{i,j,p} = CC{i}{j,p};
                    end
                end
            end
        end
    end
end

function [countMat, sycountMat] = makeSCcountMatrixLarge(roiIdxs, sz, rateTh, synTh, type)

    % read neuron info (id, connection number, size)
    Nid = []; Nstatus = [];
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize;

    % read synapse info
    Sdir = []; StoN = []; Srate = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear StoS; clear Sloc;
    Sid = uint32(1:length(StoN));

    % read synapse location in FDA
    SlocFc = [];
    load('data/synapseloc_fdacal.mat');

    % make presynapse index
    cfile = 'results/hemibrain_v1_2_synapseCell.mat';
    if exist(cfile,'file')
        load(cfile);
    else
        C = cell(sz(1),sz(2),sz(3));
        for i=1:size(SlocFc,1)
            t = ceil(SlocFc(i,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                C{t(1),t(2),t(3)} = [C{t(1),t(2),t(3)},i];
            else
                disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
            end
        end
        save(cfile,'C','-v7.3');
    end

    % use only accurate synapse more than 'rate'
    idx = find(Srate < rateTh); % use only accurate synapse more than 'rate'
    Sdir(idx) = 0; 
    clear Srate;

    roimax = length(roiIdxs);
    nfile = ['results/cache-' type '_L_Nin_Nout.mat'];
    if exist(nfile,'file')
        load(nfile);
    else
        Nin = cell(roimax,1); Nout = cell(roimax,1);
        for i=1:roimax
            if isempty(roiIdxs{i}), continue; end
            disp(['Nin_Nout ' num2str(i) ' / ' num2str(roimax)]);
    
            % find post synapses in that Cube ROI.
    %        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
    %        C{x,y,z}
            D = C(roiIdxs{i});
            sids = [];
            for j=1:length(D)
                sids = [sids, D{j}];
            end
            sids = unique(sids);
            idx = find(Sdir(sids)==2); % get post-synapse ids in this ROI
            postsids = sids(idx);
            outnids = unique(StoN(postsids)); % get (post-synapse) traced & orphan body-ids
            % get pre and post neurons in ROI(i)
            logis = ismember(Nid,outnids);
            Nout{i} = Nid(logis & Nstatus==1); % output traced neuron in ROI(i)
    
            idx = find(Sdir(sids)==1); % get pre-synapse ids in this ROI
            presids = sids(idx);
            innids = unique(StoN(presids)); % get (pre-synapse) traced & orphan body-ids
            logis = ismember(Nid,innids);
            Nin{i} = Nid(logis & Nstatus==1); % input traced neuron in ROI(i)
        end
        save(nfile,'Nin','Nout','-v7.3');
    end

    PS = cell(roimax,1); LOC = cell(roimax,1);
    for i=1:roimax
        if isempty(Nout{i}), continue; end
        disp(['pre-syn_loc ' num2str(i) ' / ' num2str(roimax)]);

        % ROI(i) output all cells to pre-synapses for other ROIs
        logi = ismember(StoN,Nout{i}); % find synapses which belong to ROI(i) output neurons
        outsids = Sid(logi);
        idx = find(Sdir(outsids)==1); % get pre-synapse ids of output neurons
        PS{i} = outsids(idx);
        LOC{i} = SlocFc(PS{i},:); % get 3D location in FDA Cal template.
    end
    clear SlocFc; clear Sdir; clear Nid; clear Nstatus;

    % set pool num. this calculation takes time. we need big pool num.
%    delete(gcp('nocreate')); % shutdown pools
%    parpool(32);

    sycountMat = sparse(roimax,roimax);
    countMat = sparse(roimax,roimax);
    parfor i=1:roimax
        if isempty(PS{i}), continue; end
        disp(['get synaptic connections of ROI ' num2str(i) ' / ' num2str(roimax)]);
        tic;
            
        % get connected synapse counts in each ROI (from ROI to connected ROI)
        conSlocFc = LOC{i};
        Vs = zeros(sz(1),sz(2),sz(3),'uint16');
        N = cell(sz(1),sz(2),sz(3));
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                Vs(t(1),t(2),t(3)) = Vs(t(1),t(2),t(3)) + 1;
                N{t(1),t(2),t(3)} = [N{t(1),t(2),t(3)},StoN(PS{i}(j))];
            else
                disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
            end
        end
        B = cell(roimax,1);
        D = cell(roimax,1);
        for j=1:roimax
            if isempty(roiIdxs{j}), continue; end
            B{j} = Vs(roiIdxs{j});
            D{j} = N(roiIdxs{j});
        end

        CX = cell(roimax,1);
        X = sparse(1,roimax);
        S = sparse(1,roimax);
        for j=1:roimax
            if isempty(roiIdxs{j}), continue; end
            S(j) = nansum(B{j},1);

            % get connected neuron counts in each ROI (from ROI to connected ROI)
            nids = [];
            for k=1:length(D{j})
                nids = [nids, D{j}{k}];
            end
            numsyn = groupcounts(nids'); % number of synapse in each neuron
            nids = unique(nids);
            CX{j} = nids(numsyn >= synTh);
            X(j) = length(CX{j}); % synapse number threshold for neurons in ROI(j)
        end
        sycountMat(i,:) = S;
        countMat(i,:) = X;
        t = toc; disp(['t=' num2str(t)]);
    end
    delete(gcp('nocreate')); % shutdown pools
    clear StoN; clear C;
    clear PS; clear LOC;
end

function [weightMat] = makeSCweightMatrixLarge(sycountMat, type)
    roimax = size(sycountMat,1);
    Nout = {}; Nin = {};
    nfile = ['results/cache-' type '_L_Nin_Nout.mat'];
    load(nfile);

    delete(gcp('nocreate')); % shutdown pools
    parpool(32);

    % calculate weight matrix (full, neurons, others)
    weightMat = sparse(roimax,roimax);
    parfor i=1:roimax
        if isempty(Nout{i}), continue; end
        disp(['ROI ' num2str(i) ' / ' num2str(roimax)]);
        outnids = Nout{i};
        S = sparse(1,roimax);
        for j=1:roimax
            if isempty(Nin{j}), continue; end
            innids = Nin{j};
            % find input neuron rate from ROI(i)
            logi = ismember(innids,outnids);
            S(j) = double(sum(logi)) / length(innids); % in-weight (from i to j)
        end
        weightMat(i,:) = S;
    end
    delete(gcp('nocreate')); % shutdown pools
    weightMat = weightMat .* sycountMat;
end
