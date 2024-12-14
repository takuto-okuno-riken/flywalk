% make structural connectivity matrix from neuprint connectivity list (csv).
% countMat: only ROI was transformed. countMat2: synapse points were transformed and re-counted in all ROI. 
% this script should run after makeVoxelROIatlas.m

function makeStructConnectivity
    SCVER = 5;
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from neuprint connectivity list.
    % list data was acquired by c.fetch_roi_connectivity() of neuprint python api.
%%{
    % primary ROIs
    primaryIds = [103	107	20	111	59	68	65	78	34	4	49	51	62	106	87	47	100	24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	74	42	80	1	102	63	95	56];
    roiNum = 114;
    turnerIds = [103	107	20	111	59	68	65	78  49	51	62	106	87	47 100 24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	80	1	102	63	95	56];

    % load matfile
    fname = 'data/flyemroi_connectlist.mat'; scver = 1;
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

        [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'hemiroi');

        save([fname(1:end-4) '_cnids.mat'],'Cnids','-v7.3');
        scver = 3;
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
        dcm = ncountMat(:,:,2) - full(countMat3);
        disp(['ncountMat-countMat3 : '  num2str(sum(abs(dcm),'all'))]);
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
    figure; imagesc(log(CM)); colorbar; title('hemibrain neurons matrix');
    figure; imagesc(log(WM)); colorbar; title('hemibrain all synapse weight matrix');
%}
    % primary, R/L, name order
    load('data/flyemroi.mat');
    ids = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    labelNames = roiname(ids,1);

    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagescLabel(log10(CM), labelNames, 'hemibrain neurons matrix');
    figure; imagescLabel(log10(SM), labelNames, 'hemibrain all synapses matrix');
    WMi = nweightMat(ids,ids,2);
    WMo = outweightMat(ids,ids,2);
    figure; imagescLabel(log(WMi'), labelNames, 'hemibrain all input-rate matrix');
    figure; imagescLabel(log(WMo), labelNames, 'hemibrain all output-rate matrix');

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
    figure; imagesc(log10(CM)); colorbar; title('hemibrain neurons matrix');
    figure; imagesc(log10(SM)); colorbar; title('hemibrain synapses matrix');
%    CM2b = ncountMat(ids,ids,2); SMb = sycountMat(ids,ids,2); WM2b = weightMat2(ids,ids,2); WMob = outweightMat(ids,ids,2);
%    WM3b = WM2b ./ SMb; WM3b(SMb==0) = 0; % pure ROI-input neuron connection weight
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from neuprint connectivity list (flyem hemibrain) by FlyWire EM data.
    % extract ROI ids from hemibrain mask
%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    fname = 'data/flyemroi_fw_connectlist.mat';
    if exist(fname,'file')
        load(fname);
    else
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

        [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, 'hemiroi_fw');

        countMat = []; weightMat = []; scver = 5;
    end
    if scver <= SCVER
        % set same order of FlyEM hemibrain.
        flyemname = ['data/flyemroi_connectlist.mat'];
        cl = load(flyemname);
        primaryIds = cl.primaryIds;
        clear cl;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end
    ids = primaryIds; CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title('hemibrain fw neurons matrix');
    figure; imagesc(log(SM)); colorbar; title('hemibrain fw synapses matrix');
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson matrix csv.
    % extract ROI ids from hemibrain mask
    % branson matrix csv was acquired from Turner et al. (2021).
%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
    fname = 'data/bransonhemi_connectlist.mat';
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

        [countMat2, sycountMat, weightMat2, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson');
        scver = 5;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end
%{
    ids = primaryIds; CM = ncountMat(ids,ids); WM = weightMat(ids,ids);
    figure; imagesc(log(CM)); colorbar; title('branson neurons matrix');
    figure; imagesc(log(WM)); colorbar; title('branson synapse weight matrix');
%}
    ids = primaryIds; CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title('branson neurons matrix');
    figure; imagesc(log(SM)); colorbar; title('branson synapses matrix');
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    % extract ROI ids from hemibrain mask
%{
    fname = 'data/hemibranson7065_connectlist.mat';
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

        [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson7065');

        countMat = []; weightMat = []; scver = 5;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','primaryIds','roiNum','scver','-v7.3');
    end

    ids = primaryIds;
    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title(['branson7065 neurons matrix']);
    figure; imagesc(log(SM)); colorbar; title(['branson7065 synapses matrix']);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    % extract ROI ids from hemibrain mask
%{
    fname = 'data/wirebranson7065_connectlist.mat';
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

        [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, 'wirebranson7065');

        countMat = []; weightMat = []; scver = 5;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','primaryIds','roiNum','scver','-v7.3');
    end

    ids = primaryIds;
    CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title(['wirebranson7065 neurons matrix']);
    figure; imagesc(log(SM)); colorbar; title(['wirebranson7065 synapses matrix']);
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas.
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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
    
            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas by FlyWire EM data
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k) '_fw'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % set same order of FlyEM hemibrain.
            flyemname = ['data/' lower(idstr(1:end-3)) '_connectlist.mat'];
            cl = load(flyemname);
            primaryIds = cl.primaryIds;
            clear cl;
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas (flyWire based clustering).
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['wireBranson7065km' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/wireBranson7065km' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 k-means atlas (flyWire based clustering) by FlyWire EM data
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['wireBranson7065km' num2str(k) '_fw'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/wireBranson7065km' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for cube ROI.
    % extract ROI ids from hemicube4 mask
%{
    for k=[12 8 4]
        idstr = ['hemiCube' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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
        figure; imagesc(log(CM2)); colorbar; title([idstr ' neurons 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for piece ROI.
    % extract ROI ids from hemipiece3 mask
    for k=[12 8 4 3 2]
        idstr = ['hemiPiece' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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
            save(fname,'countMat','weightMat','countMat2','sycountMat','weightMat2','outweightMat','primaryIds','roiNum','scver','-v7.3');
        end
    
        ids = primaryIds;
        CM2 = countMat2(ids,ids); SM = sycountMat(ids,ids);
        figure; imagesc(log(CM2)); colorbar; title([idstr ' neurons 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for fanshape body (FB) and others.
    % extract voxel ids from fanshape-body (FB), eliptic-body (EB), and other atlas.
%{
%    roiids = {[101],[57],[57,51],[51,62,20,111,100]}; % FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
    roiids = {1	5	7	27	30	32	43	52	54	57	59	63	65	67	78	82	89	93	95	100	101	106	113};

    for k = 1:length(roiids)
        idstr = num2str(roiids{k}(1));
        for j=2:length(roiids{k}), idstr=[idstr '-' num2str(roiids{k}(j))]; end

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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
    
            [ncountMat, sycountMat, nweightMat, outweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, ['hemiroi' idstr]);
    
            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = countMat2(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title(['hemiroi ' idstr ' neurons 2 matrix']);
        figure; imagesc(log(SM)); colorbar; title(['hemiroi ' idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from synapse list for all EM ROI voxels (except fibers).
%%{
    clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

        save([fname(1:end-4) '_wm.mat'],'weightMat','primaryIds','roiNum','scver','-v7.3');
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means atlas.
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000 30000]
        idstr = ['hemiCmkm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means (smoothing) atlas.
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000]
        idstr = ['hemiCmkm' num2str(k) 'r1w1'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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
    
            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from k-means atlas by FlyWire EM data
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiCmkm' num2str(k) '_fw'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % set same order of FlyEM hemibrain.
            flyemname = ['data/' lower(idstr(1:end-3)) '_connectlist.mat'];
            cl = load(flyemname);
            primaryIds = cl.primaryIds;
            clear cl;
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas.
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000 30000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from distance based k-means atlas by FlyWire EM data.
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiDistKm' num2str(k) '_fw'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrixFw(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % set same order of FlyEM hemibrain.
            flyemname = ['data/' lower(idstr(1:end-3)) '_connectlist.mat'];
            cl = load(flyemname);
            primaryIds = cl.primaryIds;
            clear cl;
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM)); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM)); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from fully random cluster atlas.
    % extract ROI ids from hemibrain mask
%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiRand' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
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

            [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM),[8 10]); colorbar; title([idstr ' neurons 2 matrix']);
        figure; imagesc(log(SM),[8 12]); colorbar; title([idstr ' synapses matrix']);
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from fully random voxel atlas.
    % extract ROI ids from hemibrain mask

    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000]
        idstr = ['hemiVrand' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];

        clear countMat2; clear ncountMat; clear sycountMat; clear weightMat2; scver = 1;
        if exist(fname,'file')
            load(fname);
        else
            atlV = niftiread(['atlas/hemiVrand' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            if k <= 1000
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
            else
                [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));
                nweightMat = []; outweightMat = []; syweightMat = [];
            end

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            % reorder by tree clustering
            cm2 = ncountMat(:,:,2); cm2(isnan(cm2)) = 0;
            eucD = pdist(cm2,'euclidean');
            Z = linkage(eucD,'ward');
            [H,T,outperm] = dendrogram(Z,k);
            primaryIds = outperm; % use default leaf order
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end

        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log(CM),[8 10]); colorbar; title([idstr ' neurons matrix']);
        figure; imagesc(log(SM),[8 12]); colorbar; title([idstr ' synapses matrix']);
    end
end

function imagescLabel(mat, labelNames, titlestr)
    imagesc(mat); colorbar; daspect([1 1 1]); title(titlestr);
    set(gca,'XTick',1:size(mat,1));
    set(gca,'YTick',1:size(mat,1));
    set(gca,'XTickLabel',labelNames);
    set(gca,'YTickLabel',labelNames);
end
