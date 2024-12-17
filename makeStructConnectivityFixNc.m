% make structural connectivity matrix from neuprint connectivity list (csv).
% countMat: only ROI was transformed. countMat2: synapse points were transformed and re-counted in all ROI. 
% this script should run after makeVoxelROIatlas.m

function makeStructConnectivityFixNc
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
    fname = 'data/hemiroi_connectlist.mat'; scver = 1;
    load(fname);

    if scver <= SCVER
        info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
        aV = niftiread(info); aV(:) = 0;
        cnt = 1;
        roiIdxs = {};
        listing = dir(['atlas/hemiroi/*.nii.gz']);
        for i=1:length(listing)
            V = niftiread(['atlas/hemiroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            idx = find(V>0);
            roiIdxs{i} = idx;
            if any(ismember(turnerIds,i))
                aV(idx) = cnt; cnt = cnt + 1;
            end
            sz = size(V);
        end
%        niftiwrite(aV,'atlas/hemiFlyem52atlasCal.nii','Compressed',true);

        [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'hemiroi');
    end

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
%}

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson matrix csv.
    % extract ROI ids from hemibrain mask
    % branson matrix csv was acquired from Turner et al. (2021).
%%{
    fname = 'data/bransonhemi_connectlist.mat';
    load(fname);

    if scver <= SCVER % recount sycountMat to make compatibility with FlyWire
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
        [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson'); % recount post-synaptic connection of terget ROI
        scver = 5;
    end
    if scver <= SCVER
        scver = scver + 0.1;
        save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
    end
    ids = primaryIds; CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
    figure; imagesc(log(CM)); colorbar; title('branson neurons matrix');
    figure; imagesc(log(SM)); colorbar; title('branson synapses matrix');
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson 7065 atlas.
    % extract ROI ids from hemibrain mask
%%{
    fname = 'data/hemibranson7065_connectlist.mat';
    load(fname);

    if scver <= SCVER % recount sycountMat to make compatibility with FlyWire
        atlV = niftiread('atlas/hemiBranson7065atlasCal.nii.gz');
        roimax = max(atlV(:));
        sz = size(atlV);

        roiIdxs = {};
        for i=1:roimax
            roiIdxs{i} = find(atlV==i);
        end

        primaryIds = 1:roimax;
        roiNum = length(primaryIds);

        [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, 'branson7065');

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
    % make structural connectivity matrix from branson 7065 k-means atlas.
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiBranson7065km' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];
        load(fname);

        if scver <= SCVER % recount sycountMat to make compatibility with FlyWire
            atlV = niftiread(['atlas/hemiBranson7065km' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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
    % make structural connectivity matrix from k-means atlas.
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000 30000]
        idstr = ['hemiCmkm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];
        load(fname);

        if scver <= SCVER 
            atlV = niftiread(['atlas/hemiCmkm' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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
%%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000]
        idstr = ['hemiCmkm' num2str(k) 'r1w1'];
        fname = ['data/' lower(idstr) '_connectlist.mat'];
        load(fname);

        if scver <= SCVER 
            atlV = niftiread(['atlas/hemiCmkm' num2str(k) 'r1w1atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);
    
            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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
    % make structural connectivity matrix from distance based k-means atlas.
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000 5000 10000 15000 20000 30000]
        idstr = ['hemiDistKm' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];
        load(fname);

        if scver <= SCVER 
            atlV = niftiread(['atlas/hemiDistKm' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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
    % make structural connectivity matrix from fully random cluster atlas.
    % extract ROI ids from hemibrain mask
%%{
    for k=[20 30 50 100 200 300 500 1000]
        idstr = ['hemiRand' num2str(k)];
        fname = ['data/' lower(idstr) '_connectlist.mat'];
        load(fname);

        if scver <= SCVER 
            atlV = niftiread(['atlas/hemiRand' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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
        load(fname);

        if scver <= SCVER 
            atlV = niftiread(['atlas/hemiVrand' num2str(k) 'atlasCal.nii.gz']);
            roimax = max(atlV(:));
            sz = size(atlV);
    
            roiIdxs = {};
            for i=1:roimax
                roiIdxs{i} = find(atlV==i);
            end
    
            primaryIds = 1:roimax;
            roiNum = length(primaryIds);

            [ncountMat, sycountMat] = makeSCcountMatrix(roiIdxs, sz, 0.8, 0, lower(idstr));

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

