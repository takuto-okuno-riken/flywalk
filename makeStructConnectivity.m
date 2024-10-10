% make structural connectivity matrix from neuprint connectivity list & branson matrix csv.

function makeStructConnectivity
    % ---------------------------------------------------------------------
    % make structural connectivity matrix from neuprint connectivity list.
    % list data was acquired by c.fetch_roi_connectivity() of neuprint python api.

    % primary ROIs
    primaryIds = [103	107	20	111	59	68	65	78	34	4	49	51	62	106	87	47	100	24	27	43	38	5	57	22	89	101	97	75	50	58	41	113	10	2	32	66	45	30	67	19	76	31	82	93	54	52	8	7	74	42	80	1	102	63	95	56];
    roiNum = 114;

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

    % find primary ROI ids from list 
    % !caution! neurite does not return full primary ROI connections
%    id1 = unique(connectlist(:,1));
%    id2 = unique(connectlist(:,2));
%    id3 = unique([id1, id2]);

    ids = primaryIds; CM = countMat(ids,ids); WM = weightMat(ids,ids);
    figure; imagesc(log(CM)); colorbar; title('himibrain cell count matrix');
    figure; imagesc(log(WM)); colorbar; title('himibrain synapse weight matrix');

    % ---------------------------------------------------------------------
    % make structural connectivity matrix from branson matrix csv.
    % extract ROI ids from hemibrain mask
    % branson matrix csv was acquired from Turner et al. (2021).

    fname = 'data/branson_connectlist.mat';
    if exist(fname,'file')
        load(fname);
    else
        % read hemibrain mask
        Vm = niftiread('data/jrc2018f_flyemhemibrainCal_invFDACal.nii.gz');
        Vm(Vm>0) = 1;
        Vm(Vm<1) = 0;
    
        % read branson atlas (FDA registered)
        Vb = niftiread('data/JRC2018_branson_atlasCal_invFDACal.nii.gz');
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
        save(fname,'connectlist','countMat','weightMat','primaryIds','roiNum');
    end

    ids = primaryIds; CM = countMat(ids,ids); WM = weightMat(ids,ids);
    figure; imagesc(log(CM)); colorbar; title('branson cell count matrix');
    figure; imagesc(log(WM)); colorbar; title('branson synapse weight matrix');
end
