% make structural connectivity matrix with various conditions.
% for FlyEM, synapse number threshold at a connection and synaptic accuracy rate are changed.
% for FlyWire, synapse number threshold at a connection. 
% this script should run after makeVoxelROIatlas.m

function makeStructConnectivityWithConditions
    SCVER = 5;
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of hemibrain primary ROI atlas.

    % primary, R/L, name order
    load('data/flyemroi.mat');
    primaryIds = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
    roiNum = length(primaryIds);
    labelNames = roiname(primaryIds,1);
%{
    rateThs = [50 60 70 80 90];
    synThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(rateThs)
        rateTh = rateThs(r);
        for j=1:length(synThs)
            synTh = synThs(j);

            fname = ['data/flyemroi_hb' num2str(synTh) 'sr' num2str(rateTh) '_connectlist.mat'];
            if exist(fname,'file')
                load(fname);
            else
                roiIdxs = {};
                listing = dir(['atlas/flyemroi/*.nii.gz']);
                for i=1:length(listing)
                    V = niftiread(['atlas/flyemroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
                    idx = find(V>0);
                    roiIdxs{i} = idx;
                    sz = size(V);
                end
        
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrix(roiIdxs, sz, rateTh/100, synTh, ['hemiroi_hb' num2str(synTh) 'sr' num2str(rateTh)]);
        
                countMat = []; weightMat = []; scver = 5;
            end
            if scver <= SCVER
                scver = scver + 0.1;
                save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
            end
        
            ids = primaryIds;
            CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
            figure; imagesc(log10(CM)); colorbar; title(['hemibrain scTh' num2str(synTh) 'srTh' num2str(rateTh) ' neurons matrix']);
            figure; imagesc(log10(SM)); colorbar; title(['hemibrain scTh' num2str(synTh) 'srTh' num2str(rateTh) ' synapses matrix']);
        end
    end
%}
    % ---------------------------------------------------------------------
    % make structural connectivity matrix of hemibrain primary ROI atlas (flyem hemibrain) by FlyWire EM data.
%%{
    rateThs = [50 60 70 80 90 100];
    synThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(rateThs)
        rateTh = rateThs(r);
        for j=1:length(synThs)
            synTh = synThs(j);
    
            clear ncountMat; clear sycountMat; clear nweightMat; scver = 1;
            fname = ['data/flyemroi_fw' num2str(synTh) 'sr' num2str(rateTh) '_connectlist.mat'];
            if exist(fname,'file')
                load(fname);
            else
                roiIdxs = {};
                listing = dir(['atlas/flyemroi/*.nii.gz']);
                for i=1:length(listing)
                    V = niftiread(['atlas/flyemroi/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
                    idx = find(V>0);
                    roiIdxs{i} = idx;
                    sz = size(V);
                end
        
                [ncountMat, sycountMat, nweightMat, outweightMat, syweightMat] = makeSCcountMatrixFw(roiIdxs, sz, rateTh, synTh, ['hemiroi_fw' num2str(synTh) 'sr' num2str(rateTh)]);
        
                countMat = []; weightMat = []; scver = 5;
            end
            if scver <= SCVER
                scver = scver + 0.1;
                save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
            end
            ids = primaryIds;
            CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
            figure; imagesc(log10(CM)); colorbar; title(['flywire scTh' num2str(synTh) ' neurons matrix']);
            figure; imagesc(log10(SM)); colorbar; title(['flywire scTh' num2str(synTh) ' synapses matrix']);
        end
    end

    % ---------------------------------------------------------------------
    % make structural connectivity matrix of hemibrain primary ROI atlas (flyem hemibrain) by taking mean from FlyEM and FlyWire.
    rateTh = 80;
    synThs = [0 5 10 20 30 50 100];
    for j=1:length(synThs)
        synTh = synThs(j);

        fname = ['data/flyemroi_avg' num2str(synTh) '_connectlist.mat'];
        if exist(fname,'file')
            load(fname);
        else
            t1name = ['data/flyemroi_hb' num2str(synTh) 'sr' num2str(rateTh) '_connectlist.mat'];
            t1 = load(t1name);
            t2name = ['data/flyemroi_fw' num2str(synTh) '_connectlist.mat'];
            t2 = load(t2name);
            ncountMat = (t1.ncountMat + t2.ncountMat) / 2;
            nweightMat = (t1.nweightMat + t2.nweightMat) / 2;
            sycountMat = (t1.sycountMat + t2.sycountMat) / 2;
            syweightMat = (t1.syweightMat + t2.syweightMat) / 2;
            outweightMat = (t1.outweightMat + t2.outweightMat) / 2;

            countMat = []; weightMat = []; scver = 5;
        end
        if scver <= SCVER
            scver = scver + 0.1;
            save(fname,'countMat','weightMat','ncountMat','nweightMat','sycountMat','outweightMat','syweightMat','primaryIds','roiNum','scver','-v7.3');
        end
        ids = primaryIds;
        CM = ncountMat(ids,ids,2); SM = sycountMat(ids,ids,2);
        figure; imagesc(log10(CM)); colorbar; title(['fly average scTh' num2str(synTh) ' neurons matrix']);
        figure; imagesc(log10(SM)); colorbar; title(['fly average scTh' num2str(synTh) ' synapses matrix']);
    end
end
