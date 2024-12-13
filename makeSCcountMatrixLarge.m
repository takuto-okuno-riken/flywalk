% make SC neuron & synapse count matrix for huge matrix by hemibrain FlyEM structure data.

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
