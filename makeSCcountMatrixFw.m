% make SC neuron & synapse count matrix by FlyWire structure data.

function [countMat, sycountMat, weightMat, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrixFw(roiIdxs, sz, rateTh, synTh, type)

    % read neuron info (id, type)
%    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)
%    clear Nid; clear Nscore; clear Ntype; 

    % read synapse info
    Sid = []; postNidx = []; preNidx = [];
    load('data/flywire783_synapse.mat');
    score = (cleftScore >= rateTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % read synapse location in FDA
    SpostlocFc = [];
    load('data/flywire783i_sypostloc_fdacal.mat');

    % make presynapse index
    cfile = 'results/flywire783i_sypostCell.mat';
    if exist(cfile,'file')
        load(cfile);
    else
        C = cell(sz(1),sz(2),sz(3));
        for i=1:size(SpostlocFc,1)
            if ~valid(i), continue; end % use only valid pre-post synapse set
            t = ceil(SpostlocFc(i,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                C{t(1),t(2),t(3)} = [C{t(1),t(2),t(3)},i];
            else
                disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
            end
        end
        save(cfile,'C','-v7.3');
    end

    isweight = (nargout >= 3);
    isoutweight = (nargout >= 4);
    issyweight = (nargout >= 5);
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

            % find post synapses in that ROI.
            D = C(roiIdxs{i});
            sididx = [];
            for j=1:length(D)
                sididx = [sididx, D{j}]; % get pre-post synapse set in this ROI
            end
            logis = ismember(Sidx, sididx); % get valid post-synapse ids in this ROI
            rsidx = Sidx(logis & score);

            nidx = postNidx(rsidx);
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            outnidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) traced root-ids

            % get pre and post neurons in ROI(i)
            Nout{i}{2} = outnidx; % output traced neuron in ROI(i)
    
            if isweight
                nidx = preNidx(sididx); 
                numsyn = groupcounts(nidx); % number of synapse in each neuron
                nidx = unique(nidx);
                innidx = nidx(numsyn >= synTh); % get thresholded (pre-synapse) traced root-ids
                Nin{i}{2} = innidx; % input traced neuron in ROI(i)
            end
%{
            % get pre and post synapses in ROI(i)
            Sout{i}{2} = sididx; % output traced synapses in ROI(i)
%}
            % syweightMat is heavy. calculate is if only it is required.
            if issyweight
                logis = ismember(preNidx,Nin{i}{2});
                logis = ismember(sididx,Sidx(logis)); % find may be slow, but Sid will consume memory.
                Sin{i}{2} = sididx(logis); % output traced synapses in ROI(i)
            end
        end
        if issyweight
            save(nfile,'Nin','Nout','Sin','-v7.3');
        else
            save(nfile,'Nin','Nout','-v7.3');
        end
    end

    roimax = length(roiIdxs);
    countMat = nan(roimax,roimax,3,'single');
    sycountMat = nan(roimax,roimax,3,'single');
    if isweight
        weightMat = nan(roimax,roimax,3,'single');
    end
    if isoutweight
        outweightMat = nan(roimax,roimax,3,'single');
    end
    if issyweight
        syweightMat = nan(roimax,roimax,3,'single');
    end

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(14);

    CC = cell(roimax,1);
%    for i=1:roimax
    parfor i=1:roimax
        if isempty(roiIdxs{i}), continue; end
        disp(['get synaptic connections of ROI ' num2str(i) ' / ' num2str(roimax)]);

        % find post synapses in that ROI.
%        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
%        C{x,y,z}

        % only traced neurons are used
        CX = cell(roimax,3);
        for p=2 
            outnidx = Nout{i}{p};

            % ROI(i) output all cells to pre-synapses for other ROIs
            logi = ismember(preNidx,outnidx); % find synapses which belong to ROI(i) output neurons
            sidx = Sidx(logi & valid & score);

            % get connected synapse counts in each ROI (from ROI to connected ROI)
            conSlocFc = SpostlocFc(sidx,:); % get (pre-post) 3D location in FDA Cal template.
            N = cell(sz(1),sz(2),sz(3));
            for j=1:size(conSlocFc,1)
                t = ceil(conSlocFc(j,:));
                if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                    N{t(1),t(2),t(3)} = [N{t(1),t(2),t(3)},preNidx(sidx(j))];
                else
                    disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
                end
            end
            conSlocFc = []; % clear memory;

            % get connected neuron counts in each ROI (from ROI to connected ROI)
            X = zeros(1,roimax,'int32');
            Y = zeros(1,roimax,'int32');
            for j=1:roimax
                D = N(roiIdxs{j});
                nidx = [];
                for k=1:length(D)
                    nidx = [nidx, D{k}];
                end
                numsyn = groupcounts(nidx'); % number of synapse in each neuron
                nidx = unique(nidx);
                CX{j,p} = nidx(numsyn >= synTh);
                X(j) = length(CX{j,p}); % synapse number threshold for neurons in ROI(j)
                ns = numsyn(numsyn >= synTh);
                if ~isempty(ns), Y(j) = nansum(ns,1); end
            end
            countMat(i,:,p) = X;
            sycountMat(i,:,p) = Y;
        end
        CC{i} = CX;
    end
    clear SpostlocFc;

    % calculate weight matrix (full, neurons, others)
    if nargout < 3, return; end
    for p=2
%        for i=1:roimax
        parfor i=1:roimax
            if isempty(Nout{i}), continue; end
            disp(['ROI ' num2str(i) ' / ' num2str(roimax)]);
            outnidx = Nout{i}{p};
            X = zeros(1,roimax,'single');
            Y = zeros(1,roimax,'single');
            SX = zeros(1,roimax,'single');
            for j=1:roimax
                if isempty(Nin{j}), continue; end
                innidx = Nin{j}{p};
                % find input neuron rate from ROI(i)
                logi = ismember(innidx,outnidx);
                X(j) = single(sum(logi)) / length(innidx); % in-weight (from i to j)
                % syweightMat is heavy. calculate is if only it is required.
                if issyweight
                    insidx = Sin{j}{p};
                    logis = ismember(preNidx,innidx(logi));
                    logis = ismember(insidx,Sidx(logis));
                    SX(j) = single(sum(logis)) / length(insidx); % in-synaptic-weight (from i to j)
                end
                % find output neuron rate from ROI(i)
                if isoutweight
                    logi2 = ismember(outnidx,innidx); % actually, same as x(j)
                    Y(j) = single(sum(logi2)) / length(outnidx); % out-weight (from i to j)
                end
            end
            weightMat(i,:,p) = X;
            if isoutweight
                outweightMat(i,:,p) = Y;
            end
            if issyweight
                syweightMat(i,:,p) = SX;
            end
        end
    end

    % pure input and output cells (full, neurons, others count)
    if nargout < 6, return; end
    Ncount = zeros(roimax,2,3,'single');
    for p=2
        for i=1:roimax
            if ~isempty(Nin{i})
                Ncount(i,1,p) = length(Nin{i}{p});
            end
            if ~isempty(Nout{i})
                Ncount(i,2,p) = length(Nout{i}{p});
            end
        end
    end

    % connected neuron ids (cell), only output neurons and others, but not full.
    if nargout < 7, return; end
    Cnids = cell(roimax,roimax,3);
    for p=2
        for i=1:roimax
            for j=1:roimax
                if ~isempty(CC{i})
                    Cnids{i,j,p} = CC{i}{j,p};
                end
            end
        end
    end
end
