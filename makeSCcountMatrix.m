% make SC neuron & synapse count matrix by hemibrain FlyEM structure data.

function [countMat, sycountMat, weightMat, outweightMat, syweightMat, Ncount, Cnids] = makeSCcountMatrix(roiIdxs, sz, rateTh, synTh, type, spiTh, epsilon, minpts, rcdistTh, isrand, calcRange)
    if nargin < 11, calcRange = {1:3, 1:2}; end
    if nargin < 10, isrand = false; end
    if nargin < 9, rcdistTh = 0; end
    if nargin < 8, minpts = 1; end
    if nargin < 7, epsilon = 3000; end
    if nargin < 6, spiTh = 0; end
    calcRange1 = calcRange{1};
    calcRange2 = calcRange{2};

    % read neuron info (id, connection number, size)
    Nid = []; Nstatus = [];
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % read synapse info
    Sdir = []; StoN = []; Srate = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;
    Sid = uint32(1:length(StoN))';
    srate = (Srate >= rateTh); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    % read synapse location in FDA
    SlocFc = [];
    load('data/hemibrain_v1_2_synapseloc_fdacal.mat');

    % read synapse separation index
    if spiTh > 0
        load(['data/hemibrain_v1_2_synapses_sepidx' num2str(synTh) 'sr' num2str(rateTh*100) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);
        splogi = (Spidx >= spiTh);
        if isrand
            rnum = sum(splogi);
            slogi = (srate & straced);
            idx = find(slogi);
            pidx = randperm(length(idx));
            slogi(idx(pidx(rnum+1:end))) = 0;
            splogi = slogi;
        end
        s1spidx = ismember(StoS(:,1),Sid(splogi));
        s2spidx = ismember(StoS(:,2),Sid(splogi));
        ssspidx = (s1spidx & s2spidx);
    else
        ssspidx = (ssrate | true); % no separation index threshold
    end

    if rcdistTh > 0
        load(['data/hemibrain_v1_2_synapses_reci'  num2str(synTh) 'sr' num2str(rateTh*100) '.mat']);
        clear SrcCloseSid;
        rclogi = ~(SrcCloseDist < rcdistTh); % nan should be ignored by <.
        if isrand
            rnum = sum(SrcCloseDist < rcdistTh);
            slogi = srate & straced;
            idx = find(slogi);
            pidx = randperm(length(idx));
            slogi(idx(pidx(rnum+1:end))) = 0;
            rclogi = ~slogi;
        end
        s1rcdist = ismember(StoS(:,1),Sid(rclogi));
        s2rcdist = ismember(StoS(:,2),Sid(rclogi));
        ssrcdist = (s1rcdist & s2rcdist);
    else
        ssrcdist = (ssrate | true); % no reciprocal distance threshold
    end

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
    Sdir(Srate < rateTh) = 0;  % use only accurate synapse more than 'rate'
    clear Srate;
    % use only more separated synapse
    if spiTh > 0
        Sdir((Spidx >= 0) & ~splogi) = 0;  
        clear Spidx;
    end
    if rcdistTh > 0
        Sdir(~rclogi) = 0;  
        clear SrcCloseDist;
    end

    if ~exist('results/cache','dir'), mkdir('results/cache'); end

    isweight = (nargout >= 3);
    isoutweight = (nargout >= 4);
    issyweight = (nargout >= 5);
    roimax = length(roiIdxs);
    nfile = ['results/cache/' type '_Nin_Nout.mat'];
    if exist(nfile,'file')
        load(nfile);
    else
        Nin = cell(roimax,1); Nout = cell(roimax,1);
        Sin = cell(roimax,1); Sout = cell(roimax,1);
        for i=1:roimax
            if isempty(roiIdxs{i}), continue; end
            disp(['Nin_Nout ' num2str(i) ' / ' num2str(roimax)]);
    
            % find post synapses in that ROI.
    %        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
    %        C{x,y,z}
            D = C(roiIdxs{i});
            sids = [];
            for j=1:length(D)
                sids = [sids, D{j}];
            end
            slogi = ismember(Sid,sids);
            postsids = Sid(slogi & Sdir==2); % get valid post-synapse ids in this ROI

            nids = StoN(postsids);
            numsyn = groupcounts(nids); % number of synapse in each neuron
            nids = unique(nids);
            outnids = nids(numsyn >= synTh); % get thresholded (post-synapse) traced & orphan body-ids

            % get pre and post neurons in ROI(i)
            Nout{i}{1} = outnids; % all output cells (including orphan, etc)
            logis = ismember(Nid,outnids);
            Nout{i}{2} = Nid(logis & Nstatus==1); % output traced neuron in ROI(i)
            logis = ismember(outnids, Nout{i}{2});
            Nout{i}{3} = outnids(~logis); % output orphan bodys
    
            if isweight
                presids = Sid(slogi & Sdir==1);  % get valid pre-synapse ids in this ROI
                nids = StoN(presids);
                numsyn = groupcounts(nids); % number of synapse in each neuron
                nids = unique(nids);
                innids = nids(numsyn >= synTh); % get thresholded (pre-synapse) traced & orphan body-ids

                Nin{i}{1} = innids; % all input cells (including orphan, etc)
                logis = ismember(Nid,innids);
                Nin{i}{2} = Nid(logis & Nstatus==1); % input traced neuron in ROI(i)
                logis = ismember(innids, Nin{i}{2});
                Nin{i}{3} = innids(~logis); % input orphan bodys
            end
%{
            % get pre and post synapses in ROI(i)
            Sout{i}{1} = postsids; % all output synapses (including orphan, etc)
            logis = ismember(StoN,Nout{i}{2});
            logis = ismember(postsids,Sid(logis)); % find may be slow, but Sid will consume memory.
            Sout{i}{2} = postsids(logis); % output traced synapses in ROI(i)
            Sout{i}{3} = postsids(~logis); % output orphan bodys' synapses
%}
            % syweightMat is heavy. calculate is if only it is required.
            if issyweight
                logis = ismember(StoN,Nin{i}{1});
                logis = ismember(presids,Sid(logis)); % find may be slow, but Sid will consume memory.
                thpresids = presids(logis);
                Sin{i}{1} = thpresids; % all thresholded input synapses (including orphan, etc)
                logis = ismember(StoN,Nin{i}{2});
                logis = ismember(thpresids,Sid(logis)); % find may be slow, but Sid will consume memory.
                Sin{i}{2} = thpresids(logis); % output thresholded traced synapses in ROI(i)
                Sin{i}{3} = thpresids(~logis); % output thresholded orphan bodys' synapses
            end
        end
        if issyweight
            save(nfile,'Nin','Nout','Sin','-v7.3');
        else
            save(nfile,'Nin','Nout','-v7.3');
        end
    end
    clear Nid; clear Nstatus;

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
    parpool(18);

    CC = cell(roimax,1);
%    for i=1:roimax
    parfor i=1:roimax
        if isempty(roiIdxs{i}), continue; end
        disp(['get synaptic connections of ROI ' num2str(i) ' / ' num2str(roimax)]);

        % find post synapses in that ROI.
%        [x,y,z] = ind2sub(sz,roiIdxs{i}(10)); % check by ind2sub
%        C{x,y,z}

        % three patterns, full (including orphan, etc), neurons, others (orphan, etc)
        CX = cell(roimax,3);
        for p=calcRange1
            outnids = Nout{i}{p};

            % ROI(i) output all cells to pre-synapses for other ROIs
            slogi = ismember(StoN,outnids); % find synapses which belong to ROI(i) output neurons
            presids = Sid(slogi & Sdir==1); % get pre-synapse ids of output neurons
            sslogi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
            if p==1
                cpresids = StoS(sslogi & ssrate & ssspidx & ssrcdist,1);
                cpostsids = StoS(sslogi & ssrate & ssspidx & ssrcdist,2);
            else
                cpresids = StoS(sslogi & ssrate & sstraced & ssspidx & ssrcdist,1);
                cpostsids = StoS(sslogi & ssrate & sstraced & ssspidx & ssrcdist,2);
            end
            [cpostsids, ia] = unique(cpostsids);
            cpresids = cpresids(ia);
            outsids = []; presids = []; % clear memory

            % get connected synapse counts in each ROI (from ROI to connected ROI)
            conSlocFc = SlocFc(cpostsids,:); % get 3D location in FDA Cal template.
            N = cell(sz(1),sz(2),sz(3));
            for j=1:size(conSlocFc,1)
                t = ceil(conSlocFc(j,:));
                if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                    N{t(1),t(2),t(3)} = [N{t(1),t(2),t(3)},StoN(cpresids(j))];
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
                nids = [];
                for k=1:length(D)
                    nids = [nids, D{k}];
                end
                numsyn = groupcounts(nids'); % number of synapse in each neuron
                nids = unique(nids);
                CX{j,p} = nids(numsyn >= synTh);
                X(j) = length(CX{j,p}); % synapse number threshold for neurons in ROI(j)
                ns = numsyn(numsyn >= synTh);
                if ~isempty(ns), Y(j) = nansum(ns,1); end
            end
            countMat(i,:,p) = X;
            sycountMat(i,:,p) = Y;
        end
        CC{i} = CX;
    end
    clear SlocFc; clear Sdir;

    % calculate weight matrix (full, neurons, others)
    if nargout < 3, return; end
    for p=calcRange2
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
                % find input neuron rate from ROI(i)
                logi = ismember(innids,outnids);
                X(j) = single(sum(logi)) / length(innids); % in-weight (from i to j)
                % syweightMat is heavy. calculate is if only it is required.
                if issyweight
                    insids = Sin{j}{p};
                    logis = ismember(StoN,innids(logi));
                    logis = ismember(insids,Sid(logis));
                    SX(j) = single(sum(logis)) / length(insids); % in-synaptic-weight (from i to j)
                end
                % find output neuron rate from ROI(i)
                if isoutweight
                    logi2 = ismember(outnids,innids); % actually, same as x(j)
                    Y(j) = single(sum(logi2)) / length(outnids); % out-weight (from i to j)
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
    for p=calcRange2
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
