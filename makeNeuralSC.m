% make neural Struct Connectivity data.

function makeNeuralSC
    % DBscan param
    epsilon = 5; % micro meter. almost 2 voxels.
    minpts = 3; % set 1, but isolated synapse will be ignored

    % check neural input & output voxels (FlyEM)
    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    scTh = 60; synTh = 5; % almost flywire codex compatible setting
%    scTh = 0; synTh = 0; % for cheking neuPRINT+ compatible
    checkNeuralInputOutputVoxels(synTh, scTh/100);

    checkNeuralInputOutputDistance('hemi', synTh, scTh);

    checkNeuralDBScan('hemi', synTh, scTh, epsilon, minpts);

    checkNeuralReciprocalConnections(synTh, scTh/100);

    % check neural input & output voxels (FlyWire)
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible

    checkNeuralInputOutputVoxelsFw(synTh, scTh);

    checkNeuralInputOutputDistance('wire', synTh, scTh);

    checkNeuralDBScan('wire', synTh, scTh, epsilon, minpts);

    checkNeuralReciprocalConnectionsFw(synTh, scTh);
end

function checkNeuralInputOutputVoxels(synTh, confTh)

    if ~exist('results/neuralsc','dir'), mkdir('results/neuralsc'); end

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralInOutVoxels.mat'];
    if exist(fname,'file'), return; end

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    Sdir = []; StoN = []; Srate = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;
    Sid = uint32(1:length(StoN))';
    srate = (Srate >= confTh); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < confTh) = 0;  % use only accurate synapse more than 'rate'
    clear Srate;

    % FlyEM read synapse location in FDA
    load('data/synapseloc_fdacal.mat');

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0;
    sz = size(Vt);

    % count pre (to post) and post synapse count
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    outIdx = cell(nlen,1);
    outCount = cell(nlen,1);
    inIdx = cell(nlen,1);
    inCount = cell(nlen,1);
    for i=1:nlen
        logi = ismember(StoN,tracedNids(i)); % find synapses which belong to target neuron
        presids = Sid(logi & Sdir==1); % output of a neuron
        sslogi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
        cpostsids = StoS(sslogi & ssrate & sstraced,2); % post-sid is unique
        if synTh > 0 && ~isempty(cpostsids)
            nids = StoN(cpostsids);
            numsyn = groupcounts(nids); % number of synapse in each neuron
            nids = unique(nids);
            nids2 = nids(numsyn >= synTh); % get thresholded (post-synapse) neuron ids
            logi2 = ismember(StoN,nids2);
            slogi = ismember(Sid,cpostsids);
            cpostsids = Sid(slogi & logi2);
        end

        postsids = Sid(logi & Sdir==2); % input of a neuron
        if synTh > 0 && ~isempty(postsids)
            sslogi = ismember(StoS(:,2),postsids); % get post-synapse to connected pre-synapse
            cpresids = StoS(sslogi & ssrate & sstraced,1); % pre-sid is not unique, but need to keep post-synapse count for FlyWire compatibility
            nids = StoN(cpresids);
            numsyn = groupcounts(nids); % number of (connected post) synapse in each (pre-synapse) neuron
            nids = unique(nids);
            nids2 = nids(numsyn >= synTh); % get thresholded (pre-synapse) neuron ids
            logi2 = ismember(StoN,nids2);
            slogi = ismember(Sid,cpresids);
            cpresids = Sid(slogi & logi2);
            sslogi2 = ismember(StoS(:,1),cpresids);
            postsids2 = StoS(sslogi2,2);  % get post-synapses to connect cpresids (including other neuron's post-synapses)
            plogi = ismember(postsids,postsids2);
            postsids = postsids(plogi);
        end
        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : process(' num2str(i) ') nid=' num2str(tracedNids(i)) ' postsids=' num2str(length(postsids)) ' presids=' num2str(length(presids))]);
    
        % get connected post-synapse counts from pre-synapses (output)
        conSlocFc = SlocFc(cpostsids,:); V = Vt; % get 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        idx = find(V>0);
        outIdx{i} = idx;
        outCount{i} = V(idx);

        % get post-synapse count of neuron (input)
        conSlocFc = SlocFc(postsids,:); V = Vt; % get 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        idx = find(V>0);
        inIdx{i} = idx;
        inCount{i} = V(idx);
    end
    save(fname,'inCount','inIdx','outCount','outIdx','tracedNids','-v7.3');
end

function checkNeuralInputOutputVoxelsFw(synTh, scoreTh)

    fname = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat'];
    if exist(fname,'file'), return; end

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load('data/flywire783_synapse.mat');
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % read synapse location in FDA
    load('data/flywire783i_sypostloc_fdacal.mat');

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0;
    sz = size(Vt);

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    tracedNidx = 1:nlen;
    outIdx = cell(nlen,1);
    outCount = cell(nlen,1);
    inIdx = cell(nlen,1);
    inCount = cell(nlen,1);
    for i=1:nlen
        logi = ismember(preNidx,i); % find pre-synapses which belong to target neurons
        if synTh > 0                % for checking flywire codex compatible
            nidx=postNidx(logi);        % get connected post-synapse neurons
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            outnidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(postNidx,outnidx);
            presidx = Sidx(logi & logi2 & valid & score);
        else
            presidx = Sidx(logi & valid & score);
        end

        logi = ismember(postNidx,i); % find post-synapses which belong to target neurons
        if synTh > 0                 % for checking flywire codex compatible
            nidx=preNidx(logi);         % get connected pre-synapse neurons
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            innidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(preNidx,innidx);
            postsidx = Sidx(logi & logi2 & valid & score);
        else
            postsidx = Sidx(logi & valid & score);
        end
        disp(['wire' num2str(synTh) 'sr' num2str(scoreTh) ' : process(' num2str(i) ') nid=' num2str(Nid(i)) ' postsids=' num2str(length(postsidx)) ' presids=' num2str(length(presidx))]);

        % get connected post-synapse counts from APL pre-synapses (output)
        conSlocFc = SpostlocFc(presidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        idx = find(V>0);
        outIdx{i} = idx;
        outCount{i} = V(idx);

        % get post-synapse count of APL neuron (input)
        conSlocFc = SpostlocFc(postsidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        idx = find(V>0);
        inIdx{i} = idx;
        inCount{i} = V(idx);
    end
    save(fname,'inCount','inIdx','outCount','outIdx','tracedNidx','-v7.3');
end

function checkNeuralInputOutputDistance(scname, synTh, confTh)
    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(confTh) '_neuralInOutDistance.mat'];
    if exist(fname,'file'), return; end

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % load input output voxel info
    niofname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(confTh) '_neuralInOutVoxels.mat'];
    load(niofname);

    D = cell(length(inIdx),1);
    Didx = cell(length(inIdx),2);
    for i=1:length(inIdx)
        if isempty(inIdx{i}) || isempty(outIdx{i}), continue; end

        scinidx = inIdx{i};
        scoutidx = outIdx{i};
        vinidx = V(scinidx); scinidx(vinidx==0)=[];
        voutidx = V(scoutidx); scoutidx(voutidx==0)=[];
        if length(scinidx)==0 || length(scoutidx)==0, continue; end

        scinlen = length(scinidx);
        scoutlen = length(scoutidx);
        disp(['dist ' scname num2str(synTh) 'sr' num2str(confTh) ' : process i=' num2str(i) ' invox=' num2str(scinlen) ' outvox=' num2str(scoutlen)]);

        Sin = zeros(scinlen,1,3,'single');
        for j=1:scinlen
            [x,y,z] = ind2sub(sz,scinidx(j));
            Sin(j,1,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
        end
        Sin = repmat(Sin,[1 scoutlen 1]);
        Sout = zeros(1,scoutlen,3,'single');
        for j=1:scoutlen
            [x,y,z] = ind2sub(sz,scoutidx(j));
            Sout(1,j,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
        end
        Sout = repmat(Sout,[scinlen 1 1]);
        D{i} = sqrt(sum((Sin - Sout).^2,3));

        if scinlen > 2
            eucD = pdist(D{i},'euclidean');
            Z = linkage(eucD,'single');
            [T,ids1] = dendrogramNoplot(Z,scinlen);
        else
            ids1 = 1:scinlen;
        end
        if scoutlen > 2
            eucD = pdist(D{i}(ids1,:)','euclidean');
            Z = linkage(eucD,'single');
            [T,ids2] = dendrogramNoplot(Z,scoutlen);
        else
            ids2 = 1:scoutlen;
        end
        Didx{i,1} = ids1;
        Didx{i,2} = ids2;
%        figure; imagesc(D{i}(ids1,ids2)); colorbar;
    end
    save(fname,'D','Didx','-v7.3');
end

function checkNeuralDBScan(scname, synTh, confTh, epsilon, minpts)
    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(confTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat'];
    if exist(fname,'file'), return; end

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    switch(scname)
    case 'hemi'
        % FlyEM read neuron info (id, connection number, size)
        load('data/hemibrain_v1_2_neurons.mat');
        clear Nconn; clear Ncrop; clear Nsize; 
    case 'wire'
        load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)
    end

    % load input output voxel info
    niofname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(confTh) '_neuralInOutVoxels.mat'];
    load(niofname);

    DBidx = cell(length(inIdx),1);
    inlen = cell(length(inIdx),1);
    for i=1:length(inIdx)
        if isempty(inIdx{i}) || isempty(outIdx{i}), continue; end

        scinidx = inIdx{i};
        scoutidx = outIdx{i};
        vinidx = V(scinidx);
        voutidx = V(scoutidx);
        scinlen = length(scinidx);
        scoutlen = length(scoutidx);

        X = zeros(scinlen+scoutlen,3,'single');
        for j=1:scinlen
            [x,y,z] = ind2sub(sz,scinidx(j));
            X(j,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
        end
        for j=scinlen+1:scinlen+scoutlen
            [x,y,z] = ind2sub(sz,scoutidx(j-scinlen));
            X(j,:) = [x y z] .* [2.45, 2.28, 3.715]; % * voxel size
        end
        dbidx = int16(dbscan(X,epsilon,minpts));
        inlen{i} = scinlen;

        % remove out of mask voxels from cluster
        dbidx([vinidx;voutidx]==0) = -2;

        % reorder cluster number by Y axis
        mcls = max(dbidx);
        sdbidx = dbidx;
        topYs = ones(mcls,1) * 1e3; % set dummy Y
        for j=1:mcls
            Yt = X(dbidx==j,2);
            if ~isempty(Yt), topYs(j) = min(Yt); end % head top is zero direction
        end
        [ys,icls] = sort(topYs);
        for j=1:mcls
            sdbidx(dbidx==icls(j)) = j;
        end
        DBidx{i} = sdbidx;

        cls = unique(sdbidx);
        disp(['dbscan ' scname num2str(synTh) 'sr' num2str(confTh) ' : process (' num2str(i) ') nid=' num2str(Nid(i)) ' cls=' num2str(length(cls(cls>0))) ' invox=' num2str(scinlen) ' outvox=' num2str(scoutlen)]);
    end
    save(fname,'inlen','DBidx','-v7.3');
end

function checkNeuralReciprocalConnections(synTh, confTh)

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    if exist(fname,'file'), return; end
    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    Sdir = []; StoN = []; Srate = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;
    Sid = uint32(1:length(StoN))';
    srate = (Srate >= confTh); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < confTh) = 0;  % use only accurate synapse more than 'rate'
    clear Srate;

    % count pre (to post) and post synapse count
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    inNids = cell(nlen,1);
    outNids = cell(nlen,1);
    rcNids = cell(nlen,1);
    for i=1:nlen
        logi = ismember(StoN,tracedNids(i)); % find synapses which belong to target neuron
        presids = Sid(logi & Sdir==1); % output of a neuron
        sslogi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
        cpostsids = StoS(sslogi & ssrate & sstraced,2); % post-sid is unique
        if synTh > 0 && ~isempty(cpostsids)
            nids = StoN(cpostsids);
            numsyn = groupcounts(nids); % number of synapse in each neuron
            nids = unique(nids);
            nids2 = nids(numsyn >= synTh); % get thresholded (post-synapse) neuron ids
            logi2 = ismember(StoN,nids2);
            slogi = ismember(Sid,cpostsids);
            cpostsids = Sid(slogi & logi2);
        end
        outnids = StoN(cpostsids);
        outnids = unique(outnids);

        postsids = Sid(logi & Sdir==2); % input of a neuron
        sslogi = ismember(StoS(:,2),postsids); % get post-synapse to connected pre-synapse
        cpresids = StoS(sslogi & ssrate & sstraced,1); % pre-sid is not unique, but need to keep post-synapse count for FlyWire compatibility
        if synTh > 0 && ~isempty(postsids)
            nids = StoN(cpresids);
            numsyn = groupcounts(nids); % number of (connected post) synapse in each (pre-synapse) neuron
            nids = unique(nids);
            nids2 = nids(numsyn >= synTh); % get thresholded (pre-synapse) neuron ids
            logi2 = ismember(StoN,nids2);
            slogi = ismember(Sid,cpresids);
            cpresids = Sid(slogi & logi2);
        end
        innids = StoN(cpresids);
        innids = unique(innids);

        logis = ismember(outnids,innids);
        rcNids{i} = outnids(logis); % reciprocally connected neurons.
        outNids{i} = outnids;
        inNids{i} = innids;

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : reciprocal process(' num2str(i) ') nid=' num2str(tracedNids(i)) ' in/out/reci=' num2str(length(inNids{i})) '/' num2str(length(outNids{i})) '/' num2str(length(rcNids{i})) ...
            ' postsids=' num2str(length(postsids)) ' presids=' num2str(length(presids))]);
    end
    save(fname,'rcNids','-v7.3');
    save(fnameNin,'inNids','outNids','-v7.3');
end

function checkNeuralReciprocalConnectionsFw(synTh, scoreTh)

    fname = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralReciprocalConnections.mat'];
    if exist(fname,'file'), return; end
    fnameNin = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load('data/flywire783_synapse.mat');
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    inNidx = cell(nlen,1);
    outNidx = cell(nlen,1);
    rcNidx = cell(nlen,1);
    for i=1:nlen
        logi = ismember(preNidx,i); % find pre-synapses which belong to target neurons
        if synTh > 0                % for checking flywire codex compatible
            nidx=postNidx(logi);        % get connected post-synapse neurons    
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            outnidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(postNidx,outnidx);
            outnidx = postNidx(logi & logi2 & valid & score);
        else
            outnidx = postNidx(logi & valid & score); % get connected neuron index
        end
        outnidx = unique(outnidx);

        logi = ismember(postNidx,i); % find post-synapses which belong to target neurons
        if synTh > 0                 % for checking flywire codex compatible
            nidx=preNidx(logi);         % get connected pre-synapse neurons
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            innidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(preNidx,innidx);
            innidx = preNidx(logi & logi2 & valid & score);
        else
            innidx = preNidx(logi & valid & score);
        end
        innidx = unique(innidx);

        logis = ismember(outnidx,innidx);
        rcNidx{i} = outnidx(logis); % reciprocally connected neurons.
        outNidx{i} = outnidx;
        inNidx{i} = innidx;

        disp(['wire' num2str(synTh) 'sr' num2str(scoreTh) ' : reciprocal process(' num2str(i) ') nid=' num2str(Nid(i)) ' in/out/reci=' num2str(length(inNidx{i})) '/' num2str(length(outNidx{i})) '/' num2str(length(rcNidx{i}))]);
    end
    save(fname,'rcNidx','-v7.3');
    save(fnameNin,'inNidx','outNidx','-v7.3');
end
