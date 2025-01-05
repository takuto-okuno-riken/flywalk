% make neural Struct Connectivity data.

function makeNeuralSC
    % DBscan param
    epsilon = 5; % micro meter. almost 2 voxels.
    minpts = 3; % set 1, but isolated synapse will be ignored

    % check neural input & output voxels (FlyEM)
    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    scTh = 60; synTh = 5; % almost flywire codex compatible setting
%    scTh = 80; synTh = 5; % high confidence & connection setting.
%    scTh = 0; synTh = 0; % for cheking neuPRINT+ compatible
    checkNeuralInputOutputVoxels(synTh, scTh/100);

    checkNeuralInputOutputDistance('hemi', synTh, scTh);

    checkNeuralDBScan('hemi', synTh, scTh, epsilon, minpts);

    checkNeuralReciprocalConnections(synTh, scTh/100);
%{
    Cnames = {'SMP352'};
    Cnids = {[266187532, 266528078, 266528086, 296194535, 328593903, 5813009926]};
    checkNeuralNamedReciprocalConnections(Cnames, Cnids, synTh, scTh/100)
%}
    checkNeuralAutoConnections(synTh, scTh/100);
    
    checkNeuralTriFeedforward(synTh, scTh/100);

    checkNeuralTriUnicycle(synTh, scTh/100);

%    checkNeuralQuadFeedforward(synTh, scTh/100);

    % check neural input & output voxels (FlyWire)
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
%    scTh = 130; synTh = 5; % high confidence & connection setting.

    checkNeuralInputOutputVoxelsFw(synTh, scTh);

    checkNeuralInputOutputDistance('wire', synTh, scTh);

    checkNeuralDBScan('wire', synTh, scTh, epsilon, minpts);

    checkNeuralReciprocalConnectionsFw(synTh, scTh);

    checkNeuralAutoConnectionsFw(synTh, scTh);
    
    checkNeuralTriFeedforwardFw(synTh, scTh);

    checkNeuralTriUnicycleFw(synTh, scTh);
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
    load('data/hemibrain_v1_2_synapseloc_fdacal.mat');

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
    rcpreSids = cell(nlen,1);
    rcpostSids = cell(nlen,1);
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

        % find reciprocal synapses of neuron(i)
        if ~isempty(rcNids{i})
            rclogi = ismember(StoN,rcNids{i}); % pre & post synapses of reciprocal neurons
            sslogi1 = ismember(StoS(:,1),Sid(rclogi)); % get pre-synapse
            cpostsids1 = StoS(sslogi1 & ssrate & sstraced,2); % post-sid is unique
            sslogi2 = ismember(StoS(:,2),Sid(rclogi)); % get post-synapse
            cpresids2 = StoS(sslogi2 & ssrate & sstraced,1); % pre-sid is not unique
            logis = ismember(postsids,cpostsids1);
            rcpostSids{i} = postsids(logis);
            logis = ismember(presids,cpresids2);
            rcpreSids{i} = presids(logis);
        end

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : reciprocal process(' num2str(i) ') nid=' num2str(tracedNids(i)) ' in/out/reci=' num2str(length(inNids{i})) '/' num2str(length(outNids{i})) '/' num2str(length(rcNids{i})) ...
            ' post/pre=' num2str(length(postsids)) '/' num2str(length(presids)) ' rcpost/rcpre=' num2str(length(rcpostSids{i})) '/' num2str(length(rcpreSids{i}))]);
    end
    save(fname,'rcNids','rcpostSids','rcpreSids','-v7.3');
    save(fnameNin,'inNids','outNids','-v7.3');
end

function checkNeuralNamedReciprocalConnections(Cnames, Cnids, synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    load(fname);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize;
    tracedNids = Nid(Nstatus==1);
    Nidx = 1:length(tracedNids);

    for i=1:length(Cnids)
        logis = ismember(tracedNids,Cnids{i});
        idxs = Nidx(logis);
        rcnids = [];
        for j=1:length(idxs)
            nid = tracedNids(idxs(j));
            rcn = rcNids{idxs(j)};
            disp([Cnames{i} ' ' num2str(nid) ' has reciprocal neurons : ' num2str(rcn')]);
            rcnids = [rcnids; rcn];
        end
        numnid = groupcounts(rcnids);
        rcnids = unique(rcnids);
        rcnids = rcnids(numnid>1); nnid = numnid(numnid>1);
        str = [];
        for j=1:length(rcnids), str=[str num2str(rcnids(j)) ' (' num2str(nnid(j)) '), ']; end
        disp([Cnames{i} ' multiple reciprocal neurons : ' str]);
    end
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
    rcpreSidx = cell(nlen,1);
    rcpostSidx = cell(nlen,1);
    for i=1:nlen
        prlogi = ismember(preNidx,i); % find pre-synapses which belong to target neurons
        if synTh > 0                % for checking flywire codex compatible
            nidx=postNidx(prlogi);        % get connected post-synapse neurons    
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            outnidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(postNidx,outnidx);
            outnidx = postNidx(prlogi & logi2 & valid & score);
        else
            outnidx = postNidx(prlogi & valid & score); % get connected neuron index
        end
        outnidx = unique(outnidx);

        pologi = ismember(postNidx,i); % find post-synapses which belong to target neurons
        if synTh > 0                 % for checking flywire codex compatible
            nidx=preNidx(pologi);         % get connected pre-synapse neurons
            numsyn = groupcounts(nidx); % number of synapse in each neuron
            nidx = unique(nidx);
            innidx = nidx(numsyn >= synTh); % get thresholded (post-synapse) neuron index
            logi2 = ismember(preNidx,innidx);
            innidx = preNidx(pologi & logi2 & valid & score);
        else
            innidx = preNidx(pologi & valid & score);
        end
        innidx = unique(innidx);

        logis = ismember(outnidx,innidx);
        rcNidx{i} = outnidx(logis); % reciprocally connected neurons.
        outNidx{i} = outnidx;
        inNidx{i} = innidx;

        % find reciprocal synapses
        if ~isempty(rcNidx{i})
            rclogi1 = ismember(preNidx,rcNidx{i}); % pre-synapse of reciprocal neurons
            rclogi2 = ismember(postNidx,rcNidx{i}); % post-synapse of reciprocal neurons
            rcpostSidx{i} = Sidx(prlogi & rclogi2 & valid & score);
            rcpreSidx{i} = Sidx(pologi & rclogi1 & valid & score); 
        end

        disp(['wire' num2str(synTh) 'sr' num2str(scoreTh) ' : reciprocal process(' num2str(i) ') nid=' num2str(Nid(i)) ' in/out/reci=' num2str(length(inNidx{i})) '/' num2str(length(outNidx{i})) '/' num2str(length(rcNidx{i})) ...
            'rcpost/rcpre=' num2str(length(rcpostSidx{i})) '/' num2str(length(rcpreSidx{i}))]);
    end
    save(fname,'rcNidx','rcpostSidx','rcpreSidx','-v7.3');
    save(fnameNin,'inNidx','outNidx','-v7.3');
end

function checkNeuralAutoConnections(synTh, confTh)

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    if ~exist(fname,'file'), return; end
    load(fname);
    if exist('rcNids','var'), return; end

    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];
    load(fnameNin);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % count pre (to post) and post synapse count
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    autoNids = [];
    for i=1:nlen
        outnids = outNids{i};
        if ~isempty(outnids)
            if any(outnids==tracedNids(i))
                autoNids = [autoNids tracedNids(i)];
            end
        end
    end
    save(fname,'rcNids','autoNids','-v7.3');
end

function checkNeuralAutoConnectionsFw(synTh, scoreTh)

    fname = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralReciprocalConnections.mat'];
    if ~exist(fname,'file'), return; end
    load(fname);
    if exist('rcNidx','var'), return; end

    fnameNin = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(fnameNin);

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    autoNidx = [];
    for i=1:nlen
        outnidx = outNidx{i};
        if ~isempty(outnidx)
            if any(outnidx==i)
                autoNidx = [autoNidx i];
            end
        end
    end
    save(fname,'rcNidx','autoNidx','-v7.3');
end

function checkNeuralTriFeedforward(synTh, confTh)

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralTriFeedforward.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % extract pure output neurons (no reciprocal)
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    poutNids = cell(nlen,1);
    for i=1:nlen
        rcnids = rcNids{i};
        outnids = outNids{i};
        if ~isempty(outnids)
            logis = ismember(outnids,rcnids);
            poutNids{i} = outnids(~logis); % pure out Nids
        end
    end

    % extract triangle feed forward neurons
    tripffNids = cell(nlen,1); % pure output
    trieffNids = cell(nlen,1); % extra (including reciprocal)
    for i=1:nlen
        poutnids = poutNids{i}; % pure output
        if isempty(poutnids), continue; end

        tlogi = zeros(length(poutnids),1,'logical'); % init zero logis
        tlogi2 = zeros(length(poutnids),1,'logical'); % init zero logis
        for j=1:length(poutnids)
            logis = ismember(tracedNids,poutnids(j));
            poutnids2 = poutNids(logis); % pure output. this should extract one cell
            outnids2 = outNids(logis);   % all output (including reciprocal)
%            idx = find(tracedNids==poutnids(j)); % find version.
%            poutnids2 = {poutNids{idx}};
            if ~isempty(outnids2)
                plogi = ismember(poutnids,poutnids2{1});
                alogi = ismember(poutnids,outnids2{1});
                tlogi = tlogi | plogi;
                tlogi2 = tlogi2 | alogi;
            end
        end
        tripffNids{i} = poutnids(tlogi);
        trieffNids{i} = poutnids(tlogi2 & ~tlogi);

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : tri feedforward process(' num2str(i) ') nid=' num2str(tracedNids(i)) ...
            ' out/pout/pff/eff=' num2str(length(outNids{i})) '/' num2str(length(poutNids{i})) '/' num2str(length(tripffNids{i})) '/' num2str(length(trieffNids{i}))]);
    end
    save(fname,'tripffNids','trieffNids','-v7.3');
end

function checkNeuralTriFeedforwardFw(synTh, scoreTh)

    fname = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralTriFeedforward.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralReciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % extract pure output neurons (no reciprocal)
    nlen = length(Nid);
    Nidx = (1:nlen)';
    poutNidx = cell(nlen,1);
    for i=1:nlen
        rcnidx = rcNidx{i};
        outnidx = outNidx{i};
        if ~isempty(outnidx)
            logis = ismember(outnidx,rcnidx);
            poutNidx{i} = outnidx(~logis); % pure out Nidx
        end
    end

    % extract triangle feed forward neurons
    tripffNidx = cell(nlen,1); % pure output
    trieffNidx = cell(nlen,1); % extra (including reciprocal)
    for i=1:nlen
        poutnidx = poutNidx{i};
        if isempty(poutnidx), continue; end

        tlogi = zeros(length(poutnidx),1,'logical'); % init zero logis
        tlogi2 = zeros(length(poutnidx),1,'logical'); % init zero logis
        for j=1:length(poutnidx)
            logis = ismember(Nidx,poutnidx(j));
            poutnidx2 = poutNidx(logis); % pure output. this should extract one cell
            outnidx2 = outNidx(logis);   % all output (including reciprocal)
            if ~isempty(outnidx2)
                plogi = ismember(poutnidx,poutnidx2{1});
                alogi = ismember(poutnidx,outnidx2{1});
                tlogi = tlogi | plogi;
                tlogi2 = tlogi2 | alogi;
            end
        end
        tripffNidx{i} = poutnidx(tlogi);
        trieffNidx{i} = poutnidx(tlogi2 & ~tlogi);

        disp(['wire' num2str(synTh) 'sr' num2str(scoreTh) ' : tri feedforward process(' num2str(i) ') nid=' num2str(Nid(i)) ...
            ' out/pout/pff/eff=' num2str(length(outNidx{i})) '/' num2str(length(poutNidx{i})) '/' num2str(length(tripffNidx{i})) '/' num2str(length(trieffNidx{i}))]);
    end
    save(fname,'tripffNidx','trieffNidx','-v7.3');
end

function checkNeuralTriUnicycle(synTh, confTh)

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralTriUnicycle.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % extract pure output neurons (no reciprocal)
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    poutNids = cell(nlen,1);
    pinNids = cell(nlen,1);
    for i=1:nlen
        rcnids = rcNids{i};
        outnids = outNids{i};
        innids = inNids{i};
        if ~isempty(outnids)
            logis = ismember(outnids,rcnids);
            poutNids{i} = outnids(~logis); % pure output Nids
        end
        if ~isempty(innids)
            logis = ismember(innids,rcnids);
            pinNids{i} = innids(~logis); % pure input Nids
        end
    end

    % extract triangle unicicle neurons
    tripucNids = cell(nlen,1);
    trieucNids = cell(nlen,1);
    for i=1:nlen
        poutnids = poutNids{i};
        if isempty(poutnids), continue; end
        pinnids = pinNids{i};
        if isempty(pinnids), continue; end

        tlogi = zeros(length(pinnids),1,'logical'); % init zero logis
        tlogi2 = zeros(length(pinnids),1,'logical'); % init zero logis
        for j=1:length(poutnids)
            logis = ismember(tracedNids,poutnids(j));
            poutnids2 = poutNids(logis); % pure output. this should extract one cell
            outnids2 = outNids(logis);   % all output (including reciprocal)
            if ~isempty(outnids2)
                plogi = ismember(pinnids,poutnids2{1});
                alogi = ismember(pinnids,outnids2{1});
                tlogi = tlogi | plogi;
                tlogi2 = tlogi2 | alogi;
            end
        end
        tripucNids{i} = pinnids(tlogi);
        trieucNids{i} = pinnids(tlogi2 & ~tlogi);

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : tri unicycle process(' num2str(i) ') nid=' num2str(tracedNids(i)) ...
            ' in/pin/puc/euc=' num2str(length(inNids{i})) '/' num2str(length(pinNids{i})) '/' num2str(length(tripucNids{i})) '/' num2str(length(trieucNids{i}))]);
    end
    save(fname,'tripucNids','trieucNids','-v7.3');
end

function checkNeuralTriUnicycleFw(synTh, scoreTh)

    fname = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralTriUnicycle.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralReciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % extract pure output neurons (no reciprocal)
    nlen = length(Nid);
    Nidx = (1:nlen)';
    poutNidx = cell(nlen,1);
    pinNidx = cell(nlen,1);
    for i=1:nlen
        rcnidx = rcNidx{i};
        outnidx = outNidx{i};
        innidx = inNidx{i};
        if ~isempty(outnidx)
            logis = ismember(outnidx,rcnidx);
            poutNidx{i} = outnidx(~logis); % pure output Nids
        end
        if ~isempty(innidx)
            logis = ismember(innidx,rcnidx);
            pinNidx{i} = innidx(~logis); % pure input Nids
        end
    end

    % extract triangle feed forward neurons
    tripucNidx = cell(nlen,1);
    trieucNidx = cell(nlen,1);
    for i=1:nlen
        poutnidx = poutNidx{i};
        if isempty(poutnidx), continue; end
        pinnidx = pinNidx{i};
        if isempty(pinnidx), continue; end

        tlogi = zeros(length(pinnidx),1,'logical'); % init zero logis
        tlogi2 = zeros(length(pinnidx),1,'logical'); % init zero logis
        for j=1:length(poutnidx)
            logis = ismember(Nidx,poutnidx(j));
            poutnidx2 = poutNidx(logis); % pure output. this should extract one cell
            outnidx2 = outNidx(logis);   % all output (including reciprocal)
            if ~isempty(outnidx2)
                plogi = ismember(pinnidx,poutnidx2{1});
                alogi = ismember(pinnidx,outnidx2{1});
                tlogi = tlogi | plogi;
                tlogi2 = tlogi2 | alogi;
            end
        end
        tripucNidx{i} = pinnidx(tlogi);
        trieucNidx{i} = pinnidx(tlogi2 & ~tlogi);

        disp(['wire' num2str(synTh) 'sr' num2str(scoreTh) ' : tri unicycle process(' num2str(i) ') nid=' num2str(Nid(i)) ...
            ' in/pin/puc/euc=' num2str(length(inNidx{i})) '/' num2str(length(pinNidx{i})) '/' num2str(length(tripucNidx{i})) '/' num2str(length(trieucNidx{i}))]);
    end
    save(fname,'tripucNidx','trieucNidx','-v7.3');
end

function checkNeuralQuadFeedforward(synTh, confTh)

    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralQuadFeedforward.mat'];
    if exist(fname,'file'), return; end
    outNids = {}; rcNids = {};
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralReciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % extract pure output neurons (no reciprocal)
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);
    poutNids = cell(nlen,1);
    for i=1:nlen
        rcnids = rcNids{i};
        outnids = outNids{i};
        if ~isempty(outnids)
            logis = ismember(outnids,rcnids);
            poutNids{i} = outnids(~logis); % pure out Nids
        end
    end

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    % extract quad feed forward neurons
    quadpffNids = cell(nlen,1); % pure output
    quadeffNids = cell(nlen,1); % extra (including reciprocal)
    parfor i=1:nlen
        poutnids = poutNids{i}; % pure output
        if isempty(poutnids), continue; end
        outnids = outNids{i}; 
        outlogi = ismember(tracedNids,outnids);

        tlogi = zeros(length(tracedNids),1,'logical'); % init zero logis
        tlogi2 = zeros(length(tracedNids),1,'logical'); % init zero logis
        for j=1:length(poutnids)
            logis = ismember(tracedNids,poutnids(j));
            poutnids2 = poutNids(logis); % pure output. this should extract one cell
            outnids2 = outNids(logis);   % all output (including reciprocal)
            if isempty(outnids2), continue; end
            plogi2 = ismember(tracedNids,poutnids2{1});
            alogi2 = ismember(tracedNids,outnids2{1});

            for k=j+1:length(poutnids)
                logis = ismember(tracedNids,poutnids(k));
                poutnids3 = poutNids(logis); % pure output. this should extract one cell
                outnids3 = outNids(logis);   % all output (including reciprocal)
                if isempty(outnids3), continue; end
                % find common target nids of poutnids(j) and (k)
                plogi3 = ismember(tracedNids,poutnids3{1});
                alogi3 = ismember(tracedNids,outnids3{1});
                tlogi = tlogi | (plogi2 & plogi3);
                tlogi2 = tlogi2 | (alogi2 & alogi3);
            end
        end
        quadpffNids{i} = tracedNids(tlogi & ~outlogi);
        quadeffNids{i} = tracedNids(tlogi2 & ~tlogi & ~outlogi);

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : quad feedforward process(' num2str(i) ') nid=' num2str(tracedNids(i)) ...
            ' out/pout/pff/eff=' num2str(length(outNids{i})) '/' num2str(length(poutNids{i})) '/' num2str(length(quadpffNids{i})) '/' num2str(length(quadeffNids{i}))]);
    end
    save(fname,'quadpffNids','quadeffNids','-v7.3');
end
