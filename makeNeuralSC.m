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
    conf = getSCconfig('hemi', synTh, scTh);

    checkNeuralInputOutputVoxelsFw(conf);

%    checkNeuralInputOutputDistance(conf);

%    checkNeuralDBScan(conf, epsilon, minpts);

    checkNeuralReciprocalConnectionsFw(conf);

    checkDistanceReciprocalConnectionsFw(conf);
%{
    Cnames = {'SMP352'};
    Cnids = {[266187532, 266528078, 266528086, 296194535, 328593903, 5813009926]};
    checkNeuralNamedReciprocalConnections(Cnames, Cnids, synTh, scTh/100)
%}
    checkNeuralAutoConnectionsFw(conf);
%%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                checkNeuralTriFeedforwardFw(conf, (p12==1), (p13==1), (p23==1));
                checkDistanceTriConnectionsFw(conf,'Feedforward', (p12==1), (p13==1), (p23==1));
            end
        end
    end
%}
%%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                checkNeuralTriUnicycleFw(conf, (p12==1), (p13==1), (p23==1));
                checkDistanceTriConnectionsFw(conf,'Unicycle', (p12==1), (p13==1), (p23==1));
            end
        end
    end
%}
%    checkNeuralQuadFeedforward(synTh, scTh/100);

    % check neural input & output voxels (FlyWire)
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
%    scTh = 130; synTh = 5; % high confidence & connection setting.
    conf = getSCconfig('wire', synTh, scTh);

    checkNeuralInputOutputVoxelsFw(conf);

%    checkNeuralInputOutputDistance(conf);

%    checkNeuralDBScan(conf, epsilon, minpts);

    checkNeuralReciprocalConnectionsFw(conf);

    checkDistanceReciprocalConnectionsFw(conf);

    checkNeuralAutoConnectionsFw(conf);
%%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                checkNeuralTriFeedforwardFw(conf, (p12==1), (p13==1), (p23==1));
                checkDistanceTriConnectionsFw(conf,'Feedforward', (p12==1), (p13==1), (p23==1));
            end
        end
    end
%}
%%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                checkNeuralTriUnicycleFw(conf, (p12==1), (p13==1), (p23==1));
                checkDistanceTriConnectionsFw(conf,'Unicycle', (p12==1), (p13==1), (p23==1));
            end
        end
    end
%}
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

function checkNeuralInputOutputVoxelsFw(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat'];
    if exist(fname,'file'), return; end

    % FlyWire read neuron info
    Nid = [];
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    preNidx = []; postNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % read synapse location in FDA
    SpostlocFc = [];
    load(conf.sypostlocFdaFile);

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0;
    sz = size(Vt);

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    tracedNidx = 1:nlen;
    outIdx = cell(nlen,1);
    outCount = cell(nlen,1);
    inIdx = cell(nlen,1);
    inCount = cell(nlen,1);
    parfor i=1:nlen
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
        disp([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : process(' num2str(i) ') nid=' num2str(Nid(i)) ' postsids=' num2str(length(postsidx)) ' presids=' num2str(length(presidx))]);

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

function checkNeuralInputOutputDistance(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutDistance.mat'];
    if exist(fname,'file'), return; end

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % load input output voxel info
    niofname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat'];
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
        disp(['dist ' scname num2str(synTh) 'sr' num2str(scoreTh) ' : process i=' num2str(i) ' invox=' num2str(scinlen) ' outvox=' num2str(scoutlen)]);

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

function checkNeuralDBScan(conf, epsilon, minpts)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat'];
    if exist(fname,'file'), return; end

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % read neuron info (id, connection number, size)
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % load input output voxel info
    niofname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat'];
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
        disp(['dbscan ' conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : process (' num2str(i) ') nid=' num2str(Nid(i)) ' cls=' num2str(length(cls(cls>0))) ' invox=' num2str(scinlen) ' outvox=' num2str(scoutlen)]);
    end
    save(fname,'inlen','DBidx','-v7.3');
end

function checkNeuralReciprocalConnections(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
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
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
    load(fname);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize;
    tracedNids = Nid(Nstatus==1);
    Nidx = 1:length(tracedNids);

    % FlyEM read synapse info
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;

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

function checkDistanceReciprocalConnections(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalDistances.mat'];
    if exist(fname,'file'), return; end
    rcfname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
    rcNids = {}; rcpreSids = {}; rcpostSids = {};
    load(rcfname);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize;
    tracedNids = Nid(Nstatus==1);

    % FlyEM read synapse info
    Sloc = []; StoS = []; StoN = [];
    load('data/hemibrain_v1_2_synapses.mat');

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    nlen = length(rcNids);
    rcpreCloseDist = cell(nlen,1);
    rcpreCloseSids = cell(nlen,1);
    rcpostCloseDist = cell(nlen,1);
    rcpostCloseSids = cell(nlen,1);
    parfor i=1:nlen
        nids = rcNids{i};
        if isempty(nids), continue; end

        prelocs = single(Sloc(rcpreSids{i},:)); % hemibrain original space
        postlocs = single(Sloc(rcpostSids{i},:));
        preP = reshape(prelocs,[size(prelocs,1) 1 3]);
        postP = reshape(postlocs,[1 size(postlocs,1) 3]);
        preP = repmat(preP,[1 size(postP,2) 1]);
        postP = repmat(postP,[size(prelocs,1) 1 1]);
        D = sqrt(sum((preP-postP).^2,3)) * 8; % unit is nano meter

        % rcpre/rcpost synapse pair should be same reciprocal neurons.
        rcpresids = rcpreSids{i};
        rclen1 = length(rcpresids);
        rcnidlogis1 = cell(rclen1,1) ;
        for k=1:rclen1
            sslogi = (StoS(:,1)==rcpresids(k));
            postnids = StoN(StoS(sslogi,2));
            rcnidlogis1{k} = ismember(tracedNids,postnids);
        end
        rcpostsids = rcpostSids{i};
        rclen2 = length(rcpostsids);
        [sslogi,ssidx] = ismember(StoS(:,2),rcpostsids);
        idx = ssidx(sslogi);
        [~,idx2] = sort(idx);
        prenids2 = StoN(StoS(sslogi,1));
        prenids = prenids2(idx2); % sort to original order.
        L = zeros(rclen1,rclen2,'logical');
        for k=1:rclen1
            for l=1:rclen2
                L(k,l) = any(rcnidlogis1{k} & (tracedNids==prenids(l)));
            end
        end
        D(~L) = nan; % rcpre/rcpost synapse pair should be same reciprocal neurons.

        % find minimum distance of pre and post
        [minD, idxD] = min(D,[],2,'omitnan');
        rcpreCloseDist{i} = minD;
        rcpreCloseSids{i} = rcpostSids{i}(idxD);
        [minD, idxD] = min(D,[],1,'omitnan');
        rcpostCloseDist{i} = minD;
        rcpostCloseSids{i} = rcpreSids{i}(idxD);
        disp(['find closest reciprocal synapses (' num2str(i) ') prenid=' num2str(tracedNids(i))]);
%{
        [md,idx] = sort(minD);
        for j=1:length(idx)
            if j>3, break; end
            k = idx(j);
            pt = prelocs(k,:);
            logis = ismember(StoS(:,2),rcpreCloseSids{i}(k));
            presid = StoS(logis,1);
            disp(['      postnid=' num2str(StoN(presid)) ' preloc=' num2str(pt(1)) ',' num2str(pt(2)) ',' num2str(pt(3)) ...
                ' dist=' num2str(md(j)) 'nm (' num2str(md(j)/8) 'voxel)'])
        end
%}
    end
    save(fname,'rcpreCloseDist','rcpreCloseSids','rcpostCloseDist','rcpostCloseSids','-v7.3');
end

function checkNeuralReciprocalConnectionsFw(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat'];
    if exist(fname,'file'), return; end
    fnameNin = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];

    % FlyWire read neuron info
    Nid = [];
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    preNidx = []; postNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    inNidx = cell(nlen,1);
    outNidx = cell(nlen,1);
    rcNidx = cell(nlen,1);
    rcpreSidx = cell(nlen,1);
    rcpostSidx = cell(nlen,1);
    parfor i=1:nlen
        prlogi = ismember(preNidx,i); % find pre-synapses which belong to neuron(i)
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

        pologi = ismember(postNidx,i); % find post-synapses which belong to neuron(i)
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
            rcprlogi = ismember(preNidx,rcNidx{i}); % pre-synapse of reciprocal neurons
            rcpologi = ismember(postNidx,rcNidx{i}); % post-synapse of reciprocal neurons
            rcpreSidx{i} = Sidx(prlogi & rcpologi & valid & score);
            rcpostSidx{i} = Sidx(pologi & rcprlogi & valid & score); 
        end

        disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : reciprocal process(' num2str(i) ') nid=' num2str(Nid(i)) ' in/out/reci=' num2str(length(inNidx{i})) '/' num2str(length(outNidx{i})) '/' num2str(length(rcNidx{i})) ...
            ' rcpost/rcpre=' num2str(length(rcpostSidx{i})) '/' num2str(length(rcpreSidx{i}))]);
    end
    save(fname,'rcNidx','rcpostSidx','rcpreSidx','-v7.3');
    save(fnameNin,'inNidx','outNidx','-v7.3');
end

function checkDistanceReciprocalConnectionsFw(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalDistances.mat'];
    if exist(fname,'file'), return; end
    rcfname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat'];
    rcpreSidx = {}; rcpostSidx = {}; rcNidx = {};
    load(rcfname);

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % FlyWire read synapse locations
    Spostloc = []; Spreloc = [];
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    nlen = length(rcNidx);
    rcCloseDist = cell(nlen,1);
    rcCloseSidx = cell(nlen,1);
    for i=1:nlen
        nids = rcNidx{i};
        if isempty(nids), continue; end

        prelocs = single(Spreloc(rcpreSidx{i},:)) ./ conf.swcSize .* conf.voxelSize;    % reciprocal pre-synapse on Nid(i); FlyWire original space (unit is nano meter) (int32)
        postlocs = single(Spostloc(rcpostSidx{i},:)) ./ conf.swcSize .* conf.voxelSize; % reciprocal post-synapse on Nid(i); 
        preP = reshape(prelocs,[size(prelocs,1) 1 3]);
        postP = reshape(postlocs,[1 size(postlocs,1) 3]);
        preP = repmat(preP,[1 size(postP,2) 1]);
        postP = repmat(postP,[size(preP,1) 1 1]);
        D = sqrt(sum((preP-postP).^2,3));
        clear preP; clear postP;

        rcNidx1 = postNidx(rcpreSidx{i});
        rcNidx2 = preNidx(rcpostSidx{i})';
        rcNidx1 = repmat(rcNidx1,[1 size(rcNidx2,2)]);
        rcNidx2 = repmat(rcNidx2,[size(rcNidx1,1) 1]);
        D(rcNidx1~=rcNidx2) = nan; % rcpre/rcpost synapse pair should be same reciprocal neurons.
        clear rcNidx1; clear rcNidx2;

        % find minimum distance of pre and post (should be same reciprocal target neuron)
        if size(D,1) > size(D,2)
            [minD, idxD] = min(D,[],2,'omitnan');
            rcCloseDist{i} = minD;
            rcCloseSidx{i} = rcpostSidx{i}(idxD);
        else
            [minD2, idxD2] = min(D,[],1,'omitnan');
            rcCloseDist{i} = minD2;
            rcCloseSidx{i} = rcpreSidx{i}(idxD2);
        end
        disp(['find closest reciprocal synapses (' num2str(i) ') nid=' num2str(Nid(i))]);
    end
    save(fname,'rcCloseDist','rcCloseSidx','-v7.3');
end

function checkNeuralAutoConnections(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
    if ~exist(fname,'file'), return; end
    load(fname);
    if exist('autoNids','var'), return; end

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
    save(fname,'rcNids','rcpostSids','rcpreSids','autoNids','-v7.3');
end

function checkNeuralAutoConnectionsFw(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat'];
    if ~exist(fname,'file'), return; end
    load(fname);
    if exist('autoNidx','var'), return; end

    fnameNin = ['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(fnameNin);

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

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
    save(fname,'rcNidx','rcpostSidx','rcpreSidx','autoNidx','-v7.3');
end

function checkNeuralTriFeedforward(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_triangleFeedforward.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neural_Nin_Nout.mat'];
    outNids = {};
    load(fnameReci);
    load(fnameNin);

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

    % extract triangle feed forward neurons
    triffNids = cell(nlen,1);
    triffSids = cell(nlen,1);
    parfor i=1:nlen
        poutnids = poutNids{i}; % pure output
        if isempty(poutnids), continue; end

        % connected post-synapses
        logis = (StoN==tracedNids(i));
        sslogi = ismember(StoS(:,1),Sid(logis));
        cpostsids = StoS(sslogi & ssrate & sstraced,2); % get connected post-synapse from tracedNids(i)

        tricount = 0;
        ffNids = cell(length(poutnids),2); % pure output nids, extra nids (including reciprocal)
        ffSids = cell(length(poutnids),5); % (i)->tri(j) post-sids, pure (i)->ff(j) post-sids, pure tri(j)->ff(j), extra ..., extra ...
        for j=1:length(poutnids)
            logis = (tracedNids==poutnids(j));
            poutnids2 = poutNids(logis); % pure output. this should extract one cell
            outnids2 = outNids(logis);   % all output (including reciprocal)
%            idx = find(tracedNids==poutnids(j)); % find version.
%            poutnids2 = {poutNids{idx}};
            if ~isempty(outnids2)
                plogi = ismember(poutnids,poutnids2{1});
                pffnids = poutnids(plogi);
                alogi = ismember(poutnids,outnids2{1});
                effnids = poutnids(alogi & ~plogi);
                ffNids{j,1} = pffnids;
                ffNids{j,2} = effnids;
                if ~isempty(pffnids) || ~isempty(effnids)
                    % find synapses
                    poutlogi = ismember(StoN,poutnids(j)); % get post-synapse of tracedNids(i)->tri neuron(j)
                    cpslogi = ismember(cpostsids,Sid(poutlogi));
                    ffSids{j,1} = cpostsids(cpslogi);
                    sslogi = ismember(StoS(:,1),Sid(poutlogi));
                    tricpostsids = StoS(sslogi & ssrate & sstraced,2); % get connected post-synapse from tri neuron(j)

                    if ~isempty(pffnids)
                        poutlogi = ismember(StoN,pffnids); % get post-synapse of tracedNids(i)->pure ff neuron
                        cpslogi1 = ismember(cpostsids,Sid(poutlogi));
                        cpslogi2 = ismember(tricpostsids,Sid(poutlogi));
                        ffSids{j,2} = cpostsids(cpslogi1);
                        ffSids{j,3} = tricpostsids(cpslogi2);
                    end
                    if ~isempty(effnids)
                        poutlogi = ismember(StoN,effnids); % get post-synapse of tracedNids(i)->extra ff neuron
                        cpslogi1 = ismember(cpostsids,Sid(poutlogi));
                        cpslogi2 = ismember(tricpostsids,Sid(poutlogi));
                        ffSids{j,4} = cpostsids(cpslogi1);
                        ffSids{j,5} = tricpostsids(cpslogi2);
                    end
                    tricount = tricount + 1;
                end
            end
        end
        if tricount > 0
            triffNids{i} = ffNids;
            triffSids{i} = ffSids;
        end

        disp(['hemi' num2str(synTh) 'sr' num2str(confTh*100) ' : tri feedforward process(' num2str(i) ') nid=' num2str(tracedNids(i)) ...
            ' out/pout/tri=' num2str(length(outNids{i})) '/' num2str(length(poutNids{i})) '/' num2str(tricount)]);
    end
    save(fname,'triffNids','triffSids','-v7.3');
end

function checkNeuralTriFeedforwardFw(conf, isPure12, isPure13, isPure23)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;
    if isPure12, pstr='Pu'; else pstr='Rc'; end
    if isPure13, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end
    if isPure23, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangleFeedforward' pstr '.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    outNidx = {}; rcNidx = {};
    load(fnameReci);
    load(fnameNin);

    % FlyWire read neuron info
    Nid = []; 
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    postNidx = []; preNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

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

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    % extract triangle feed forward neurons
    triNidx = cell(nlen,1);
    triSidx = cell(nlen,1);
%    for i=1:500 %nlen
    parfor i=1:nlen
        if isPure12
            trinidx = poutNidx{i}; % pure output
        else
            trinidx = rcNidx{i};   % reciprocal output
        end
        if isempty(trinidx), continue; end

        % connected post-synapses
        cpostsidx = Sidx(preNidx==i & valid & score); % find connected post-synapses from Nid(i)

        tricount = 0;
        ffNidx = cell(length(trinidx),1); % ff nids
        ffSidx = cell(length(trinidx),3); % (i)->tri(j) sidx, (i)->ff(j), tri(j)->ff(j)
        for j=1:length(trinidx)           % triangle neuron loop
            logis = (Nidx==trinidx(j));
            if isPure23
                trioutnidx = poutNidx(logis); % tri output nidxs. this should extract one cell
            else
                trioutnidx = rcNidx(logis);   % tri output nidxs. this should extract one cell
            end
%            idx = find(Nidx==outnids(j)); % find version.
%            poutnidx2 = {poutNidx{idx}};
            if ~isempty(trioutnidx)
                if isPure13
                    onlogi = ismember(poutNidx{i},trioutnidx{1});
                    ffnidx = poutNidx{i}(onlogi);
                else
                    onlogi = ismember(rcNidx{i},trioutnidx{1});
                    ffnidx = rcNidx{i}(onlogi);
                end
                ffNidx{j,1} = ffnidx;
                if ~isempty(ffnidx)
                    % find synapses
                    trinsidx = Sidx(postNidx==trinidx(j) & valid & score);
                    cpslogi = ismember(cpostsidx,trinsidx);                   % get post-synapse of Nid(i)->tri neuron(j)
                    ffSidx{j,1} = cpostsidx(cpslogi);
                    tricpostsidx = Sidx(preNidx==trinidx(j) & valid & score); % get connected post-synapse from tri neuron(j)

                    fflogi = ismember(postNidx, ffnidx);
                    ffinsidx = Sidx(fflogi & valid & score);
                    cpslogi1 = ismember(cpostsidx,ffinsidx);    % get post-synapse of Nid(i)->ff neuron
                    cpslogi2 = ismember(tricpostsidx,ffinsidx); % get post-synapse of tri(j)->ff neuron
                    ffSidx{j,2} = cpostsidx(cpslogi1);
                    ffSidx{j,3} = tricpostsidx(cpslogi2);
                    tricount = tricount + 1;
                end
            end
        end
        if tricount > 0
            triNidx{i} = ffNidx;
            triSidx{i} = ffSidx;
        end

        disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : tri feedforward ' pstr ' process(' num2str(i) ') nid=' num2str(Nid(i)) ...
            ' out/pout/tri=' num2str(length(outNidx{i})) '/' num2str(length(poutNidx{i})) '/' num2str(tricount)]);
    end
    save(fname,'triNidx','triSidx','-v7.3');
end

function checkDistanceTriConnectionsFw(conf, trType, isPure12, isPure13, isPure23)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;
    if isPure12, pstr='Pu'; else pstr='Rc'; end
    if isPure13, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end
    if isPure23, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangle' trType pstr 'Dists.mat'];
    if exist(fname,'file'), return; end
    fnameTri = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangle' trType pstr '.mat'];
    fnameNin = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    triNidx = {}; triSidx = {};
    load(fnameTri);
    load(fnameNin);

    % FlyWire read neuron info
    Nid = [];
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    postNidx = []; preNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % FlyWire read synapse locations
    Spostloc = []; Spreloc = [];
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    nlen = length(triNidx);
    triCloseDist = cell(nlen,1);
    triCloseSidx = cell(nlen,1);
%    for i=1:500
    parfor i=1:nlen
        trinidx = triNidx{i}; % tri nidx
        if isempty(trinidx), continue; end
        trisidx = triSidx{i};

        D = cell(size(trinidx,1),1);
        S = cell(size(trinidx,1),1);
        for j=1:size(trinidx,1)
            ffnidx = trinidx{j};
            if ~isempty(ffnidx)
                sylocs1 = single(Spreloc(trisidx{j,1},:)) ./ conf.swcSize .* conf.voxelSize;      % tri pre-synapse on Nid(i)
                switch(trType)
                case 'Feedforward'
                    sylocs2 = single(Spreloc(trisidx{j,2},:)) ./ conf.swcSize .* conf.voxelSize;  % ff pre-synapse on Nid(i); 
                case 'Unicycle'
                    sylocs2 = single(Spostloc(trisidx{j,2},:)) ./ conf.swcSize .* conf.voxelSize; % uc post-synapse on Nid(i); 
                end
                sylocs3 = single(Spreloc(trisidx{j,3},:)) ./ conf.swcSize .* conf.voxelSize;      % ff pre-synapse on Tri(j)
                syP1 = reshape(sylocs1,[size(sylocs1,1) 1 3]);
                syP2 = reshape(sylocs2,[1 size(sylocs2,1) 3]);
                syP3 = reshape(sylocs3,[size(sylocs3,1) 1 3]);
                syP1 = repmat(syP1,[1 size(syP2,2) 1]);
                syP12 = repmat(syP2,[size(syP1,1) 1 1]);
                D1 = sqrt(sum((syP1-syP12).^2,3));
                syP1 = []; syP12 = []; % clear memory

                % Nid(i), Tri(j) are fixed. ff/uc are not.
                if size(D1,1) > size(D1,2)
                    [minD, idxD] = min(D1,[],2,'omitnan');
                    triCdist12 = minD;
                    triCsidx1  = trisidx{j,1};
                    triCsidx2a = trisidx{j,2}(idxD);
                    [minD2, idxD2] = min(D1,[],1,'omitnan');
                    triCdist12c = minD2;
                    triCsidx1c  = trisidx{j,1}(idxD2);
                    triCsidx2c  = trisidx{j,2};
                else
                    [minD2, idxD2] = min(D1,[],1,'omitnan');
                    triCdist12 = minD2;
                    triCsidx1  = trisidx{j,1}(idxD2);
                    triCsidx2a = trisidx{j,2};
                end
                D1 = []; % clear memory

                syP23 = repmat(syP2,[size(syP3,1) 1 1]);
                syP3 = repmat(syP3,[1 size(syP2,2) 1]);
                D2 = sqrt(sum((syP3-syP23).^2,3));
                syP23 = []; syP3 = []; % clear memory

                switch(trType)
                case 'Feedforward'
                    syN23 = postNidx(trisidx{j,2})';
                case 'Unicycle'
                    syN23 = preNidx(trisidx{j,2})';
                end                   
                syN3 = postNidx(trisidx{j,3});
                syN23 = repmat(syN23,[size(syN3,1) 1]);
                syN3 = repmat(syN3,[1 size(syN23,2)]);
                D2(syN3~=syN23) = nan; % post synapse pair should be same neuron.
                syN23 = []; syN3 = [];

                % Nid(i), Tri(j) are fixed. ff/uc are not.
                if size(D2,1) > size(D2,2)
                    [minD, idxD] = min(D2,[],2,'omitnan');
                    triCdist23 = minD;
                    triCsidx2b = trisidx{j,2}(idxD);
                    triCsidx3  = trisidx{j,3};
                    [minD2, idxD2] = min(D2,[],1,'omitnan');
                    triCdist23c = minD2;
                    triCsidx2c  = trisidx{j,2};
                    triCsidx3c  = trisidx{j,3}(idxD2);
                else
                    [minD2, idxD2] = min(D2,[],1,'omitnan');
                    triCdist23 = minD2;
                    triCsidx2b = trisidx{j,2};
                    triCsidx3  = trisidx{j,3}(idxD2);
                end
                D2 = []; % clear memory

                if size(triCdist12,1) > size(triCdist12,2)
                    if size(triCdist23,1) > size(triCdist23,2)
                        % checking both side vector
                        Dj = zeros(length(triCdist12),1,'single');
                        Sj = zeros(length(triCdist12),3,'int32');
                        for k=1:length(triCdist12)
                            logis = (triCsidx2a(k)==triCsidx2c);
                            Dj(k) = triCdist12(k) + triCdist23c(logis);
                            Sj(k,:) = [triCsidx1(k) triCsidx2c(logis) triCsidx3c(logis)];
                        end
                        Dj2 = zeros(length(triCdist23),1,'single');
                        Sj2 = zeros(length(triCdist23),3,'int32');
                        for k=1:length(triCdist23)
                            logis = (triCsidx2c==triCsidx2b(k));
                            Dj2(k) = triCdist12c(logis) + triCdist23(k);
                            Sj2(k,:) = [triCsidx1c(logis) triCsidx2c(logis) triCsidx3(k)];
                        end
                        D{j} = [Dj; Dj2];
                        S{j} = [Sj; Sj2];
                    else
                        Dj = zeros(length(triCdist12),1,'single');
                        Sj = zeros(length(triCdist12),3,'int32');
                        for k=1:length(triCdist12)
                            logis = (triCsidx2a(k)==triCsidx2b);
                            Dj(k) = triCdist12(k) + triCdist23(logis);
                            Sj(k,:) = [triCsidx1(k) triCsidx2b(logis) triCsidx3(logis)];
                        end
                        D{j} = Dj;
                        S{j} = Sj;
                    end
                else
                    % find close tri(j)-ff/uc synapse pairs
                    if size(triCdist23,1) > size(triCdist23,2)
                        Dj = zeros(length(triCdist23),1,'single');
                        Sj = zeros(length(triCdist23),3,'int32');
                        for k=1:length(triCdist23)
                            logis = (triCsidx2a==triCsidx2b(k));
                            Dj(k) = triCdist12(logis) + triCdist23(k);
                            Sj(k,:) = [triCsidx1(logis) triCsidx2a(logis) triCsidx3(k)];
                        end
                        D{j} = Dj;
                        S{j} = Sj;
                    else
                        D{j} = triCdist12(:) + triCdist23(:);
                        if sum(triCsidx2a-triCsidx2b) ~= 0, disp([num2str(i) '-' num2str(j) ' 2a == 2b error']); end
                        S{j} = [triCsidx1(:), triCsidx2a(:), triCsidx3(:)]; % 2a == 2b
                    end
                end
                % check neural connections.
                for k=1:size(S{j},1)
                    sidx = S{j}(k,:);
                    prnidx = preNidx(sidx);
                    nidx = postNidx(sidx);
                    if strcmp(trType,'Unicycle')
                        tmp = nidx(2); nidx(2) = prnidx(2); prnidx(2) = tmp;
                    end
                    if prnidx(1) ~= i || prnidx(2) ~= i || prnidx(3) ~= nidx(1) || nidx(2) ~= nidx(3)
                        disp(['bad ' num2str(i) '-' num2str(j) '-' num2str(k) ' and nidx=' num2str(i) ',' num2str(nidx(2)) ',' num2str(nidx(3))]);
                    end
                end
            end
        end
        triCloseDist{i} = D;
        triCloseSidx{i} = S;
        disp(['find closest triangle synapses (' num2str(i) ') nid=' num2str(Nid(i))]);
    end
    save(fname,'triCloseDist','triCloseSidx','-v7.3');
end

function checkNeuralTriUnicycle(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_triangleUnicycle.mat'];
    if exist(fname,'file'), return; end
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
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

function checkNeuralTriUnicycleFw(conf, isPure12, isPure13, isPure23)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;
    if isPure12, pstr='Pu'; else pstr='Rc'; end
    if isPure13, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end
    if isPure23, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangleUnicycle' pstr '.mat'];
    if exist(fname,'file'), return; end
    outNidx = {}; inNidx = {}; rcNidx = {};
    fnameReci = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat'];
    fnameNin = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(fnameReci);
    load(fnameNin);

    % FlyWire read neuron info
    Nid = []; 
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    postNidx = []; preNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

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

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    triNidx = cell(nlen,1);
    triSidx = cell(nlen,1);
%    for i=1:500 %nlen
    parfor i=1:nlen
        if isPure12
            trinidx = poutNidx{i}; % pure output
        else
            trinidx = rcNidx{i};   % reciprocal output
        end
        if isempty(trinidx), continue; end

        % connected post-synapses
        postsidx = Sidx(postNidx==i & valid & score); % find post-synapses on Nid(i)
        cpostsidx = Sidx(preNidx==i & valid & score); % find connected post-synapses from Nid(i)

        tricount = 0;
        ucNidx = cell(length(trinidx),1); % uc nids
        ucSidx = cell(length(trinidx),3); % (i)->tri(j) sidx, (i)->ff(j), tri(j)->ff(j)
        for j=1:length(trinidx)           % triangle neuron loop
            logis = (Nidx==trinidx(j));
            if isPure23
                trioutnidx = poutNidx(logis); % tri output nidxs. this should extract one cell
            else
                trioutnidx = rcNidx(logis);   % tri output nidxs. this should extract one cell
            end
%            idx = find(Nidx==poutnids(j)); % find version.
%            poutnidx2 = {poutNidx{idx}};
            if ~isempty(trioutnidx)
                if isPure13
                    onlogi = ismember(pinNidx{i},trioutnidx{1});
                    ucnidx = pinNidx{i}(onlogi);
                else
                    onlogi = ismember(rcNidx{i},trioutnidx{1});
                    ucnidx = rcNidx{i}(onlogi);
                end
                ucNidx{j,1} = ucnidx;
                if ~isempty(ucnidx)
                    % find synapses
                    trinsidx = Sidx(postNidx==trinidx(j) & valid & score);    % get post-synapse of Nid(i)->tri neuron(j)
                    cpslogi = ismember(cpostsidx,trinsidx);
                    ucSidx{j,1} = cpostsidx(cpslogi);
                    tricpostsidx = Sidx(preNidx==trinidx(j) & valid & score); % get connected post-synapse from tri neuron(j)

                    uclogi = ismember(postNidx, ucnidx);
                    ucinsidx = Sidx(uclogi & valid & score);
                    uclogi2 = ismember(preNidx, ucnidx);
                    ucoutsidx = Sidx(uclogi2 & valid & score);
                    cpslogi1 = ismember(postsidx,ucoutsidx);    % get pre-synapse of Nid(i)<-uc neuron
                    cpslogi2 = ismember(tricpostsidx,ucinsidx); % get post-synapse of tri(j)->uc neuron
                    ucSidx{j,2} = postsidx(cpslogi1);
                    ucSidx{j,3} = tricpostsidx(cpslogi2);
                    tricount = tricount + 1;
                end
            end
        end
        if tricount > 0
            triNidx{i} = ucNidx;
            triSidx{i} = ucSidx;
        end

        disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : tri unicycle process(' num2str(i) ') nid=' num2str(Nid(i)) ...
            ' in/pin/tri=' num2str(length(inNidx{i})) '/' num2str(length(pinNidx{i})) '/' num2str(tricount)]);
    end
    save(fname,'triNidx','triSidx','-v7.3');
end

function checkNeuralQuadFeedforward(synTh, confTh)
    fname = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_neuralQuadFeedforward.mat'];
    if exist(fname,'file'), return; end
    outNids = {}; rcNids = {};
    fnameReci = ['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh*100) '_reciprocalConnections.mat'];
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
