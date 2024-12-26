% make neural Struct Connectivity data.

function makeNeuralSC
    % check neural input & output voxels (FlyEM)
    checkNeuralInputOutputVoxels();

    % check neural input & output voxels (FlyWire)
    checkNeuralInputOutputVoxelsFw();
end

function checkNeuralInputOutputVoxels()
    hrateTh = 0.8; % FlyEM hemibrain synapse rate threshold

    fname = 'results/hemi_neuralInOutVoxels.mat';
    if exist(fname,'file'), return; end

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    Sdir = []; StoN = []; Srate = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;
    Sid = uint32(1:length(StoN));
    srate = (Srate >= hrateTh); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < hrateTh) = 0;  % use only accurate synapse more than 'rate'
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
    for i=1:length(tracedNids)
        logi = ismember(StoN,tracedNids(i)); % find synapses which belong to target neuron
        presids = Sid(logi & Sdir==1);
        postsids = Sid(logi & Sdir==2);
        logi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
        cpostsids = StoS(logi & ssrate & sstraced,2);
        [cpostsids, ia] = unique(cpostsids);
        disp(['process i=' num2str(i) ' postsids=' num2str(length(postsids)) ' presids=' num2str(length(presids))]);
    
        % get connected post-synapse counts from APL pre-synapses (output)
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

        % get post-synapse count of APL neuron (input)
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

function checkNeuralInputOutputVoxelsFw()
    wrateTh = 130; % FlyWire synapse score threshold

    fname = 'results/wire_neuralInOutVoxels.mat';
    if exist(fname,'file'), return; end

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load('data/flywire783_synapse.mat');
    score = (cleftScore >= wrateTh);
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
    for i=1:length(Nid)
        nidx = i;
        logi = ismember(preNidx,nidx); % find pre-synapses which belong to APL neurons
        presidx = Sidx(logi & valid & score);
        logi = ismember(postNidx,nidx); % find post-synapses which belong to APL neurons
        postsidx = Sidx(logi & valid & score);
        disp(['process i=' num2str(i) ' postsids=' num2str(length(postsidx)) ' presids=' num2str(length(presidx))]);
    
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

