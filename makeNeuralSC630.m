% make neural Struct Connectivity data.

function makeNeuralSC630
    % check neural input & output voxels (FlyWire)
    scTh = 50; synTh = 0; % for checking flywire codex compatible
    conf = getSCconfig('wire630', synTh, scTh);

    checkNeuralReciprocalConnectionsFw(conf);

    checkNeuralNetworkPropertiesFw(conf);
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

function checkNeuralNetworkPropertiesFw(conf)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;

    fname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_network_prop.mat'];
    if exist(fname,'file'), return; end
    rcfname = ['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat'];
    load(rcfname);

    % FlyWire read neuron info
    Nid = [];
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    preNidx = []; postNidx = [];
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % count pre (to post) and post synapse count
    nlen = length(Nid);
    S = zeros(nlen,nlen,'logical');
    for i=1:nlen
        outnidx = outNidx{i};
        innidx = inNidx{i};
        S(i,outnidx) = 1;
        S(innidx,i) = 1;
    end
    L = sum(S,'all'); % total number of links
    dens = L / (nlen*(nlen-1));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : ' num2str(nlen) ' neurons, ' num2str(L) ' connections, density=' num2str(dens)]);

    E=sparse(S);
    save('d:\work\gs\erdir2.mat','E');

    % reciprocity
    [R, gR, count] = calcReciprocity(S);
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : reci ' num2str(count) ' neurons, global reciprocity=' num2str(gR)]);

    % clustering coefficient
    [gC, tri, alltri] = calcGlobalClusteringCoeff(sparse(S));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : tri ' num2str(count) ' neurons, global clustering coeff=' num2str(gC) ', tri=' num2str(tri) ', alltri=' num2str(alltri)]);

    % compared with ER
    E = generateERgraph(nlen, dens);
    EL = sum(E,'all'); % total number of links
    Edens = EL / (nlen*(nlen-1));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : ER ' num2str(nlen) ' neurons, ' num2str(EL) ' connections, density=' num2str(Edens)]);

    [R, EgR, count] = calcReciprocity(E);
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : ER reci ' num2str(count) ' neurons, global reciprocity=' num2str(EgR) ', x ER=' num2str(gR/EgR)]);

    [EgC, tri, alltri] = calcGlobalClusteringCoeff(sparse(E));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : ER tri ' num2str(count) ' neurons, global clustering coeff=' num2str(EgC) ', x ER=' num2str(gC/EgC)]);

    % compared with CFG
    G = generateCFGgraph(S);
    GL = sum(G,'all'); % total number of links
    Edens = GL / (nlen*(nlen-1));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : CFG ' num2str(nlen) ' neurons, ' num2str(GL) ' connections, density=' num2str(Edens)]);

    [R, GgR, count] = calcReciprocity(G);
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : CFG reci ' num2str(count) ' neurons, global reciprocity=' num2str(GgR) ', x CFG=' num2str(gR/GgR)]);

    [GgC, tri, alltri] = calcGlobalClusteringCoeff(sparse(G));
    disp([scname num2str(synTh) 'sr' num2str(scoreTh) ' : CFG tri ' num2str(count) ' neurons, global clustering coeff=' num2str(GgC) ', x CFG=' num2str(gC/GgC)]);
end

