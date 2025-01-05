% analyze neural functional connectivity matrix (1st and 2nd level analysis).
% this script should run after makeNeuralSC.m, extractNeuralTimeseries.m first.

function analyzeNeuralFC
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};

    % using subjects (flys). sbj 7 shows NaN row in FC matrix
    sbjids = [1 2 3 4 5 6 8 9];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % analyze neural functional connectivity (FlyEM)
    confTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    confTh = 60; synTh = 5; % almost flywire codex compatible setting
%    analyzeNeuralFc(preproc, hpfTh, smooth, nuisance, sbjids, 'hemi', synTh, scTh); % FlyEM hemi brain.

    % checkin named neuron fiber. SMP136, etc.
%    Cnames = {'SMP136','SMP352'};
%    Cnids = {[359292188, 360673271, 423563311, 546325197],[266187532, 266528078, 266528086, 296194535, 328593903, 5813009926]};
    Cnames = {'SMP352','SMP352reci'};
    Cnids = {[266187532, 266528078, 266528086, 296194535, 328593903, 5813009926],...
        [298262644,328943204,452029745,5813019513,5813020684]};
    checkingNamedFiberSynapses(Cnames, Cnids, 'hemi', synTh, confTh)

    % checking fiber neuron & synapse
    checkingFiberNeurons(preproc, hpfTh, smooth, nuisance, sbjids, 'hemi', synTh, confTh); % FlyEM hemi brain.

    checkingFiberSynapses('hemi', synTh, confTh);

    % analyze neural functional connectivity (FlyWire)
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
    analyzeNeuralFc(preproc, hpfTh, smooth, nuisance, sbjids, 'wire', synTh, scTh); % FlyWire whole brain.

end

function analyzeNeuralFc(preproc, hpfTh, smooth, nuisance, sbjids, scname, synTh, confTh)
    % load whole brain mask
    mV = int32(niftiread('template/thresholded_FDACal_mask.nii.gz')); % mask should have same transform with 4D nifti data
    midx = find(mV>0); mV(:) = 0;
    mV(midx) = 1:length(midx);

    % load neural SC result (from makeNeuralSC.m)
    inIdx = {}; outIdx = {};
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(confTh) '_neuralInOutVoxels.mat']);

    if ~exist('results/neuralfc','dir'), mkdir('results/neuralfc'); end

    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                % output file
                outfname = ['results/neuralfc/' smooth{k} hpfstr nuisance{n} preproc scname num2str(synTh) 'sr' num2str(confTh) '-fc.mat'];
                if exist(outfname,'file'), continue; end

                % load ROI time-series (from extractNeuralTimeseries.m)
                listing = dir(['results/neuralts/' smooth{k} hpfstr nuisance{n} preproc 'sub-fly-*-ts.mat']);
                CX = {};
                for i=1:length(sbjids)
                    name = listing(sbjids(i)).name;
                    disp(['loading ' name ' ...']);
                    t = load(['results/neuralts/' name]);
                    CX{i} = t.X;
                    clear t;
                end
                cmlen = length(CX);

                % calc FC in each neuron
                mFz = cell(length(inIdx),1);
                for i=1:length(inIdx)
%                parfor i=1:length(inIdx)  % CX is too big
                    scinidx = inIdx{i};
                    scoutidx = outIdx{i};
                    if length(scinidx)==0 || length(scoutidx)==0
                        disp(['empty input or output voxel i=' num2str(i)]);
                        continue;
                    end
                    minidx = mV(scinidx); minidx(minidx==0)=[];
                    moutidx = mV(scoutidx); moutidx(moutidx==0)=[];
                    if length(minidx)==0 || length(moutidx)==0
                        disp(['empty input or output voxel (out of mask) i=' num2str(i)]);
                        continue;
                    end

                    disp(['process i=' num2str(i)]);
                    F3 = zeros(length(minidx),length(moutidx),cmlen,'single');
                    for j=1:cmlen
                        F3(:,:,j) = corr(CX{j}(minidx,:)',CX{j}(moutidx,:)');
                    end
                    F3z = atanh(F3); % z transformed (better FC-SC corr).
                    mFz{i} = nanmean(F3z,3);
%                    figure; imagesc(mFz{i}); colorbar;
                end
                save(outfname, 'mFz', '-v7.3');
            end
        end
    end
end

function checkingFiberNeurons(preproc, hpfTh, smooth, nuisance, sbjids, scname, synTh, confTh)
    roiNum = 1000;
    distPer = 70; % distance threshold (percent)

    idstr = [scname 'DistKm' num2str(roiNum)];
    if synTh==0 && confTh==80
        roitype = idstr;
    else
        roitype = [idstr '_hb' num2str(synTh) 'sr' num2str(confTh)];
    end

    % load ROI SC info
    fname = ['results/sc/' lower(roitype) '_connectlist.mat']; % default setting
    load(fname);

    % load ROI distance matrix
    fname = ['results/sc/' lower(idstr) '_dist.mat'];
    load(fname);
%    figure; histogram(distMat); title('histogram of ROI distances')
    distTh = prctile(distMat(:),distPer);

    % load ROI neuron input and output
    fname = ['results/cache/' lower(roitype) '_Nin_Nout.mat']; % default setting
    load(fname);

    % load ROI FC info
    hpfstr = '';
    if hpfTh(1) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
    T2 = [];
    fname = ['results/fc/' smooth{1} hpfstr nuisance{1} preproc roitype '-func.mat'];
    load(fname);

    % load ROI time-series (from extractROItimeseries.m)
    str = split(roitype,'_');
    load(['results/ts/' smooth{1} hpfstr nuisance{1} preproc str{1} '-ts.mat']); % FlyWire SC calculation.
    ids = primaryIds;
    F3 = [];
    for i=1:length(sbjids)
        CM = single(corr(CX{sbjids(i)}'));
        D = CM(ids,ids);
        F3 = cat(3,F3,D);
    end
    F3z = atanh(F3); % z transformed (better FC-SC corr).
    mFz = nanmean(F3z,3);
    clear CX;

    % find long distance and connections
    ldcfname = ['results/sc/' lower(roitype) '_longDist' num2str(distPer) '_neurons.mat'];
    if exist(ldcfname,'file')
        load(ldcfname)
    else
        Nio = cell(roiNum,roiNum);
        Nii = cell(roiNum,roiNum);
        Noo = cell(roiNum,roiNum);
        for i=1:roiNum
            disp(['long dist process i=' num2str(i) '/' num2str(roiNum)]);
            outnids1 = Nout{i}{2};
            innids1 = Nin{i}{2};
            for j=1:roiNum
                if distMat(i,j) < distTh, continue; end
                outnids2 = Nout{j}{2};
                innids2 = Nin{j}{2};
                iologi = ismember(outnids1,innids2); 
                Nio{i,j} = outnids1(iologi); % neurons from ROI(i) output to ROI(j) input
    
                % checking in-in, out-out
                if j<=i, continue; end
                iilogi = ismember(innids1,innids2);
                Nii{i,j} = innids1(iilogi);
                oologi = ismember(outnids1,outnids2);
                Noo{i,j} = outnids1(oologi);
            end
        end
        save(ldcfname,'Nio','Nii','Noo','-v7.3');
    end

    t2cfname = ['results/cache/' lower(roitype) '_longDist' num2str(distPer) '_T2.mat'];
    if ~exist(t2cfname,'file') % T2 cache for later
        dM = distMat;
        dM(dM<distTh)=0; dM(dM>0)=1;
        save(t2cfname,'mFz','T2','dM','-v7.3');
    end

    Cio = zeros(roiNum,roiNum);
    Cii = zeros(roiNum,roiNum);
    Coo = zeros(roiNum,roiNum);
    for i=1:roiNum
        for j=1:roiNum
            Cio(i,j) = numel(Nio{i,j});
            Cii(i,j) = numel(Nii{i,j});
            Coo(i,j) = numel(Noo{i,j});
        end
    end
    ids = primaryIds;
    figure; imagesc(log10(Cio(ids,ids))); colorbar; daspect([1 1 1]); title([roitype ' long distance neuron count in/out']); xlabel('input ROI'); ylabel('output ROI');
    figure; imagesc(log10(Cii(ids,ids))); colorbar; daspect([1 1 1]); title([roitype ' long distance neuron count in/in']); xlabel('input ROI'); ylabel('input ROI');
    figure; imagesc(log10(Coo(ids,ids))); colorbar; daspect([1 1 1]); title([roitype ' long distance neuron count out/out']); xlabel('output ROI'); ylabel('output ROI');

    CC = {Cio,Coo}; Ctype={'in/out','out/out'};
    for k=1
        logis = (distMat>=distTh);
        T2d = T2(logis); mFzd = mFz(logis); ldN = CC{k}(logis);
        T2d0 = T2d(ldN==0); T2d1 = T2d(ldN>0); ldN1 = ldN(ldN>0);
        mFzd0 = mFzd(ldN==0); mFzd1 = mFzd(ldN>0);
        [h, pt] = ttest2(T2d0,T2d1);
        [h, pz] = ttest2(mFzd0,mFzd1);

%        figure; histogram(T2d0); hold on; histogram(T2d1); hold off; title([roitype ' ' Ctype{k} ' long distance T-val fiber vs. no-fiber p=' num2str(pt)]);
        X = nan(max(numel(T2d0),numel(T2d1)),2); X(1:numel(T2d0),1) = T2d0; X(1:numel(T2d1),2) = T2d1;
        figure; boxplot(X); title([roitype ' ' Ctype{k} ' long distance T-val fiber vs. no-fiber p=' num2str(pt)]); ylabel('T-value'); xticklabels({'no-fiber','fiber'})
        [r, pr] = corr(T2d1,ldN1);
        figure; scatter(ldN1,T2d1,12,'x'); title([roitype ' ' Ctype{k} ' long distance T-val vs. neuron count r=' num2str(r)]); ylabel('T-value'); xlabel('neuron count');
%{
        X = nan(max(numel(mFzd0),numel(mFzd1)),2); X(1:numel(mFzd0),1) = mFzd0; X(1:numel(mFzd1),2) = mFzd1;
        figure; boxplot(X); title([roitype ' ' Ctype{k} ' long distance mFC(z) fiber vs. no-fiber p=' num2str(pz)]); ylabel('m-FC(z)'); xticklabels({'no-fiber','fiber'})
        [r, pr] = corr(mFzd1,ldN1);
        figure; scatter(ldN1,mFzd1,12,'x'); title([roitype ' ' Ctype{k} ' long distance mFC(z) vs. neuron count r=' num2str(r)]); ylabel('m-FC(z)'); xlabel('neuron count');
%}    
        G = {[1 10],[10 20],[20 40],[40 60],[60 80],[80 100],[100 200]};
        Glabels = {};
        GT2 = [];
        for i=1:length(G)
            X = T2d1(G{i}(1)<= ldN1 & ldN1 < G{i}(2));
%            X = mFzd1(G{i}(1)<= ldN1 & ldN1 < G{i}(2));
            if i>1, GT2(:,i) = nan; end
            GT2(1:numel(X),i) = X;
            Glabels{i} = [num2str(G{i}(1)) '-' num2str(G{i}(2))];
        end
        figure; boxplot(GT2); title([roitype ' ' Ctype{k} ' long distance T-val vs. neuron count r=' num2str(r)]); ylabel('T-value'); xlabel('neuron count'); xticklabels(Glabels);
    end
end

function checkingNamedFiberSynapses(Cnames, Cnids, scname, synTh, confTh)

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    StoN = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    Sid = uint32(1:length(StoN))';
    srate = (Srate >= confTh/100); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < confTh/100) = 0;  % use only accurate synapse more than 'rate'
    clear Srate;

    % FlyEM read synapse location in FDA
    load('data/hemibrain_v1_2_synapseloc_fdacal.mat');

    for i=1:length(Cnames)
        % fiber neuron synapse matrix
        symfname = ['results/fb/' scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} '_synMat.mat'];
        if exist(symfname,'file')
            load(symfname)
        else
            nids = Cnids{i};
            [fbsyMat, fbsyC] = calcFiberSynapseMat(nids, StoN, Sdir, StoS, ssrate, sstraced);
            save(symfname,'nids','fbsyMat','fbsyC','-v7.3');
        end
        nidNum = length(Cnids{i});
        figure; imagesc(fbsyMat,[0 4]); colorbar; daspect([1 1 1]); title([scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} ' neurons synapse count']);
        ylabel('output neuron'); xlabel('input neuron'); xticks(1:nidNum); yticks(1:nidNum); xticklabels(nidLabels(Cnids{i})); yticklabels(nidLabels(Cnids{i}));

        % output post-synapse cloud of neuron fibers for confirmation
        niifname = ['results/fb/' scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} '_postsynFDACal.nii'];
        if ~exist([niifname '.gz'],'file')
            [hbV, hbinfo] = calcFiberSynapseV(fbsyC, SlocFc);
            niftiwrite(hbV,niifname,hbinfo,'Compressed',true);
        end

        % find common input and output neuron & synapses
        symfname = ['results/fb/' scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} '_commonInOut.mat'];
        if exist(symfname,'file')
            load(symfname)
        else
            [innids, outnids, infbsyMat, outfbsyMat, Cinfbsy, Coutfbsy] = findFiberCommonInOut(Cnids{i}, Nid, StoN, Sdir, StoS, ssrate, sstraced);
            save(symfname,'innids','outnids','infbsyMat','outfbsyMat','Cinfbsy','Coutfbsy','-v7.3');
            if ~isempty(innids)
                csvfname = ['results/fb/' scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} '_commonInNeuron.csv'];
                writematrix(innids',csvfname);
            end
            if ~isempty(outnids)
                csvfname = ['results/fb/' scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} '_commonOutNeuron.csv'];
                writematrix(outnids',csvfname);
            end
        end
        figure; imagesc(infbsyMat,[0 4]); colorbar; daspect([1 1 1]); title([scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} ' neurons synapse count']);
        ylabel('out to fiber neurons'); xlabel('fiber neurons'); xticks(1:length(Cnids{i})); yticks(1:length(innids)); xticklabels(nidLabels(Cnids{i})); yticklabels(nidLabels(innids));
        figure; imagesc(outfbsyMat,[0 4]); colorbar; daspect([1 1 1]); title([scname num2str(synTh) 'sr' num2str(confTh) '_' Cnames{i} ' neurons synapse count']);
        ylabel('fiber neurons'); xlabel('out from fiber neurons'); xticks(1:length(outnids)); yticks(1:length(Cnids{i})); xticklabels(nidLabels(outnids)); yticklabels(nidLabels(Cnids{i}));
    end
end

function labels = nidLabels(nids)
    labels = cell(1,length(nids));
    for i=1:length(nids)
        labels{i} = [num2str(nids(i))];
    end
end

function checkingFiberSynapses(scname, synTh, confTh)
    roiNum = 1000;
    distPer = 70; % distance threshold (percent)

    idstr = [scname 'DistKm' num2str(roiNum)];
    if synTh==0 && confTh==80
        roitype = idstr;
    else
        roitype = [idstr '_hb' num2str(synTh) 'sr' num2str(confTh)];
    end

    if ~exist('results/fb','dir'), mkdir('results/fb'); end

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    StoN = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    Sid = uint32(1:length(StoN))';
    srate = (Srate >= confTh/100); % use only accurate synapse more than 'rate'
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < confTh/100) = 0;  % use only accurate synapse more than 'rate'
    clear Srate;

    % FlyEM read synapse location in FDA
    load('data/hemibrain_v1_2_synapseloc_fdacal.mat');

    % load ROI atlas
    atlasinfo = niftiinfo(['atlas/' roitype 'atlasCal.nii.gz']);
    aV = niftiread(atlasinfo);

    % load long distance neurons
    ldcfname = ['results/sc/' lower(roitype) '_longDist' num2str(distPer) '_neurons.mat'];
    load(ldcfname)

    Cio = zeros(roiNum,roiNum);
    Cii = zeros(roiNum,roiNum);
    Coo = zeros(roiNum,roiNum);
    for i=1:roiNum
        for j=1:roiNum
            Cio(i,j) = numel(Nio{i,j});
            Cii(i,j) = numel(Nii{i,j});
            Coo(i,j) = numel(Noo{i,j});
        end
    end

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(24);

    % let's check top 50 count fibers
    [cios, cioidx] = sort(Cio(:),'descend');
    cnt = zeros(roiNum,2);
    for ii=1:200
        [x, y] = ind2sub([roiNum roiNum], cioidx(ii));
        nids = Nio{x,y}; % neurons from ROI(x) to ROI(y)
        cnt(x,1) = cnt(x,1) + 1;
        cnt(y,2) = cnt(y,2) + 1;
        if cnt(x,1) > 2 || cnt(y,2) > 2 % ignore (too many) redundunt regions
            disp(['count over. ignore ROI pair fb' num2str(x) '-' num2str(y)]);
            continue;
        end
        disp(['process(' num2str(ii) ') ROI pair fb' num2str(x) '-' num2str(y)]);

        % output ROI info for confirmation
        niifname = ['results/fb/' idstr '_fb' num2str(x) '-' num2str(y) 'atlasCal.nii'];
        if ~exist([niifname '.gz'],'file')
            V=aV; V(:)=0; V(aV==x)=1; V(aV==y)=2; 
            niftiwrite(V,niifname,atlasinfo,'Compressed',true);
        end

        % output neuron list by csv
        csvfname = ['results/fb/' idstr '_fb' num2str(x) '-' num2str(y) '_neurons.csv'];
        if ~exist(csvfname,'file')
            writematrix(nids',csvfname);
        end

        % fiber neuron synapse matrix
        symfname = ['results/sc/' lower(roitype) '_fb' num2str(x) '-' num2str(y) '_synMat.mat'];
        if exist(symfname,'file')
            load(symfname)
        else
            [fbsyMat, fbsyC] = calcFiberSynapseMat(nids, StoN, Sdir, StoS, ssrate, sstraced);
            save(symfname,'nids','fbsyMat','fbsyC','-v7.3');
        end
        figure; imagesc(fbsyMat,[0 4]); colorbar; daspect([1 1 1]); title([roitype ' fb' num2str(x) '-' num2str(y) ' (in/out) long distance neurons synapse count']);
        ylabel('output neuron'); xlabel('input neuron');

        % output post-synapse cloud of neuron fibers for confirmation
        niifname = ['results/fb/' idstr '_fb' num2str(x) '-' num2str(y) '_postsynFDACal.nii'];
        if ~exist([niifname '.gz'],'file')
            [hbV, hbinfo] = calcFiberSynapseV(fbsyC, SlocFc);
            niftiwrite(hbV,niifname,hbinfo,'Compressed',true);
        end
    end
end

function [fbsyMat, fbsyC] = calcFiberSynapseMat(nids, StoN, Sdir, StoS, ssrate, sstraced)
    Sid = uint32(1:length(StoN))';
    nidNum = length(nids);
    fbsyMat = zeros(nidNum,nidNum);
    fbsyC = cell(nidNum,nidNum);
    parfor i=1:nidNum
        disp(['fiber synapse process i=' num2str(i)]);
        logi = ismember(StoN,nids(i)); % find synapses which belong to a target neuron(i)
        presids = Sid(logi & Sdir==1); % output of a neuron(i)
        sslogi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse in any ROI
        cpostsids = StoS(sslogi & ssrate & sstraced,2); % post-sid is unique
        if isempty(cpostsids), continue; end
        for j=1:nidNum
            if i==j, continue; end
            logi = ismember(StoN,nids(j)); % find synapses which belong to a target neuron(j)
            postsids = Sid(logi & Sdir==2); % input of a neuron(j)
            plogi = ismember(postsids,cpostsids); % get common between postsids & cpostsids in any ROI
            fbsyC{i,j} = postsids(plogi);
            fbsyMat(i,j) = length(fbsyC{i,j});
        end
    end
end

function [innids, outnids, infbsyMat, outfbsyMat, Cinfbsy, Coutfbsy] = findFiberCommonInOut(nids, Nid, StoN, Sdir, StoS, ssrate, sstraced)
    Sid = uint32(1:length(StoN))';
    nidNum = length(nids);
    inlogi = ones(length(Nid),1,'logical');
    outlogi = ones(length(Nid),1,'logical');
    Cpresids = cell(1,nidNum);
    Cpostsids = cell(1,nidNum);
    for i=1:nidNum
        disp(['fiber common in/out process i=' num2str(i)]);
        logi = ismember(StoN,nids(i)); % find synapses which belong to a target neuron(i)
        presids = Sid(logi & Sdir==1); % output of a neuron(i)
        sslogi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
        cpostsids = StoS(sslogi & ssrate & sstraced,2); % post-sid is unique
        Cpostsids{i} = cpostsids;
        outlogi = outlogi & ismember(Nid,StoN(cpostsids));

        postsids = Sid(logi & Sdir==2); % input of a neuron(i)
        sslogi = ismember(StoS(:,2),postsids); % get post-synapse to connected pre-synapse
        cpresids = StoS(sslogi & ssrate & sstraced,1); % post-sid is unique
        Cpresids{i} = cpresids;
        inlogi = inlogi & ismember(Nid,StoN(cpresids));
    end
    innids = Nid(inlogi);
    outnids = Nid(outlogi);
    infbsyMat = zeros(length(innids),nidNum);
    Cinfbsy = cell(length(innids),nidNum);
    for i=1:length(innids)
        logi = ismember(StoN,innids(i)); % find synapses which belong to a target neuron(i)
        presids = Sid(logi & Sdir==1); 
        if isempty(presids), continue; end
        for j=1:nidNum
            slogi = ismember(presids, Cpresids{j});
            Cinfbsy{i,j} = presids(slogi);
            infbsyMat(i,j) = length(Cinfbsy{i,j});
        end
    end
    outfbsyMat = zeros(nidNum,length(outnids));
    Coutfbsy = cell(nidNum,length(outnids));
    for j=1:length(outnids)
        logi = ismember(StoN,outnids(j)); % find synapses which belong to a target neuron(i)
        postsids = Sid(logi & Sdir==2); 
        if isempty(postsids), continue; end
        for i=1:nidNum
            slogi = ismember(postsids, Cpostsids{j});
            Coutfbsy{i,j} = postsids(slogi);
            outfbsyMat(i,j) = length(Coutfbsy{i,j});
        end
    end
end

function [hbV, hbinfo] = calcFiberSynapseV(fbsyC, SlocFc)
    fbsy = [];
    nidNum = size(fbsyC,1);
    for i=1:nidNum
        for j=1:nidNum
            if ~isempty(fbsyC{i,j})
                fbsy = [fbsy; fbsyC{i,j}];
            end
        end
    end
    fbsy = unique(fbsy);
    SpostlocFc = SlocFc(fbsy,:);

    hbinfo = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    hbV = niftiread(hbinfo); hbV(:) = 0;
    sz = size(hbV);

    for i=1:size(SpostlocFc,1)
        t = ceil(SpostlocFc(i,:));
        if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
            hbV(t(1),t(2),t(3)) = hbV(t(1),t(2),t(3)) + 1;
        else
            disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
        end
    end
end
