% analyze and plot Neural FC result.
% this script can run after analyzeNeuralFC.m

function plotNeuralFC
    % DBscan param
    epsilon = 3000; % nanometer
    minpts = 1; % set 1, no isolated synapse

    % check NeuralFC of FlyEM
    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    confTh = 60; synTh = 5; % almost flywire codex compatible setting
    conf = getSCconfig('hemi',synTh,scTh);

%    showNeuralFCFw(conf, epsilon, minpts); % no use

    showNeuralDBScanFw(conf, epsilon, minpts); % ext figure.4-1
%    showNeuralDBScanSyCloudFw(conf, epsilon, minpts, [0 0.1]); % ext figure.4-2
%    showNeuralDBScanSyCloudFw(conf, epsilon, minpts, [0.9 1]); % ext figure.4-2

%    showReciprocalDistanceGraphFw(conf); % ext figure.4-1
    showReciprocalDistanceSyCloudFw(conf, 2000); % ext figure.4-2

%    showReciprocalDistanceFw(conf, uint64(425790257), 2000); % APL-R
%    showReciprocalDistanceFw(conf, [], []); % ext figure.4-1
%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                showTriangleDistanceFw(conf,'Feedforward', (p12==1), (p13==1), (p23==1));
            end
        end
    end
%}
%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                showTriangleDistanceFw(conf,'Unicycle', (p12==1), (p13==1), (p23==1))
            end
        end
    end
%}
    % check NeuralFC of FlyWire
    scTh = 140; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
    conf = getSCconfig('wire',synTh,scTh);

%    showNeuralFCFw(conf, epsilon, minpts); % no use

%    showNeuralDBScanFw(conf, epsilon, minpts); % figure.4
%    showNeuralDBScanSyCloudFw(conf, epsilon, minpts, [0 0.1]); % ext figure.4-2
%    showNeuralDBScanSyCloudFw(conf, epsilon, minpts, [0.9 1]); % ext figure.4-2
%    showNeuralDBScanSpidxFw(conf, epsilon, minpts, int64(720575940644632087)); % WAGN figure.5

%    showReciprocalDistanceGraphFw(conf); % figure.4
%    showReciprocalDistanceSyCloudFw(conf, 2000); % ext figure.4-2

%    showReciprocalDistanceFw(conf, int64(720575940628908548), 2000); % CT1-R (2um is good threshold based on histogram) (ext figure.4-2).
%    showReciprocalDistanceFw(conf, int64(720575940613583001), 2000); % APL-R (ext figure.4-2)
 
    showReciprocalDistanceFw(conf, [], []); % show all neuron, top 3 closest
%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                showTriangleDistanceFw(conf,'Feedforward', (p12==1), (p13==1), (p23==1))
            end
        end
    end
%}
%{
    for p12=1:2
        for p13=1:2
            for p23=1:2
                showTriangleDistanceFw(conf,'Unicycle', (p12==1), (p13==1), (p23==1))
            end
        end
    end
%}
end

function showNeuralFC(synTh, confTh, epsilon, minpts)
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    tNid = Nid(Nstatus==1); % traced nids
%    writematrix(tNid,'tracedNid.csv'); % tempraly output

    % FlyEM read neural SC
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_neuralInOutVoxels.mat']);
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_neuralInOutDistance.mat']);
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

    % read neural FC
    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                % output file
                fcname = ['results/neuralfc/' smooth{k} hpfstr nuisance{n} preproc 'hemi' num2str(synTh) 'sr' num2str(confTh) '-fc.mat'];
                load(fcname);

                for i=1:200
                    if ~isempty(mFz{i})
                        type = ''; % neuro transmitter type is unknown
                        F = mFz{i}(:); idx = find(~isinf(F)); F = F(idx);
                        Dt = D{i}(:); Dt = sqrt(Dt(idx));
                        r = corr(F,Dt);
                        disp([num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type ' in=' num2str(sum(inCount{i})) ' (' num2str(size(mFz{i},1)) ...
                            ')  out=' num2str(sum(outCount{i})) ' ('  num2str(size(mFz{i},2)) ')']);

                        % scatter correlation between FC vs. Distance
%                        figure; scatter(F,Dt);
%                        figure; histogram2(F,Dt,20,'Normalization','probability');
%                        xlim([0 4.5]); ylim([1 9]); xlabel('m-FC(z)'); ylabel('sqrt(distance)'); title(['Synaptic FC neuron (' num2str(i) ') r=' num2str(r)]);

                        % calc DBscan synaptic cluster order
                        scinlen = inlen{i};
                        C = DBidx{i};
                        maxcls = max(C);
                        Cidx1 = C(1:scinlen); Cidx1(Cidx1<-1)=[]; % remove out of brain mask voxel
                        Cidx2 = C(scinlen+1:end); Cidx2(Cidx2<-1)=[]; % remove out of brain mask voxel
                        [s,si1]=sort(Cidx1); si1(s<0)=[]; % remove out of cluster voxels
                        [s,si2]=sort(Cidx2); si2(s<0)=[]; % remove out of cluster voxels

                        % FC matrix and SC matrix
                        figure; imagesc(mFz{i}(si1,si2)); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
%                        figure; imagesc(D{i}(si1,si2)); colorbar; title(['Sy Dist ' num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type]);

                        % calc DBscan synaptic cluster
                        CF = calcDBscanSycluster(mFz{i}, maxcls, Cidx1, Cidx2);
                        figure; imagesc(CF,[0 4]); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
                        if maxcls>1, xticks([1:maxcls]); yticks([1:maxcls]); end

                        % plot 3D points
                        plot3Dpoints(sz, inIdx{i}, outIdx{i}, C, synTh, confTh, i, tNid(i), type, false);
                    end
                end
            end
        end
    end
end

function showNeuralFCFw(conf, epsilon, minpts)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read neural SC
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat']);
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutDistance.mat']);
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

    % read neural FC
    for h=1:length(hpfTh)
        hpfstr = '';
        if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        for k=1:length(smooth)
            for n=1:length(nuisance)
                % output file
                fcname = ['results/neuralfc/' smooth{k} hpfstr nuisance{n} preproc 'wire' num2str(synTh) 'sr' num2str(scoreTh) '-fc.mat'];
                load(fcname);

                for i=1:200
                    if ~isempty(mFz{i})
                        t = Ntype(i)+1;
                        type = tlabels{t};
                        F = mFz{i}(:); idx = find(~isinf(F)); F = F(idx);
                        Dt = D{i}(:); Dt = sqrt(Dt(idx));
                        r = corr(F,Dt);
                        disp([num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type ' in=' num2str(sum(inCount{i})) ' (' num2str(size(mFz{i},1)) ...
                            ')  out=' num2str(sum(outCount{i})) ' ('  num2str(size(mFz{i},2)) ')']);

                        % scatter correlation between FC vs. Distance
%                        figure; scatter(F,Dt);
%                        figure; histogram2(F,Dt,20,'Normalization','probability');
%                        xlim([0 4.5]); ylim([1 9]); xlabel('m-FC(z)'); ylabel('sqrt(distance)'); title(['Synaptic FC neuron (' num2str(i) ') r=' num2str(r)]);

                        % calc DBscan synaptic cluster order
                        scinlen = inlen{i};
                        C = DBidx{i};
                        maxcls = max(C);
                        Cidx1 = C(1:scinlen); Cidx1(Cidx1<-1)=[]; % remove out of brain mask voxel
                        Cidx2 = C(scinlen+1:end); Cidx2(Cidx2<-1)=[]; % remove out of brain mask voxel
                        [s,si1]=sort(Cidx1); si1(s<0)=[]; % remove out of cluster voxels
                        [s,si2]=sort(Cidx2); si2(s<0)=[]; % remove out of cluster voxels

                        % FC matrix and SC matrix
                        figure; imagesc(mFz{i}(si1,si2)); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
%                        figure; imagesc(D{i}(si1,si2)); colorbar; title(['Sy Dist ' num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type]);

                        % calc DBscan synaptic cluster
                        CF = calcDBscanSycluster(mFz{i}, maxcls, Cidx1, Cidx2);
                        figure; imagesc(CF,[0 4]); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
                        if maxcls>1, xticks([1:maxcls]); yticks([1:maxcls]); end

                        % plot 3D points
                        plot3Dpoints(sz, inIdx{i}, outIdx{i}, C, synTh, scoreTh, i, Nid(i), type, true);
                    end
                end
            end
        end
    end
end

function showNeuralDBScanFw(conf, epsilon, minpts)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.
    validSidx = Sidx(valid & score);
    syVnum = 2*length(validSidx);
    disp([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' ep' num2str(epsilon) 'nm cl' num2str(minpts) ' : valid pre & post synapse num=' num2str(syVnum) ', neuron num=' num2str(length(Nid))])

    % FlyEM read synapse info
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % FlyWire read neural SC
    DBcount = {};
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

    % pre-post-synapse separate index
    load([conf.neuSepidxFile num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);
    Nspidx = double(Nspidx) / 10000;
    Nspidx(Nspidx<0) = nan;

    % whole brain mesh
    mesh = load(conf.brainMeshFile);

    nlen = length(Nid);
    clcount(clcount<0) = 0; % remove -1 case
    clAll = sum(clcount,1);

    % synapse count & calc pre-post-synapse separation index
    syCount = zeros(nlen,6,'single');
    Ldbs = nan(nlen,1,'single');
    H = nan(nlen,1,'single');
    for i=1:nlen
        if isempty(DBcount{i}), continue; end
        cmat = double(DBcount{i});
        prelogi = (cmat(:,1)==0);  % pre-synapse is zero (pure input cluster)
        postlogi = (cmat(:,2)==0); % post-synapse is zero (pure output cluster)
        syCount(i,1) = sum(cmat,'all');
        syCount(i,2:3) = sum(cmat,1);
        syCount(i,4) = sum(cmat(prelogi,:),'all');
        syCount(i,5) = sum(cmat(postlogi,:),'all');
        syCount(i,6) = size(cmat,1);     % cluster size

        % pre-post-synapse separation index based on DBscan clustering
        Ni = sum(cmat,2);
        Ldbs(i) = sum(abs(cmat(:,1)-cmat(:,2))./Ni .* (Ni./syCount(i,1))) ;   % linear version (cluster synaptic weight)

        % segregation index (Schneider-mizell et al., 2016)
        Pi = cmat(:,2) ./ Ni;
        Si = -(Pi.*log(Pi) + (1-Pi).*log(1-Pi)); % this could be NaN
        Si(isnan(Si)) = 0; % segregated 
        S = nansum(Ni .* Si) / syCount(i,1);
        Fi = sum(Pi.*Ni) / syCount(i,1);
%        Fi = syCount(i,2) / syCount(i,1);
        Snorm = -(Fi*log(Fi) + (1-Fi)*log(1-Fi));
        if isnan(Snorm)
            H(i) = 1; % segregated
        else
            H(i) = 1 - S/Snorm;
        end
    end

    % show histogram
    edges = 0:0.05:1;
    N = [];
    for i=1:length(tlabels)
        idx = find(Ntype==(i-1));
        h = histcounts(Nspidx(idx),edges);
        N = [N, h'];
    end
    % Sdbs is better than linear version with FlyWire and hemibrain
    figure; bar(N,'stacked','LineWidth',0.1); legend(tlabels); xticklabels(''); xlabel('separation index [0 1]'); ylabel('neuron count');
    title([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : pre-post-synapse separation index histogram']);
%    figure; histogram(Ldbs);

    % show histogram (segregation index)
    edges = 0:0.05:1;
    N = [];
    for i=1:length(tlabels)
        idx = find(Ntype==(i-1));
        h = histcounts(H(idx),edges);
        N = [N, h'];
    end
    % Sdbs is better than linear version with FlyWire and hemibrain
    figure; bar(N,'stacked','LineWidth',0.1); legend(tlabels); xticklabels(''); xlabel('segregation index [0 1]'); ylabel('neuron count');
    title([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : pre-post-synapse segregation index histogram']);
%    figure; histogram(Ldbs);

    Nspidx(isnan(Nspidx)) = -1;
    NsywMixScore(isnan(NsywMixScore)) = -1;
    NsywSepScore(isnan(NsywSepScore)) = -1;

    syAll = sum(syCount,1);
    disp(['total cluster (synapse)=' num2str(clAll(1)) ' (' num2str(syAll(1)) '), totalIn=' num2str(clAll(2)) ' (' num2str(syAll(4)) '), totalOut=' num2str(clAll(3)) ' (' num2str(syAll(5)) ')']);
    disp(['in-cluster rate=' num2str(clAll(2)/clAll(1)) ', out-cluster rate=' num2str(clAll(3)/clAll(1))]);
    disp(['cluster belong synapse rate=' num2str(syAll(1)/syVnum) ', in-cluster synapse rate=' num2str(syAll(4)/syVnum) ', out-cluster synapse rate=' num2str(syAll(5)/syVnum)]);

    nnum = length(find(Nspidx>=0));
    spnum = length(find(Nspidx==1));
    disp(['full separation neuron rate=' num2str(spnum/nnum) ' (' num2str(spnum) '/' num2str(nnum) ')']);

    % show pure separation neurons
%{
    s1idx = find(Nspidx==1);
    [syCnt1,dIdx] = sort(syCount(s1idx,1),'descend');
    for i=1:10
        k = s1idx(dIdx(i));
        nid = Nid(k);
%        if exist('Ncrop','var') && Ncrop(k)==1, continue; end % ignore cropped body.

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(nid) '.swc'], true);

        % get all connected synapses
        prelogi = ismember(preNidx,k);
        postlogi = ismember(postNidx,k);
        loc1 = Spreloc(prelogi & valid & score,:);
        loc2 = Spostloc(postlogi & valid & score,:);

        tstr = [num2str(i) ') k=' num2str(k) ' nid=' num2str(nid) ' (' tlabels{Ntype(k)+1} ') sycount=' num2str(syCount(k,1)) ' (' num2str(syCount(k,2)) '/' num2str(syCount(k,3)) ';' num2str(syCount(k,6)) ')' ...
            ', nSPidx=' num2str(Nspidx(k)) ', syMIXscore=' num2str(NsywMixScore(k)) ', Ldbs=' num2str(Ldbs(k))];
        disp(tstr);
        figure; plotBrainSwc(swc, mesh, conf.brainMeshView, loc1, loc2, tstr);
    end
%}
    % show most mixed and large synapse neurons
%%{
    [dsyMdbs,dIdx] = sort(NsywMixScore,'descend');
    for i=1:50
        k = dIdx(i);
        nid = Nid(k);
        if exist('Ncrop','var') && Ncrop(k)==1, continue; end % ignore cropped body.

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(nid) '.swc'], true);

        % get all connected synapses
        prelogi = ismember(preNidx,k);
        postlogi = ismember(postNidx,k);
        loc1 = Spreloc(prelogi & valid & score,:);
        loc2 = Spostloc(postlogi & valid & score,:);

        tstr = [num2str(i) ') k=' num2str(k) ' nid=' num2str(nid) ' (' tlabels{Ntype(k)+1} ') sycount=' num2str(syCount(k,1)) ' (' num2str(syCount(k,2)) '/' num2str(syCount(k,3)) ';' num2str(syCount(k,6)) ')' ...
            ', nSPidx=' num2str(Nspidx(k)) ', syMIXscore=' num2str(NsywMixScore(k)) ', Ldbs=' num2str(Ldbs(k))];
        disp(tstr);
        figure; plotBrainSwc(swc, mesh, conf.brainMeshView, loc1, loc2, tstr, [conf.scname 'syMix' num2str(i) '.png']);
    end
%}
    % show most separated and large synapse neurons
%%{
    [dsySdbs,dIdx] = sort(NsywSepScore,'descend');
    for i=1:50
        k = dIdx(i);
        nid = Nid(k);
        if exist('Ncrop','var') && Ncrop(k)==1, continue; end % ignore cropped body.

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(nid) '.swc'], true);

        % get all connected synapses
        prelogi = ismember(preNidx,k);
        postlogi = ismember(postNidx,k);
        loc1 = Spreloc(prelogi & valid & score,:);
        loc2 = Spostloc(postlogi & valid & score,:);

        tstr = [num2str(i) ') k=' num2str(k) ' nid=' num2str(nid) ' (' tlabels{Ntype(k)+1} ') sycount=' num2str(syCount(k,1)) ' (' num2str(syCount(k,2)) '/' num2str(syCount(k,3)) ';' num2str(syCount(k,6)) ')' ...
            ', nSPidx=' num2str(Nspidx(k)) ', sySPscore=' num2str(NsywSepScore(k)) ', Ldbs=' num2str(Ldbs(k))];
        disp(tstr);
        figure; plotBrainSwc(swc, mesh, conf.brainMeshView, loc1, loc2, tstr, [conf.scname 'sySP' num2str(i) '.png']);
    end
%}
end

function plotBrainSwc(swc, mesh, vp, loc1, loc2, tstr, savefile)
    if nargin < 7, savefile = []; end
    plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid off; axis image; alpha(.1);
    if ~isempty(vp), view(vp(1),vp(2)); end
    xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]'); title(tstr);
    hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),8,'red','filled'); hold off;  % pre. output
    hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),8,'blue','filled'); hold off; % post. input
    alpha(.5);
    if ~isempty(mesh)
        hold on; patch('Faces',mesh.F,'Vertices',mesh.V,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'EdgeColor','none','LineStyle','none'); hold off;
    end
    if ~isempty(savefile)
        h = gca;
        h.XAxis.Visible = 'off';
        h.YAxis.Visible = 'off';
        h.ZAxis.Visible = 'off';
        title([]);
        switch(savefile(1:4))
        case 'hemi'
            set(gcf,'Position',[10 10 800 800]);
        case 'wire'
            set(gcf,'Position',[10 10 1500 700]);
        end
        saveas(gca,savefile)
    end
end

function showNeuralDBScanSpidxFw(conf, epsilon, minpts, nids)
    if nargin < 4, nids = []; end

    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.
    validSidx = Sidx(valid & score);
    syVnum = 2*length(validSidx);
    disp([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' ep' num2str(epsilon) 'nm cl' num2str(minpts) ' : valid pre & post synapse num=' num2str(syVnum) ', neuron num=' num2str(length(Nid))])

    % FlyEM read synapse info
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % FlyWire read neural SC
    DBcount = {};
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

    % pre-post-synapse separate index
    load([conf.neuSepidxFile num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);
    Nspidx = double(Nspidx) / 10000;
    Nspidx(Nspidx<0) = nan;
    load([conf.sySepidxFile num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);
    preSpidx = single(preSpidx) / 10000;
    preSpidx(preSpidx<0) = nan;
    postSpidx = single(postSpidx) / 10000;
    postSpidx(postSpidx<0) = nan;

    % whole brain mesh
    mesh = load(conf.brainMeshFile);

    if isempty(nids), nids = Nid; end
    nlen = length(nids);
    for i=1:nlen
        nid = nids(i);
        k = find(Nid==nid);
        if exist('Ncrop','var') && Ncrop(k)==1, continue; end % ignore cropped body.

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(nid) '.swc'], true);

        % get all connected synapses
        prelogi = ismember(preNidx,k);
        postlogi = ismember(postNidx,k);
        loc1 = Spreloc(prelogi & valid & score,:);
        loc2 = Spostloc(postlogi & valid & score,:);
        spidx1 = preSpidx(prelogi & valid & score);
        spidx2 = postSpidx(postlogi & valid & score);

        tstr = [num2str(i) ') k=' num2str(k) ' nid=' num2str(nid) ' (' tlabels{Ntype(k)+1} ') spidx=' num2str(Nspidx(k))];
        disp(tstr);

        figure; plotSwc(swc, [0.7 0.7 1], 0.2, true); view(3); grid off; axis image; alpha(.1);
        xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]'); title(tstr);
        hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),12,spidx1,'filled'); hold off; % pre. output
        hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),16,spidx2,'+'); hold off; % post. input
        alpha(.5); colormap(gca,'jet');
    end
end

function showNeuralDBScanSyCloudFw(conf, epsilon, minpts, range)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/nifti/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '_spidx' num2str(range(1)) '-' num2str(range(2)) 'Synapses.nii'];
    if exist([fname '.gz'],'file'), return; end

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % pre-post-synapse separate index
    load([conf.sySepidxFile num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);
    postSpidx = double(postSpidx) / 10000;
    postSpidx(postSpidx<0) = nan;
    preSpidx = double(preSpidx) / 10000;
    preSpidx(preSpidx<0) = nan;

    % read synapse location in FDA
    load(conf.syprelocFdaFile);
    load(conf.sypostlocFdaFile);

    logi = (range(1) <= preSpidx) & (preSpidx <= range(2));
    presidx = Sidx(logi & valid & score);
    logi = (range(1) <= postSpidx) & (postSpidx <= range(2));
    postsidx = Sidx(logi & valid & score);

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0;
    sz = size(Vt);

    % output synaptic cloud of spindex within range.
    conSlocFc = SprelocFc(presidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
    for j=1:size(conSlocFc,1)
        t = ceil(conSlocFc(j,:));
        if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
            V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
        else
            disp(['out of bounds ) ' num2str(t)]);
        end
    end
    conSlocFc = SpostlocFc(postsidx,:);
    for j=1:size(conSlocFc,1)
        t = ceil(conSlocFc(j,:));
        if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
            V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
        else
            disp(['out of bounds ) ' num2str(t)]);
        end
    end
    niftiwrite(V,fname,info,'Compressed',true);
end

function showReciprocalDistanceSyCloudFw(conf, dist)
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    fname = ['results/nifti/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_rcdist' num2str(dist) 'Synapses.nii'];
    if exist([fname '.gz'],'file'), return; end

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % pre-post-synapse reciprocal distance index
    load([conf.syReciFile num2str(synTh) 'sr' num2str(scoreTh) '.mat']);

    % read synapse location in FDA
    load(conf.syprelocFdaFile);
    load(conf.sypostlocFdaFile);

    logi = (SrcpostCloseDist <= dist);
    presidx = Sidx(logi & valid & score);
    logi = (SrcpreCloseDist <= dist);
    postsidx = Sidx(logi & valid & score);

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0;
    sz = size(Vt);

    % output synaptic cloud of spindex within range.
    conSlocFc = SprelocFc(presidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
    for j=1:size(conSlocFc,1)
        t = ceil(conSlocFc(j,:));
        if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
            V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
        else
            disp(['out of bounds ) ' num2str(t)]);
        end
    end
    conSlocFc = SpostlocFc(postsidx,:);
    for j=1:size(conSlocFc,1)
        t = ceil(conSlocFc(j,:));
        if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
            V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
        else
            disp(['out of bounds ) ' num2str(t)]);
        end
    end
    niftiwrite(V,fname,info,'Compressed',true);
end

function showReciprocalDistance(synTh, confTh)
    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 
    tracedNids = Nid(Nstatus==1);
    nlen = length(tracedNids);

    % FlyEM read synapse info
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

    % FlyWire read neural SC
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_reciprocalConnections.mat']);
%    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat']);
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_reciprocalDistances.mat']);

    for i=1:nlen
        if isempty(rcpreSids{i}), continue; end
        disp(['plot closest reciprocal synapses (' num2str(i) ') nid=' num2str(tracedNids(i))]);

        % close from reciprocal pre-synapse on tracedNids(i)
        prelocs = Sloc(rcpreSids{i},:);        % pre-synapse on tracedNids(i) used by reciprocal 
        rcpreDs = rcpreCloseDist{i};
        postclocs = Sloc(rcpreCloseSids{i},:); % post-synapse on tracedNids(i) used by reciprocal 
        % close from reciprocal post-synapse on tracedNids(i)
        postlocs = Sloc(rcpostSids{i},:);      % post-synapse on tracedNids(i) used by reciprocal 
        rcpostDs = rcpostCloseDist{i};
        preclocs = Sloc(rcpostCloseSids{i},:); % pre-synapse on tracedNids(i) used by reciprocal 

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc(['swc/hemibrain_v1_2/' num2str(tracedNids(i)) '.swc'], true);
%%{
        [md, didx] = sort(rcpreDs);
        for k=1:length(didx)
            if k>3, break; end
            j = didx(k);
            postcsid = rcpreCloseSids{i}(j);
            sslogi = ismember(StoS(:,2),postcsid);
            rcnid = StoN(StoS(sslogi,1));
            disp([num2str(i) '-' num2str(k) ' nids=' num2str(tracedNids(i)) ',' num2str(rcnid) ' dist=' num2str(rcpreDs(j)) 'nm']);

            % get all connected synapses
            prelogi2 = (StoN==rcnid) & (Sdir==1);
            postlogi2 = (StoN==rcnid) & (Sdir==2);
            sslogi = ismember(StoS(:,2),Sid(postlogi2));
            presids2 = StoS(sslogi,1);
            sslogi = ismember(StoS(:,1),Sid(prelogi2));
            postsids2 = StoS(sslogi,2);
            logis = ismember(rcpreSids{i},presids2);
            loc1 = Sloc(rcpreSids{i}(logis),:);
            logis = ismember(rcpostSids{i},postsids2);
            loc2 = Sloc(rcpostSids{i}(logis),:);

            % plot neuron and synapses
            swcname = ['swc/hemibrain_v1_2/' num2str(rcnid) '.swc'];
            if ~exist(swcname,'file'), continue; end
            swc2 = loadSwc(swcname, true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' nids=' num2str(tracedNids(i)) ',' num2str(rcnid) ' dist=' num2str(rcpreDs(j)) 'nm']);
            xlabel('x [8nm]'); ylabel('y [8nm]'); zlabel('z [8nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(prelocs(j,1),prelocs(j,2),prelocs(j,3),24,'red','filled'); hold off;
            hold on; scatter3(postclocs(j,1),postclocs(j,2),postclocs(j,3),24,'blue','filled'); hold off;
        end
%}
%%{
        [md, didx] = sort(rcpostDs);
        for k=1:length(didx)
            if k>3, break; end
            j = didx(k);
            postsid = rcpostSids{i}(j); % need post-synapse to determine one reciprocal pre-synapse
            sslogi = ismember(StoS(:,2),postsid);
            rcnid = StoN(StoS(sslogi,1));
            disp([num2str(i) '-' num2str(k) ' nids=' num2str(tracedNids(i)) ',' num2str(rcnid) ' dist=' num2str(rcpostDs(j)) 'nm']);

            % get all connected synapses
            prelogi2 = (StoN==rcnid) & (Sdir==1);
            postlogi2 = (StoN==rcnid) & (Sdir==2);
            sslogi = ismember(StoS(:,2),Sid(postlogi2));
            presids2 = StoS(sslogi,1);
            sslogi = ismember(StoS(:,1),Sid(prelogi2));
            postsids2 = StoS(sslogi,2);
            logis = ismember(rcpreSids{i},presids2);
            loc1 = Sloc(rcpreSids{i}(logis),:);
            logis = ismember(rcpostSids{i},postsids2);
            loc2 = Sloc(rcpostSids{i}(logis),:);

            % plot neuron and synapses
            swcname = ['swc/hemibrain_v1_2/' num2str(rcnid) '.swc'];
            if ~exist(swcname,'file'), continue; end
            swc2 = loadSwc(swcname, true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' nids=' num2str(tracedNids(i)) ',' num2str(rcnid) ' dist=' num2str(rcpostDs(j)) 'nm']);
            xlabel('x [8nm]'); ylabel('y [8nm]'); zlabel('z [8nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(preclocs(j,1),preclocs(j,2),preclocs(j,3),24,'red','filled'); hold off;
            hold on; scatter3(postlocs(j,1),postlocs(j,2),postlocs(j,3),24,'blue','filled'); hold off;
        end
%}
        % plot all neurons and synapses related reciprocal (faster line plot). This one is heavy.
%%{
        figure;
        for k=1:length(rcNids{i})
            swcname = ['swc/hemibrain_v1_2/' num2str(rcNids{i}(k)) '.swc'];
            if ~exist(swcname,'file'), continue; end
            swc2 = loadSwc(swcname);
            hold on; plotSwc(swc2, [0.8 0.8 0.8], 0.1); alpha(.1); hold off;
        end
        hold on; plotSwc(swc, [0.5 0.5 1], 1.5); hold off; view(3); grid on; axis image;
        title(['reciprocal connections(' num2str(i) ') nids=' num2str(tracedNids(i)) ',' num2str(rcnid)]);
        xlabel('x [8nm]'); ylabel('y [8nm]'); zlabel('z [8nm]');
        hold on; scatter3(prelocs(:,1),prelocs(:,2),prelocs(:,3),18,'red'); hold off;
        hold on; scatter3(postlocs(:,1),postlocs(:,2),postlocs(:,3),18,'blue'); hold off;
%}
    end
end

function showReciprocalDistanceFw(conf, targetNid, dth)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % FlyEM read synapse info
    Spostloc = []; Spreloc = [];
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % FlyWire read neural SC
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat']);
%    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat']);
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalDistances.mat']);

    if ~isempty(targetNid)
        idx = find(Nid==targetNid);
    else
        nlen = length(Nid);
        idx = 1:10:nlen;
    end
    for i=idx
        if isempty(rcpreSidx{i}), continue; end
        disp(['plot closest reciprocal synapses (' num2str(i) ') nid=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ')']);

        precd = rcpreCloseDist{i};
        postcd = rcpostCloseDist{i};
        if length(precd) > length(postcd)
            % close from reciprocal pre-synapse on Nid(i)
            preclocs = Spreloc(rcpreSidx{i},:);         % reciprocal pre-synapse on Nid(i)
            postclocs = Spostloc(rcpreCloseSidx{i},:);     % reci post sidx on Nid(i)
            postcsidxs = rcpreCloseSidx{i};
            rcprenum = length(unique(rcpreSidx{i})); rcpostnum = length(unique(rcpreCloseSidx{i}));
            rcDs = precd;
        else
            % close from reciprocal post-synapse on Nid(i)
            preclocs = Spreloc(rcpostCloseSidx{i},:);       % reci pre sidx on Nid(i)
            postclocs = Spostloc(rcpostSidx{i},:);      % reciprocal post-synapse on Nid(i)
            postcsidxs = rcpostSidx{i};
            rcprenum = length(unique(rcpostCloseSidx{i})); rcpostnum = length(unique(rcpostSidx{i}));
            rcDs = postcd;
        end

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(Nid(i)) '.swc'], true);
        prelogi = ismember(preNidx,i);
        postlogi = ismember(postNidx,i);
        sprenum = length(Sidx(prelogi & valid & score,:));
        spostnum = length(Sidx(postlogi & valid & score,:));

        % show neuron and thresholded distance synapses
        if ~isempty(dth)
            didx = find(rcDs<dth); % thresholded (reciprocal-synapse)
            if length(precd) > length(postcd)
                rcthprenum = length(unique(rcpreSidx{i}(didx))); rcthpostnum = length(unique(rcDs(didx)));
            else
                rcthprenum = length(unique(rcDs(didx))); rcthpostnum = length(unique(rcpostSidx{i}(didx)));
            end
            
            disp([num2str(i) ' nids=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ') close reci synapses=' num2str(rcthprenum+rcthpostnum) '/' num2str(rcprenum+rcpostnum) '/' num2str(sprenum+spostnum)]);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image; alpha(.1);
            xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]'); title([num2str(i) ' nids=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ')']);
            hold on; scatter3(preclocs(didx,1),preclocs(didx,2),preclocs(didx,3),8,'red','filled'); hold off;     % pre. output
            hold on; scatter3(postclocs(didx,1),postclocs(didx,2),postclocs(didx,3),8,'blue','filled'); hold off; % post. input
            alpha(.5);

            figure; histogram(rcDs); xlim([0 10000]);
            title([num2str(i) ' nids=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ') reci distance histogram'])
        end

        % show neuron and thresholded distance synapses
        [md, didx] = sort(rcDs);
        for k=1:length(didx)
            if k>3, break; end
            j = didx(k);
            rcnidx = preNidx(postcsidxs(j));
            prelogi2 = ismember(preNidx,rcnidx);
            postlogi2 = ismember(postNidx,rcnidx);

            % get all connected synapses
            loc1 = Spreloc(prelogi & postlogi2 & valid & score,:);
            loc2 = Spostloc(postlogi & prelogi2 & valid & score,:);

            disp([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') ' ...
                'con=' num2str(size(loc1,1)+size(loc2,1)) ' dist=' num2str(rcDs(j)) 'nm']);
            disp(['    pre: ' num2str(double(preclocs(j,:))./conf.swcSize) ',  post: ' num2str(double(postclocs(j,:))./conf.swcSize)]);
            if abs(preclocs(j,3)-postclocs(j,3))<=1 && rcDs(j) < 300
                a=0; % for break
            end
%{
            % plot neuron and synapses
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(rcnidx)) '.swc'], true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') ' ...
                'con=' num2str(size(loc1,1)+size(loc2,1)) ' dist=' num2str(rcDs(j)) 'nm']);
            xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(preclocs(j,1),preclocs(j,2),preclocs(j,3),24,'red','filled'); hold off;     % pre. output
            hold on; scatter3(postclocs(j,1),postclocs(j,2),postclocs(j,3),24,'blue','filled'); hold off; % post. input
%}
        end
%{
        % plot all neurons and synapses related reciprocal (faster line plot). This one is heavy.
        figure;
        for k=1:length(rcNidx{i})
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(rcNidx{i}(k))) '.swc']);
            hold on; plotSwc(swc2, [0.8 0.8 0.8], 0.1); alpha(.1); hold off;
        end
        hold on; plotSwc(swc, [0.5 0.5 1], 1.5); hold off; view(3); grid on; axis image;
        title(['reciprocal connections(' num2str(i) ') nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx))]);
        xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
        hold on; scatter3(preclocs(:,1),preclocs(:,2),preclocs(:,3),18,'red'); hold off;
        hold on; scatter3(postclocs(:,1),postclocs(:,2),postclocs(:,3),18,'blue'); hold off;
%}
    end
end

function showReciprocalDistanceGraphFw(conf)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % FlyWire read neural SC
    rcpreSidx = {};
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat']);
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalDistances.mat']);

    allDs = []; allTypes = [];
    for i=1:length(Nid)
        if isempty(rcpreSidx{i}), continue; end

        presidxs = rcpreSidx{i};         % reciprocal pre-synapse on Nid(i)
        postsidxs = rcpostSidx{i};
        rcDs1 = rcpreCloseDist{i};
        rcDs2 = rcpostCloseDist{i};
        allDs = [allDs; rcDs1; rcDs2];

        allTypes = [allTypes; ones(length(rcDs1)+length(rcDs2),1,'uint16') * Ntype(i)];
        rcprenum = length(unique(presidxs)); rcpostnum = length(unique(postsidxs));
    end

    % show histogram (full range)
    edges = 0:10000:600000;
    N = [];
    for i=1:length(tlabels)
        idx = find(allTypes==(i-1));
        h = histcounts(allDs(idx),edges);
        N = [N, h'];
    end
    figure; bar(N,'stacked','LineWidth',0.1); legend(tlabels); xlabel('distance [*10um]'); ylabel('synapse count');
    title([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : reciprocal synapse distance histogram']);

    % show histogram (less than 20 um)
    edges = 0:1000:20000;
    N = [];
    for i=1:length(tlabels)
        idx = find(allTypes==(i-1));
        h = histcounts(allDs(idx),edges);
        N = [N, h'];
    end
    % Sdbs is better than linear version with FlyWire and hemibrain
    figure; bar(N,'stacked','LineWidth',0.1); legend(tlabels); xlabel('distance [*1um]'); ylabel('synapse count');
    title([conf.scname num2str(synTh) 'sr' num2str(scoreTh) ' : reciprocal synapse distance histogram']);
end

function showTriangleDistanceFw(conf, trType, isPure12, isPure13, isPure23)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    synTh = conf.synTh;
    scoreTh = conf.scoreTh;
    scname = conf.scname;
    if isPure12, pstr='Pu'; else pstr='Rc'; end
    if isPure13, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end
    if isPure23, pstr=[pstr 'Pu']; else pstr=[pstr 'Rc']; end

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % FlyEM read synapse info
    Spostloc = []; Spreloc = [];
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % FlyWire read neural SC
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangle' trType pstr '.mat']);
%    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat']);
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_triangle' trType pstr 'Dists.mat']);
    load(['results/neuralsc/' scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalConnections.mat']);

    nlen = length(Nid);
    for i=1:1:nlen
        triCdist = triCloseDist{i};
        triCsidx = triCloseSidx{i};
        if isempty(triCdist), continue; end
        disp(['plot closest triangle synapses (' num2str(i) ') nid=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ')']);

        pD = []; pSidx = [];
        for j=1:size(triCdist,1)
            if ~isempty(triCdist{j})
                pD = [pD; triCdist{j}];
                pSidx = [pSidx; triCsidx{j}];
            end
        end

        % load swc
        swc = loadSwc([conf.swcPath '/' num2str(Nid(i)) '.swc'], true);

        prelogi1 = ismember(preNidx,i);
        postlogi1 = ismember(postNidx,i);

        [pDs, didx] = sort(pD);
        for k=1:length(didx)
            if k>5, break; end
            j = didx(k);
            sidx = pSidx(j,:);
            prnidx = preNidx(sidx);
            nidx = postNidx(sidx);
            if strcmp(trType,'Unicycle')
                tmp = nidx(2); nidx(2) = prnidx(2); prnidx(2) = tmp;
            end
            if prnidx(1) ~= i || prnidx(2) ~= i || prnidx(3) ~= nidx(1) || nidx(2) ~= nidx(3)
                disp(['bad i=' num2str(i) ' and nidx1=' num2str(nidx(1))]);
            end

            prelocs = Spreloc(sidx,:);
            postlocs = Spostloc(sidx,:);

            prelogi2 = ismember(preNidx,nidx(1));
            postlogi2 = ismember(postNidx,nidx(1));
            prelogi3 = ismember(preNidx,nidx(3));
            postlogi3 = ismember(postNidx,nidx(3));

            % get reciprocal neuron synapses
            loc1 = zeros(0,3,'int32'); loc2 = zeros(0,3,'int32');
            t1 = Spreloc(prelogi1 & postlogi2 & valid & score,:);
            t2 = Spostloc(postlogi1 & prelogi2 & valid & score,:);
            if length(t1)>=synTh, loc1 = [loc1;t1]; end
            if length(t2)>=synTh, loc2 = [loc2;t2]; end
            switch(trType)
            case 'Feedforward'
                t1 = Spreloc(prelogi1 & postlogi3 & valid & score,:);
                t2 = Spostloc(postlogi1 & prelogi3 & valid & score,:);
            case 'Unicycle'
                t2 = Spreloc(prelogi1 & postlogi3 & valid & score,:);
                t1 = Spostloc(postlogi1 & prelogi3 & valid & score,:);
            end
            if length(t1)>=synTh, loc1 = [loc1;t1]; end
            if length(t2)>=synTh, loc2 = [loc2;t2]; end
            t1 = Spreloc(prelogi2 & postlogi3 & valid & score,:);
            t2 = Spostloc(postlogi2 & prelogi3 & valid & score,:);
            if length(t1)>=synTh, loc1 = [loc1;t1]; end
            if length(t2)>=synTh, loc2 = [loc2;t2]; end

            % check reciprocal
            if isPure12 && isPure13 && isPure23
                if size(loc2,1) > 0
                    disp('!error! reciprocal exist!'); return;
                end
            else
                if size(loc2,1) == 0
                    disp('!error! reciprocal does not exist!'); return;
                end
            end

            disp([num2str(i) '-' num2str(k) ' ' trType pstr ' nids=' num2str(Nid(i)) ',' num2str(Nid(nidx(1))) ',' num2str(Nid(nidx(3))) ...
                ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(nidx(1))+1} ',' tlabels{Ntype(nidx(3))+1} ') ' ...
                'con=' num2str(size(loc1,1)+size(loc2,1)) ' rccon=' num2str(size(loc2,1)) ' dist=' num2str(pD(j)) 'nm']);
            for p=1:3
                disp(['    pre: ' num2str(double(prelocs(p,:))./conf.swcSize) ',  post: ' num2str(double(postlocs(p,:))./conf.swcSize)]);
            end
%{
            % plot neuron and synapses
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(nidx(1))) '.swc'], true);
            swc3 = loadSwc([conf.swcPath '/' num2str(Nid(nidx(3))) '.swc'], true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' ' trType pstr ' nids=' num2str(Nid(i)) ',' num2str(Nid(nidx(1))) ',' num2str(Nid(nidx(3))) ...
                ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(nidx(1))+1} ',' tlabels{Ntype(nidx(3))+1} ') dist=' num2str(pD(j)) 'nm']);
            xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; plotSwc(swc3, [0.8 1 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(prelocs(:,1),prelocs(:,2),prelocs(:,3),24,'red','filled'); hold off;        % pre. output
            hold on; scatter3(postlocs(:,1),postlocs(:,2),postlocs(:,3),24,'blue','filled'); hold off; % post. input
%}
        end
    end
end

function CF = calcDBscanSycluster(F, maxcls, Cidx1, Cidx2)
    F(isinf(F)) = nan;
    CF = nan(maxcls,maxcls,'single');
    for i=1:maxcls
        iIdx = find(Cidx1==i);
        if isempty(iIdx), continue; end
        for j=1:maxcls
            jIdx = find(Cidx2==j);
            if isempty(jIdx), continue; end
            CF(i,j) = nanmean(F(iIdx,jIdx),'all');
        end
    end
end

function plot3Dpoints(sz, inidx, outidx, C, synTh, scoreTh, i, nid, type, isfw)
    idxs = [inidx; outidx];
    X = zeros(length(idxs),3,'single');
    for ii=1:length(idxs)
        [x y z] = ind2sub(sz,idxs(ii));
        X(ii,:) = [x y z];
    end
    if isfw, X(:,3) = sz(3)-X(:,3) + 1; end% flip Z for compatible view of FlyWire codex
    cols = cluster2color(C);
    figure; scatter3(X(:,1),X(:,2),X(:,3),10,cols); daspect([1 1 1]); xlim([1 sz(1)]); ylim([1 sz(2)]); zlim([1 sz(3)]); 
    xlabel('X'); ylabel('Y'); zlabel('Z'); title(['Sy voxles '  num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(nid) ' ' type]);
    if isfw
        view(-180,85); % for FlyWire coxex
    else
        view(0, -3);   % for neuPRINT+
    end
end

function cols = cluster2color(C)
    cols = zeros(length(C),3);
    if any(C<0)
        cols(C<0,:) = repmat([1 1 1],[sum(C<0) 1]); % white is no cluster
    end
    t = getColors(max(C)); % generate the same color in clusters
    for i=1:max(C)
        cols(C==i,:) = repmat(t(i,:),[sum(C==i) 1]);
    end
end

function cols = getColors(maxcol)
    cols = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 .5 0; 0 1 .5; .5 0 1; .5 1 0; 0 .5 1; 1 0 .5; .5 .5 0; 0 .5 .5; .5 0 .5;];
    c2 = cols * .75;
    c3 = cols * .5;
    cols = [cols; c2; c3];
    if size(cols,1) < maxcol
        cols = repmat(cols,[ceil(maxcol/size(cols,1)) 1]);
    end
end

