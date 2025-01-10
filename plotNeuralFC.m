% analyze and plot Neural FC result.
% this script can run after analyzeNeuralFC.m

function plotNeuralFC
    % DBscan param
    epsilon = 5; % micro meter. almost 2 voxels.
    minpts = 3; % set 1, but isolated synapse will be ignored

    % check NeuralFC of FlyEM
    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    confTh = 60; synTh = 5; % almost flywire codex compatible setting
    conf = getSCconfig('hemi',synTh,scTh);
%    checkNeuralFCFw(conf, epsilon, minpts);

    checkReciprocalConnectionFw(conf)

    % check NeuralFC of FlyWire
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
    conf = getSCconfig('wire',synTh,scTh);
%    checkNeuralFCFw(conf, epsilon, minpts);

    checkReciprocalConnectionFw(conf)
end

function checkNeuralFC(synTh, confTh, epsilon, minpts)
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

function checkNeuralFCFw(conf, epsilon, minpts)
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

function checkReciprocalConnection(synTh, confTh)
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
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_neuralReciprocalConnections.mat']);
%    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat']);
    load(['results/neuralsc/hemi' num2str(synTh) 'sr' num2str(confTh) '_reciprocalSynapseDistances.mat']);

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

function checkReciprocalConnectionFw(conf)
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

    % FlyEM read synapse info
    Spostloc = []; Spreloc = [];
    load(conf.sypostlocFile);
    load(conf.syprelocFile);

    % FlyWire read neural SC
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neuralReciprocalConnections.mat']);
%    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_neural_Nin_Nout.mat']);
    load(['results/neuralsc/' conf.scname num2str(synTh) 'sr' num2str(scoreTh) '_reciprocalSynapseDistances.mat']);

    nlen = length(Nid);
    for i=1:nlen
        if isempty(rcpreSidx{i}), continue; end
        disp(['plot closest reciprocal synapses (' num2str(i) ') nid=' num2str(Nid(i)) ' (' tlabels{Ntype(i)+1} ')']);

        % close from reciprocal pre-synapse on Nid(i)
        prelocs = Spreloc(rcpreSidx{i},:);         % reciprocal pre-synapse on Nid(i)
        rcpreDs = rcpreCloseDist{i};
        postclocs = Spostloc(rcpreCloseSidx{i},:); % reci post sidx on Nid(i)
        % close from reciprocal post-synapse on Nid(i)
        postlocs = Spostloc(rcpostSidx{i},:);      % reciprocal post-synapse on Nid(i)
        rcpostDs = rcpostCloseDist{i};
        preclocs = Spreloc(rcpostCloseSidx{i},:);  % reci pre sidx on Nid(i)

        % for detail plot, need to add Trees toolbox (https://www.treestoolbox.org/index.html)
        swc = loadSwc([conf.swcPath '/' num2str(Nid(i)) '.swc'], true);
        prelogi = ismember(preNidx,i);
        postlogi = ismember(postNidx,i);
%%{
        [md, didx] = sort(rcpreDs);
        for k=1:length(didx)
            if k>3, break; end
            j = didx(k);
            postcsix = rcpreCloseSidx{i}(j);
            rcnidx = preNidx(postcsix);
            prelogi2 = ismember(preNidx,rcnidx);
            postlogi2 = ismember(postNidx,rcnidx);
            disp([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') dist=' num2str(rcpreDs(j)) 'nm']);

            % get all connected synapses
            loc1 = Spreloc(prelogi & postlogi2,:);
            loc2 = Spostloc(postlogi & prelogi2,:);

            % plot neuron and synapses
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(rcnidx)) '.swc'], true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') dist=' num2str(rcpreDs(j)) 'nm']);
            xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(prelocs(j,1),prelocs(j,2),prelocs(j,3),24,'red','filled'); hold off;        % pre. output
            hold on; scatter3(postclocs(j,1),postclocs(j,2),postclocs(j,3),24,'blue','filled'); hold off; % post. input
            disp(['    pre: ' num2str(double(prelocs(j,:))./conf.swcSize) ',  post: ' num2str(double(postclocs(j,:))./conf.swcSize)]);
        end
%}
%%{
        [md, didx] = sort(rcpostDs);
        for k=1:length(didx)
            if k>3, break; end
            j = didx(k);
            rcpostsidx = rcpostSidx{i}(j);
            rcnidx = preNidx(rcpostsidx);
            prelogi2 = ismember(preNidx,rcnidx);
            postlogi2 = ismember(postNidx,rcnidx);
            disp([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') dist=' num2str(rcpostDs(j)) 'nm']);

            % get all connected synapses
            loc1 = Spreloc(prelogi & postlogi2,:);
            loc2 = Spostloc(postlogi & prelogi2,:);

            % plot neuron and synapses
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(rcnidx)) '.swc'], true);
            figure;
            plotSwc(swc, [0.7 0.7 1], 1, true); view(3); grid on; axis image;
            title([num2str(i) '-' num2str(k) ' nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx)) ' (' tlabels{Ntype(i)+1} ',' tlabels{Ntype(rcnidx)+1} ') dist=' num2str(rcpostDs(j)) 'nm']);
            xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
            hold on; plotSwc(swc2, [1 0.8 0.8], 0.2, true); hold off; alpha(.4);
            hold on; scatter3(loc1(:,1),loc1(:,2),loc1(:,3),18,'red'); hold off;
            hold on; scatter3(loc2(:,1),loc2(:,2),loc2(:,3),18,'blue'); hold off;
            hold on; scatter3(preclocs(j,1),preclocs(j,2),preclocs(j,3),24,'red','filled'); hold off;  % pre. output
            hold on; scatter3(postlocs(j,1),postlocs(j,2),postlocs(j,3),24,'blue','filled'); hold off; % post. input
            disp(['    pre: ' num2str(double(preclocs(j,:))./conf.swcSize) ',  post: ' num2str(double(postlocs(j,:))./conf.swcSize)]);
        end
%}
        % plot all neurons and synapses related reciprocal (faster line plot). This one is heavy.
%%{
        figure;
        for k=1:length(rcNidx{i})
            swc2 = loadSwc([conf.swcPath '/' num2str(Nid(rcNidx{i}(k))) '.swc']);
            hold on; plotSwc(swc2, [0.8 0.8 0.8], 0.1); alpha(.1); hold off;
        end
        hold on; plotSwc(swc, [0.5 0.5 1], 1.5); hold off; view(3); grid on; axis image;
        title(['reciprocal connections(' num2str(i) ') nids=' num2str(Nid(i)) ',' num2str(Nid(rcnidx))]);
        xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
        hold on; scatter3(prelocs(:,1),prelocs(:,2),prelocs(:,3),18,'red'); hold off;
        hold on; scatter3(postlocs(:,1),postlocs(:,2),postlocs(:,3),18,'blue'); hold off;
%}
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

