% analyze and plot Neural FC result.
% this script can run after analyzeNeuralFC.m

function plotNeuralFC
    % DBscan param
    epsilon = 5; % micro meter. almost 2 voxels.
    minpts = 3; % set 1, but isolated synapse will be ignored

    % check NeuralFC of FlyEM
    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
%    scTh = 60; synTh = 5; % almost flywire codex compatible setting
    checkNeuralFC(synTh, scTh, epsilon, minpts);

    % check NeuralFC of FlyWire
    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
%    scTh = 50; synTh = 0;
%    scTh = 50; synTh = 5; % for checking flywire codex compatible
    checkNeuralFCFw(synTh, scTh, epsilon, minpts);
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
                    if ~isempty(Dmz{i})
                        type = ''; % neuro transmitter type is unknown
                        F = Dmz{i}(:); idx = find(~isinf(F)); F = F(idx);
                        Dt = D{i}(:); Dt = sqrt(Dt(idx));
                        r = corr(F,Dt);
                        disp([num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type ' in=' num2str(sum(inCount{i})) ' (' num2str(size(Dmz{i},1)) ...
                            ')  out=' num2str(sum(outCount{i})) ' ('  num2str(size(Dmz{i},2)) ')']);

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
                        figure; imagesc(Dmz{i}(si1,si2)); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
%                        figure; imagesc(D{i}(si1,si2)); colorbar; title(['Sy Dist ' num2str(synTh) 'sr' num2str(confTh) ' (' num2str(i) ') nid=' num2str(tNid(i)) ' ' type]);

                        % calc DBscan synaptic cluster
                        CF = calcDBscanSycluster(Dmz{i}, maxcls, Cidx1, Cidx2);
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

function checkNeuralFCFw(synTh, scoreTh, epsilon, minpts)
    tlabels = {'Unknown','DA','SER','GABA','GLUT','ACH','OCT'};
    preproc = 'ar'; % for move correct, slice time correct
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};

    info = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    V = niftiread(info);
    sz = size(V);

    % FlyWire read neuron info
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read neural SC
    load(['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutVoxels.mat']);
    load(['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralInOutDistance.mat']);
    load(['results/neuralsc/wire' num2str(synTh) 'sr' num2str(scoreTh) '_neuralDBScan' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

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
                    if ~isempty(Dmz{i})
                        t = Ntype(i)+1;
                        type = tlabels{t};
                        F = Dmz{i}(:); idx = find(~isinf(F)); F = F(idx);
                        Dt = D{i}(:); Dt = sqrt(Dt(idx));
                        r = corr(F,Dt);
                        disp([num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type ' in=' num2str(sum(inCount{i})) ' (' num2str(size(Dmz{i},1)) ...
                            ')  out=' num2str(sum(outCount{i})) ' ('  num2str(size(Dmz{i},2)) ')']);

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
                        figure; imagesc(Dmz{i}(si1,si2)); colorbar; title(['Sy FC '  num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type]);
                        xlabel('Synaptic cluster (output)'); ylabel('Synaptic cluster (input)');
%                        figure; imagesc(D{i}(si1,si2)); colorbar; title(['Sy Dist ' num2str(synTh) 'sr' num2str(scoreTh) ' (' num2str(i) ') nid=' num2str(Nid(i)) ' ' type]);

                        % calc DBscan synaptic cluster
                        CF = calcDBscanSycluster(Dmz{i}, maxcls, Cidx1, Cidx2);
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

