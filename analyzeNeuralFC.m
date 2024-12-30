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

%    scTh = 80; synTh = 0; % FlyEM synapse confidence & synapse count at one neuron threshold
    scTh = 60; synTh = 5; % almost flywire codex compatible setting
    analyzeNeuralFc(preproc, hpfTh, smooth, nuisance, sbjids, 'hemi', synTh, scTh); % FlyEM hemi brain.

%    scTh = 130; synTh = 0; % FlyWire synapse score & synapse count at one neuron threshold
    scTh = 50; synTh = 0;
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
                Dmz = cell(length(inIdx),1);
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
                    D3 = zeros(length(minidx),length(moutidx),cmlen,'single');
                    for j=1:cmlen
                        D3(:,:,j) = corr(CX{j}(minidx,:)',CX{j}(moutidx,:)');
                    end
                    D3z = atanh(D3); % z transformed (better FC-SC corr).
                    Dmz{i} = nanmean(D3z,3);
%                    figure; imagesc(Dmz{i}); colorbar;
                end
                save(outfname, 'Dmz', '-v7.3');
            end
        end
    end
end
