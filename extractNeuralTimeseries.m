% make pre-processed time-series from registered calcium imaging data for neural functional connectivity.
% adding smooth and highpass filter.

function extractNeuralTimeseries
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = [0]; % high-pass filter threshold
    smooth = {'s230'};
    nuisance = {'poltcomp'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    extractPreprocTimeseries(preproc, hpfTh, smooth, nuisance)
end

function extractPreprocTimeseries(preproc, hpfTh, smooth, nuisance)
    TR = 1 / 1.879714;

    % load whole brain mask
    mV = niftiread('template/thresholded_FDACal_mask.nii.gz'); % mask should have same transform with 4D nifti data
    maskidx = find(mV>0);

    % for Nuisance Signal Regression
    csfV = []; % fly doesn't have csf
    wmF = 'template/jrc2018f_IBN_fiber_bundle_mirror_maskCal_invFDACal.nii.gz';
    wminfo = niftiinfo(wmF);
    wmV = niftiread(wminfo); % mask should have same transform with 4D nifti data
    gsF = 'template/thresholded_FDACal_mask.nii.gz';
    gsinfo = niftiinfo(gsF);
    gsV = niftiread(gsinfo); % mask should have same transform with 4D nifti data
    gsV(gsV>=1) = 1;
    gsV(gsV<1) = 0;

    if ~exist('results/neuralts','dir'), mkdir('results/neuralts'); end

    % read nii files
    listing = dir(['registered/' preproc '*green_FD_Warped.nii.gz']);
    for i=1:length(listing)
        folder = listing(i).folder;
        name = listing(i).name;
        id = split(name,'.');

        % read nii file
        V = single(niftiread([folder '/' name]));
        V(isnan(V)) = 0;

        % read rp_*.txt (6 head motion parameters)
        lpre = length(preproc);
        rpf = ['registered/rp_' id{1}(lpre+1:lpre+17) '.txt'];
        if ~exist(rpf,'file')
            disp(['file not found (skipped) : ' rpf]);
            continue;
        end
        disp(['loading : ' rpf]);
        M = readmatrix(rpf);
        Md = [zeros(1,6); diff(M,1,1)];

        for h=1:length(hpfTh)
            hpfstr = '';
            if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
    
            for k=1:length(smooth)
                for n=1:length(nuisance)
                    % output file
                    outfname = ['results/neuralts/' smooth{k} hpfstr nuisance{n} preproc id{1}(lpre+1:lpre+17) '-ts.mat'];
                    if exist(outfname,'file'), continue; end

                    % smoothing
                    Vk = V;
                    if ~isempty(smooth{k}) > 0 && (strcmp(smooth{k}(1),'s') || strcmp(smooth{k}(1),'m'))
                        % gaussian filter
                        sz = str2double(smooth{k}(2:end));
                        FWHM = [(sz/10)/(2.45/2.28) sz/10 (sz/10)/(3.715/2.28)]; % voxel size;
                        sigma = FWHM / sqrt(8*log(2));
                        filterSize = 2*ceil(2*sigma)+1;
    
                        disp(['sz=' num2str(sz) ', sigma=' num2str(sigma(1)) ', flSz=' num2str(filterSize(1))]);
                        for t=1:size(Vk,4)
                            Vt = imgaussfilt3(Vk(:,:,:,t), sigma, 'FilterSize', filterSize);
                            if strcmp(smooth{k}(1),'s')
                                Vk(:,:,:,t) = Vt;
                            else
                                Vk(:,:,:,t) = max(Vt,Vk(:,:,:,t));
                            end
                        end
                    end

                    % nuisance regression
                    Xn = []; perm = []; RiQ = []; dR2i = [];
                    if ~isempty(nuisance{n})
                        hm = 0;
                        if length(nuisance{n}) >= 3 && strcmp(nuisance{n}(1:3), '6hm')
                            hm = 6;
                            nuistr = nuisance{n}(4:end);
                        elseif length(nuisance{n}) >= 4 && strcmp(nuisance{n}(1:4), '24hm')
                            hm = 24;
                            nuistr = nuisance{n}(5:end);
                        else
                            nuistr = nuisance{n};
                        end
                        if strcmp(nuistr, 'nui')
                            % get Nuisance time-series (Global Mean, Global Signal, WM)
                            Xn = getNuisanceMeanTimeSeries(Vk, csfV, wmV, gsV);
                        elseif strcmp(nuistr, 'gm')
                            % get Nuisance time-series (Global Mean)
                            Xn = getNuisanceMeanTimeSeries(Vk, [], [], []);
                        elseif strcmp(nuistr, 'gmgs')
                            % get Nuisance time-series (Global Mean, Global Signal)
                            Xn = getNuisanceMeanTimeSeries(Vk, [], [], gsV);
                        elseif strcmp(nuistr, 'gmgsacomp')
                            % get Nuisance time-series (Global Mean, Global Signal, CSF comps, WM comps)
                            Sd = getNuisanceMeanTimeSeries(Vk, [], [], gsV);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, Sd);
                            Xn = [Sd, aComp];
                        elseif strcmp(nuistr, 'gmacomp')
                            % get Nuisance time-series (Global Mean, CSF comps, WM comps)
                            Sd = getNuisanceMeanTimeSeries(Vk, [], [], []);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, Sd);
                            Xn = [Sd, aComp];
                        elseif strcmp(nuistr, 'acomp')
                            % get Nuisance time-series (CSF comps, WM comps)
                            Xn = getNuisanceaCompCor(Vk, csfV, wmV);
                        elseif strcmp(nuistr, 'tcomp')
                            % get Nuisance time-series (high tSTD comps)
                            Xn = getNuisancetCompCor(Vk);
                        elseif strcmp(nuistr, 'tacomp')
                            % get Nuisance time-series (high tSTD comps, CSF comps, WM comps)
                            tComp = getNuisancetCompCor(Vk);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, tComp);
                            Xn = [tComp, aComp];
                        elseif strcmp(nuistr, 'pol')
                            % get Nuisance time-series (polynomial)
                            Xn = getNuisancePolynomial(size(Vk,4));
                        elseif strcmp(nuistr, 'polacomp')
                            % get Nuisance time-series (polynomial, CSF comps, WM comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, Sd);
                            Xn = [Sd, aComp];
                        elseif strcmp(nuistr, 'polacompn')
                            % get Nuisance time-series (polynomial, CSF comps, WM comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, Sd, 99, 6, true);
                            Xn = [Sd, aComp];
                        elseif strcmp(nuistr, 'poltcomp')
                            % get Nuisance time-series (polynomial, high tSTD comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            tComp = getNuisancetCompCor(Vk, Sd);
                            Xn = [Sd, tComp];
                        elseif strcmp(nuistr, 'poltcompn')
                            % get Nuisance time-series (polynomial, high tSTD comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            tComp = getNuisancetCompCor(Vk, Sd, 1000, 6, true);
                            Xn = [Sd, tComp];
                        elseif strcmp(nuistr, 'poltacomp')
                            % get Nuisance time-series (polynomial, high tSTD comps, CSF comps, WM comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            tComp = getNuisancetCompCor(Vk, Sd);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, [Sd, tComp]);
                            Xn = [Sd, tComp, aComp];
                        elseif strcmp(nuistr, 'poltacompn')
                            % get Nuisance time-series (polynomial, high tSTD comps, CSF comps, WM comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            tComp = getNuisancetCompCor(Vk, Sd, 1000, 6, true);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, [Sd, tComp], 99, 6, true);
                            Xn = [Sd, tComp, aComp];
                        elseif strcmp(nuistr, 'polgmtacomp')
                            % get Nuisance time-series (polynomial, high tSTD comps, CSF comps, WM comps)
                            Sd = getNuisancePolynomial(size(Vk,4));
                            Gm = getNuisanceMeanTimeSeries(Vk, [], [], []);
                            tComp = getNuisancetCompCor(Vk, [Sd, Gm]);
                            aComp = getNuisanceaCompCor(Vk, csfV, wmV, [Sd, Gm, tComp]);
                            Xn = [Sd, Gm, tComp, aComp];
                        end
                        if hm == 6
                            Xn = [Xn, M];
                        elseif hm == 24
                            Xn = [Xn, M, Md, M.^2, Md.^2];
                        end
                        [~, ~, perm, RiQ, dR2i] = regressPrepare(Xn);
                    end

                    sz = size(Vk);
                    Vk = reshape(Vk, [sz(1)*sz(2)*sz(3) sz(4)]);
                    Vk = Vk(maskidx,:); % voxel time-series
                    roinum = length(maskidx);

                    CY = cell(1,roinum);
                    parfor j=1:roinum
%                    for j=1:roinum
                        Z = double(Vk(j,:));
                        % nuisance regression out
                        if ~isempty(Xn)
                            Z = Z - nanmean(Z,2);
                            [~, R] = regressLinear(Z', Xn, [], [], perm, RiQ, dR2i); % 1st step OLS regression
                            Z = R';
                        end
                        % high pass filter as preprocessing step (M.W.Woolrich, 2001) type.
                        if hpfTh(h) > 0
                            Z = Z - nanmean(Z,2);
                            Z = highpass(Z',hpfTh(h),1/TR);
                            Z = Z';
                        end
                        CY{j} = single(Z);
                    end
                    X = nan(roinum,sz(4),'single');
                    for j=1:length(CY)
                        X(j,:) = CY{j};
                    end
                    clear CY;

                    Xm = nanmean(X,2);
                    X = X - Xm;
        
                    % show timeseries
%                    figure; plot(X'); title(name);

                    % save time-series
                    save(outfname,'X','Xm','-v7.3');
                end
            end
        end
    end
end
