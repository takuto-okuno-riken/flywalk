% GLM for drosophila movement behavior
% single session. with Tukey 8 window
function glmFlyMB1st
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = 0; % high-pass filter threshold
    smooth = 's230';
    nuisance = 'poltcomp';

    TR = 1 / 1.879714;
    tuM = 8; % tukey window size
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hpfstr = '';
    if hpfTh > 0, hpfstr = ['hf' num2str(round(1/hpfTh))]; end

    path = 'results/glm/';
    if ~exist(path, 'dir')
        mkdir(path)
    end

    % load movement behavior data
    load('data/behavior.mat');

    % load background nii
    backinfo = niftiinfo('template/thresholded_FDACal.nii.gz');
    backV = niftiread(backinfo);

    % load mask nii
    maskinfo = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    maskV = niftiread(maskinfo);

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

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(18);

    % read calcium image nii files
    listing = dir(['registered/' preproc '*green_FD_Warped.nii.gz']);
    for i=1:length(listing)
        folder = listing(i).folder;
        name = listing(i).name;
        subject = name(7:13);

        % check session processed or not
        betaBmat = [path smooth hpfstr nuisance preproc subject '-Tukey' num2str(tuM) '.mat'];
        if exist(betaBmat,'file')
            disp(['file found : ' betaBmat]);
            checkTukeyRange([], [], path, [smooth hpfstr nuisance preproc], maskV, backV, subject, tuM)
            continue;
        end
        
        % read nii file
        niifile = [folder '/' name];
        disp(['loading : ' niifile]);
        V = single(niftiread(niifile));

        if ~isempty(smooth) > 0 && strcmp(smooth(1),'s')
            % gaussian filter
            sz = str2double(smooth(2:end));
            FWHM = [(sz/10)/(2.45/2.28) sz/10 (sz/10)/(3.715/2.28)]; % voxel size;
            sigma = FWHM / sqrt(8*log(2));
            filterSize = 2*ceil(2*sigma)+1;

            disp(['sz=' num2str(sz) ', sigma=' num2str(sigma(1)) ', flSz=' num2str(filterSize(1))]);
            for t=1:size(V,4)
                Vt = imgaussfilt3(V(:,:,:,t), sigma, 'FilterSize', filterSize);
                V(:,:,:,t) = Vt;
            end
        end

        % extract voxel time-series by mask
        disp('apply mask ...');
        aIdx = find(maskV(:) > 0);
        A = reshape(V,[],size(V,4));
        Z = A(aIdx,:);
        Z = Z - nanmean(Z,2);
        Z = Z';
        
        % high pass filter as preprocessing step (M.W.Woolrich, 2001) type.
        if hpfTh > 0
            disp(['apply highpass filter (' num2str(hpfTh) ' Hz) : tfMRI and design matrix...']);
            Z = highpass(Z,hpfTh,1/TR);
        end
        % apply nuisance removal
        Xn = []; nuistr = '';
        if ~isempty(nuisance)
            hm = 0;
            if length(nuisance) >= 3 && strcmp(nuisance(1:3), '6hm')
                hm = 6;
                nuistr = nuisance(4:end);
            elseif length(nuisance) >= 4 && strcmp(nuisance(1:4), '24hm')
                hm = 24;
                nuistr = nuisance(5:end);
            else
                nuistr = nuisance;
            end
            if strcmp(nuistr, 'nui')
                % get Nuisance time-series (Global Mean, Global Signal, WM)
                Xn = getNuisanceMeanTimeSeries(V, csfV, wmV, gsV);
            elseif strcmp(nuistr, 'gm')
                % get Nuisance time-series (Global Mean)
                Xn = getNuisanceMeanTimeSeries(V, [], [], []);
            elseif strcmp(nuistr, 'gmgs')
                % get Nuisance time-series (Global Mean, Global Signal)
                Xn = getNuisanceMeanTimeSeries(V, [], [], gsV);
            elseif strcmp(nuistr, 'poltcomp')
                % get Nuisance time-series (polynomial, high tSTD comps)
                Sd = getNuisancePolynomial(size(V,4));
                tComp = getNuisancetCompCor(V, Sd);
                Xn = [Sd, tComp];
            end
            if hm == 6
                Xn = [Xn, M];
            elseif hm == 24
                Xn = [Xn, M, Md, M.^2, Md.^2];
            end
        end

        % get design matrix
        X = MBcal(i,:)';
        % K*W*Y = K*W*B*X + K*W*e case. (SPM type. K.J.Friston, 2000)
        % actually, FLS (FEAT) do this as default option
        if hpfTh > 0
            X = highpass(X,hpfTh,1/TR);
        end
        if length(nuistr)>=3 && strcmp(nuistr(1:3),'pol')
            X = [X Xn]; % Nuisance should be raw
        else
            X = [X Xn ones(size(X,1),1)]; % Nuisance should be raw
        end

        figure; imagesc(X, [-0.1, 1.1]); colorbar;

        % check Tukey range
        checkTukeyRange(Z, X, path, [smooth hpfstr nuisance preproc], maskV, backV, subject, tuM);
    end
end


function checkTukeyRange(Z, X, path, prefix, maskV, tempV, subject, tuM)
    % contrast image params
    contnames = {'movement'};
    contrasts = {}; % no nuisanse
    Pth = 0.001; % pvalue threshold

    % ---------------------------------------------------------------------
    % AR estimation with Tukey-Taper of frequency domain
    for tuM = tuM
        betaBmat = [path prefix subject '-Tukey' num2str(tuM) '.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            [B2, RSS, df, X2is, tRs, R] = calcGlmTukey(Z, X, tuM);

            [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, maskV);

            % output beta matrix
            save(betaBmat,'B2','RSS','X2is','tRs','recel','FWHM','df','-v7.3');
        end
    
        % GLM contrast images
        contrasts{1} = zeros(size(B2,2),1); contrasts{1}(1) = 1;
        Ts = calcGlmContrastImage(contrasts, B2, RSS, X2is, tRs);

        % plot contrast image
        thParam = {df, Pth};
        clParam = {60, FWHM}; % clustering parameter for GLM contrast
        plotGlmContrastImage(contnames, Ts, thParam, clParam, maskV, true, false, tempV, ...
            ['GLM ' prefix subject 'CTukey' num2str(tuM)]);
    end
end

