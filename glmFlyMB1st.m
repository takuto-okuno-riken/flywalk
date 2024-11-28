% GLM for drosophila movement behavior
% single session. with Tukey 8 window
function glmFlyMB1st
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = 0; % high-pass filter threshold
    smooth = 's40';
    nuisance = '';

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

    % set pool num. this calculation takes time. we need big pool num.
    delete(gcp('nocreate')); % shutdown pools
    parpool(32);

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
            continue;
        end
        
        % read nii file
        niifile = [folder '/' name];
        disp(['loading : ' niifile]);
        V = single(niftiread(niifile));

        if ~isempty(smooth) > 0 && strcmp(smooth(1),'s')
            % gaussian filter
            sz = str2double(smooth(2:3));
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
        Xn = [];
        if ~isempty(nuisance)
%{
            % get Nuisance time-series (CSF, WM, Global Signal, Global Mean)
            Xn = getNuisanceMeanTimeSeries(V, csfV, wmV, gsV);
            % get Nuisance time-series (Global Mean, CSF comps, WM comps)
            Sd = getNuisanceMeanTimeSeries(V, [], [], []);
            aComp = getNuisanceaCompCor(V, csfV, wmV, Sd);
            Xn = [Sd, aComp];
%}
        end

        % get design matrix
        X = MBcal(i,:)';
        % K*W*Y = K*W*B*X + K*W*e case. (SPM type. K.J.Friston, 2000)
        % actually, FLS (FEAT) do this as default option
        if hpfTh > 0
            X = highpass(X,hpfTh,1/TR);
        end
        X = [X Xn]; % Nuisance should be raw

        figure; imagesc([X ones(size(X,1),1)], [-0.1, 1.1]); colorbar;

        % check Tukey range
        checkTukeyRange(Z, X, path, [smooth hpfstr nuisance preproc], maskV, backV, subject, tuM);
    end
end


function checkTukeyRange(Z, X, path, prefix, maskV, tempV, subject, tuM)
    % contrast image params
    contnames = {'movement'};
    contrasts = {[1 0]'}; % no nuisanse
    Pth = 0.001; % pvalue threshold

    % ---------------------------------------------------------------------
    % AR estimation with Tukey-Taper of frequency domain
    for tuM = tuM
        betaBmat = [path prefix subject '-Tukey' num2str(tuM) '.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            Xt = [X, ones(size(X,1),1)];
            [B2, RSS, df, X2is, tRs, R] = calcGlmTukey(Z, Xt, tuM);

            [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, maskV);

            % output beta matrix
            save(betaBmat,'B2','RSS','X2is','tRs','recel','FWHM','df','-v7.3');
        end
    
        % GLM contrast images
        Ts = calcGlmContrastImage(contrasts, B2, RSS, X2is, tRs);

        % plot contrast image
        thParam = {df, Pth};
        clParam = {60, FWHM}; % clustering parameter for GLM contrast
        plotGlmContrastImage(contnames, Ts, thParam, clParam, maskV, true, false, tempV, ...
            ['GLM ' prefix subject 'CTukey' num2str(tuM)]);
    end
end

