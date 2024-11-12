% GLM for drosophila movement behavior
% multi-subject. mixed-effects with Tukey-Taper.
function glmFlyMB2nd
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = 0; % high-pass filter threshold
    smooth = 's40';
    nuisance = '';

    tuM = 8; % tukey window size
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hpfstr = '';
    if hpfTh > 0, hpfstr = ['hf' num2str(round(1/hpfTh))]; end

    % load background nii
    backNii = 'data/thresholded_FDACal.nii.gz';
    backinfo = niftiinfo(backNii);
    backV = niftiread(backinfo);

    % load mask nii
    maskinfo = niftiinfo('thresholded_FDACal_mask.nii.gz');
    maskV = niftiread(maskinfo);

    % read calcium image nii files

    % contrast image params
    contnames = {'movement'};
    contrasts = {[1 0]'};
    Pth = 0.001; % pvalue threshold
    rangePlus = [nan 15];
    rangeMinus = [nan 15];

    % calc 2nd-level estimation
    B1 = [];
    X2 = [];
    FWHMs = [];
    % read calcium image nii files
    listing = dir(['registered/' preproc '*green_FD_Warped.nii.gz']);
    for i=1:length(listing)
        name = listing(i).name;
        subject = name(7:13);

        betaBmat = [path smooth hpfstr nuisance preproc subject '-Tukey' num2str(tuM) '.mat'];
        if ~exist(betaBmat,'file')
            disp(['file not found. please calc individual sessions first : ' betaBmat])
            continue;
        end

        % load beta volumes
        f = load(betaBmat);
        % 2nd-level Y vector
        B2 = f.B2(:,[1,2]); % include design and intercept (we need more than 8 length for tukey taper)
        B1 = [B1; B2'];
        FWHMs = [FWHMs; f.FWHM];

        % 2nd-level design matrix
        X2 = [X2; eye(size(B2,2))];
    end
    B1(isnan(B1)) = 0; % there might be nan
    FWHMs = mean(FWHMs,1); % let's take the mean of FWHM.

    for tuM = tuM
        betaBmat = [path smooth hpfstr nuisance preproc subject '-Tukey' num2str(tuM) 'full.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            % calc 2nd-level estimation
            [B, RSS, df, X2is, tRs, R] = calcGlmTukey(B1, X2, tuM);

            [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, atlasV);

            % output beta matrix
            save(betaBmat,'B','RSS','X2is','tRs','recel','FWHM','df','-v7.3');
        end

        % GLM contrast images
        Ts = calcGlmContrastImage(contrasts, B, RSS, X2is, tRs);

        % GLM contrast image
        thParam = {df, Pth};
        clParam = {69, FWHMs}; % clustering parameter for GLM contrast
%        thParam = {df, 0.01};
%        clParam = {64, FWHMs}; % clustering parameter for GLM contrast
        [Tth, Vts, Vfs, Tmaxs, Tcnts] = plotGlmContrastImage(contnames, Ts, thParam, clParam, maskV, true, false, backV, ...
            ['glmFlyMB ' '2nd-mix-Tukey' num2str(tuM) 'full'], rangePlus, rangeMinus, [], [], []);

        % save T-value NIfTI volume
        saveContrastNii(backNii,contnames,Vts,path,[cubename prefix 'D_2nd-mix-Tukey' num2str(tuM) 'th']);
%        saveContrastNii(backNii,contnames,Vfs,path,[cubename prefix 'D_2nd-mix-Tukey' num2str(tuM) 'full']);
    end
end

%%
function saveContrastNii(tfmri, contnames, V2s, path, outname)
    info = niftiinfo(tfmri);
    info.ImageSize = info.ImageSize(1:3);
    info.PixelDimensions = info.PixelDimensions(1:3);
    info.raw.dim(1) = 3;
    info.raw.dim(5) = 1;
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    for j=1:length(contnames)
        fname = [path outname '_' contnames{j} '.nii'];
        niftiwrite(V2s{j},fname,info,'Compressed',true);
    end
end

