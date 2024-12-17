% extract ROI time-series from registered calcium imaging data.
% adding smooth and highpass filter
% this script should run after makeVoxelROIatlas.m, makeStructConnectivity.m

function extractROItimeseries
    %%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-process
    preproc = 'ar'; % for move correct, slice time correct
%    preproc = 'r'; % for move correct only

    % output time-series (smoothing, highpass filter, nuisance removal)
    hpfTh = [0]; % high-pass filter threshold
%    hpfTh = [0, 0.1, 0.05, 0.025, 0.02, 0.01, 0.009, 0.008, 0.005, 0.001]; % high-pass filter threshold
%    smooth = {'', 's10', 's20', 's30', 's40', 's50', 's60', 's70', 's80'};
%    smooth = {'s90', 's100', 's110', 's120', 's130', 's140', 's150', 's160', 's170', 's180', 's190', 's200', 's210', 's220', 's230', 's240', 's250', 's260', 's270', 's280', 's290', 's300'};
%    smooth = {'s30', 's80'};
     smooth = {'s150','s230','s300'};
%    smooth = {''};
    nuisance = {'','gm','gmgs','nui','6hm','6hmgm','6hmgmgs','6hmnui','24hm','24hmgm','24hmgmgs','24hmnui', ... %12
        'acomp','gmacomp','gmgsacomp','tcomp','tacomp', ... %17
        '6hmacomp','6hmgmacomp','6hmgmgsacomp','6hmtcomp','6hmtacomp', ... %22
        '24hmacomp','24hmgmacomp','24hmgmgsacomp','24hmtcomp','24hmtacomp', ... %27
        'pol','polacomp','poltcomp','poltacomp','polgmtacomp', ...
        '6hmpol','6hmpolacomp','6hmpoltcomp','6hmpoltacomp','6hmpolgmtacomp', };
%    nuisance = {'6hmtacomp'}; % good for bransonhemi, branson7065km50
    nuisance = {'','poltcomp'}; % good for DistKm(synapse)
%    nuisance = {''};

    % ROI name
%    roitypes = {'hemiroi','bransonhemi'};
%    roitypes = {'hemiBranson7065'};
%    roitypes = {'hemiBranson7065km20','hemiBranson7065km30','hemiBranson7065km50','hemiBranson7065km100','hemiBranson7065km200', ...
%        'hemiBranson7065km300','hemiBranson7065km500','hemiBranson7065km1000'};
    roitypes = {'wireBranson7065km20','wireBranson7065km30','wireBranson7065km50','wireBranson7065km100','wireBranson7065km200', ...
        'wireBranson7065km300','wireBranson7065km500','wireBranson7065km1000'};
%    roitypes = {'hemiCmkm20','hemiCmkm30','hemiCmkm50','hemiCmkm100','hemiCmkm200', ...
%        'hemiCmkm300','hemiCmkm500','hemiCmkm1000'};
%    roitypes = {'hemiCmkm20r1w1','hemiCmkm30r1w1','hemiCmkm50r1w1','hemiCmkm100r1w1','hemiCmkm200r1w1', ...
%        'hemiCmkm300r1w1','hemiCmkm500r1w1','hemiCmkm1000r1w1'};
%    roitypes = {'hemiDistKm20','hemiDistKm30','hemiDistKm50','hemiDistKm100','hemiDistKm200', ...
%        'hemiDistKm300','hemiDistKm500','hemiDistKm1000'};
%    roitypes = {'hemiRand20','hemiRand30','hemiRand50','hemiRand100','hemiRand200', ...
%        'hemiRand300','hemiRand500','hemiRand1000'};
%    roitypes = {'hemiVrand20','hemiVrand30','hemiVrand50','hemiVrand100','hemiVrand200', ...
%        'hemiVrand300','hemiVrand500','hemiVrand1000'};
%    roitypes = {'hemiCube12','hemiCube8','hemiCube4'};
%    roitypes = {'hemiPiece12','hemiPiece8','hemiPiece4'};
%    roitypes = {'hemiPiece3','hemiPiece2'};
    % neuropil FB, EB, EB-bL(L), bL-b'L-aL-a'L-BU(L)
%    roitypes = {'hemiRoi101','hemiRoi57','hemiRoi57-51','hemiRoi51-62-20-111-100'};
%    roitypes = {'hemiRoi1','hemiRoi5','hemiRoi7','hemiRoi27','hemiRoi30','hemiRoi32','hemiRoi43','hemiRoi52', ...
%        'hemiRoi54','hemiRoi57','hemiRoi59','hemiRoi63','hemiRoi65','hemiRoi67','hemiRoi78','hemiRoi82', ...
%        'hemiRoi89','hemiRoi93','hemiRoi95','hemiRoi100','hemiRoi101','hemiRoi106','hemiRoi113'};
%    roitypes = {'hemiroi','hemiBranson7065km50','hemiCmkm50','hemiCmkm50r1w1','hemiDistKm50','hemiRand50','hemiVrand50'};
%    roitypes = {'hemiBranson7065km200','hemiCmkm200','hemiCmkm200r1w1','hemiDistKm200','hemiRand200','hemiVrand200'};
    roitypes = {'hemiCmkm5000','hemiDistKm5000','hemiCmkm10000','hemiDistKm10000'};
%    roitypes = {'hemiCmkm50','hemiDistKm50','hemiCmkm100','hemiDistKm100','hemiCmkm500','hemiDistKm500'}; % for large smoothing size & no nuisanse, poltcomp
%    roitypes = {'hemiroi','hemiDistKm50'}; % for all nuisanse & s30, s80

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    extractTsROItype(roitypes, preproc, hpfTh, smooth, nuisance)
end

function roiIdxs = getRoiIdxs(roitype)
    roiIdxs = {};
    switch(roitype)
    case 'hemiroi'
        listing = dir(['atlas/' roitype '/*.nii.gz']);
        for i=1:length(listing)
            V = niftiread(['atlas/' roitype '/roi' num2str(i) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
            roiIdxs{i} = find(V>0);
        end
    case 'bransonhemi'
        % read hemibrain mask
        Vm = niftiread('template/jrc2018f_flyemhemibrainCal_invFDACal.nii.gz');
        Vm(Vm>0) = 1;
        Vm(Vm<1) = 0;
        V = niftiread('atlas/JRC2018_branson_atlasCal_invFDACal.nii.gz'); % ROI mask should have same transform with 4D nifti data
        V = V .* Vm;
        roimax = max(V(:));
        for i=1:roimax
            roiIdxs{i} = find(V==i);
        end
    otherwise
        V = niftiread(['atlas/' roitype 'atlasCal.nii.gz']); % ROI mask should have same transform with 4D nifti data
        roitype = lower(roitype);
        roimax = max(V(:));
        for i=1:roimax
            roiIdxs{i} = find(V==i);
        end
    end
end

function extractTsROItype(roitypes, preproc, hpfTh, smooth, nuisance)
    TR = 1 / 1.879714;

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

    if ~exist('results/ts','dir'), mkdir('results/ts'); end

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

        for ii=1:length(roitypes)
            disp(['get ROI index : ' roitypes{ii}]);
            roiIdxs = getRoiIdxs(roitypes{ii});
            
            for h=1:length(hpfTh)
                hpfstr = '';
                if hpfTh(h) > 0, hpfstr = ['hf' num2str(round(1/hpfTh(h)))]; end
        
                for k=1:length(smooth)
                    for n=1:length(nuisance)
                        CX = {}; CXm = {};
                        % load file
                        outfname = ['results/ts/' smooth{k} hpfstr nuisance{n} preproc roitypes{ii} '-ts.mat'];
                        if exist(outfname,'file')
                            load(outfname);
                        end
                        if length(CX) < i
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
                                elseif strcmp(nuistr, 'poltcomp')
                                    % get Nuisance time-series (polynomial, high tSTD comps)
                                    Sd = getNuisancePolynomial(size(Vk,4));
                                    tComp = getNuisancetCompCor(Vk, Sd);
                                    Xn = [Sd, tComp];
                                elseif strcmp(nuistr, 'poltacomp')
                                    % get Nuisance time-series (polynomial, high tSTD comps, CSF comps, WM comps)
                                    Sd = getNuisancePolynomial(size(Vk,4));
                                    tComp = getNuisancetCompCor(Vk, Sd);
                                    aComp = getNuisanceaCompCor(Vk, csfV, wmV, [Sd, tComp]);
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

                            roinum = length(roiIdxs);
                            CZ = cell(1,roinum);
                            for j=1:roinum
                                midx = roiIdxs{j};
                                if ~isempty(midx)
                                    CZ{j} = Vk(midx,:);
                                end
                            end
                            clear Vk;

                            CY = cell(1,roinum);
                            parfor j=1:roinum
%                            for j=1:roinum
                                if ~isempty(CZ{j})
                                    Z = CZ{j};
                                    % nuisance regression out
                                    if ~isempty(Xn)
                                        disp(['apply nuisance regression out (' nuisance{n} ') : ROI(' num2str(j) ') time-series ...']);
                                        Z = Z - nanmean(Z,2);
                                        for m=1:size(Z,1)
                                            [~, R] = regressLinear(Z(m,:)', Xn, [], [], perm, RiQ, dR2i); % 1st step OLS regression
                                            Z(m,:) = R;
                                        end
                                    end
                                    % high pass filter as preprocessing step (M.W.Woolrich, 2001) type.
                                    if hpfTh(h) > 0
                                        disp(['apply highpass filter (' num2str(hpfTh(h)) ' Hz) : ROI(' num2str(j) ') time-series ...']);
                                        Z = Z - nanmean(Z,2);
                                        Z = highpass(Z',hpfTh(h),1/TR);
                                        Z = Z';
                                    end
                                    x = nanmean(Z,1);
                                else
                                    x = nan(1,sz(4));
                                end
                                CY{j} = x;
                            end
                            X = nan(length(CY),sz(4));
                            for j=1:length(CY)
                                X(j,:) = CY{j};
                            end
                            clear CY;

                            Xm = nanmean(X,2);
                            X = X - Xm;
    
                            CX{i} = X;
                            CXm{i} = Xm;
                
                            % show timeseries
    %                        figure; plot(X'); title(name);

                            % save time-series
                            save(outfname,'CX','CXm','-v7.3');
                        end
                    end
                end
            end
        end
    end
end
