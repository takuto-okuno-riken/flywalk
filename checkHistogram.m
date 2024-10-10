% checking registered calcium imaging histogram, etc.

function checkHistogram
    % load FDA cal mask
    mV = niftiread('data/thresholded_FDACal_mask.nii.gz');
    midx = find(mV>0);

    % read nii files
    listing = dir('registered/ar*green_FD_Warped.nii.gz');
    M = []; S = []; names = {};
    for i=1:length(listing)
        folder = listing(i).folder;
        name = listing(i).name;
        id = split(name,'.');

        % read nii file
        info = niftiinfo([folder '/' name]);
        V = single(niftiread(info));
        sz = size(V);
        V = reshape(V, [sz(1)*sz(2)*sz(3) sz(4)]);

        % thresholding to reduce image size
        Vi = V(midx,:);
        m = mean(Vi,'all');
        s = std(Vi,1,'all');
        M = [M, m];
        S = [S, s];
        names{i} = id{1};
        disp([id{1} ' : m=' num2str(m) ', s=' num2str(s)]);

        % show histogram
        figure; histogram(Vi); title(name);
    end
    n = length(midx);
    save('results/thresholded_FDACal_histogram.mat','M','S','names','n');
end
