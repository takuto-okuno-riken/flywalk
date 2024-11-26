% make voxel ROI atlas. Some types of atlas, such as cube ROI type, piece ROI type and
% neuropil specific voxel (ROI) type.

function makeVoxelROIatlas
    name = 'hemi';

    % read hemibrain mask
    info = niftiinfo('data/jrc2018f_flyemhemibrainCal_invFDACal.nii.gz');
    Vm = niftiread(info);
    Vm(Vm>0) = 1;
    Vm(Vm<1) = 0;

    % read FDA mask
    info = niftiinfo('data/thresholded_FDACal_mask.nii.gz');
    Vfm = niftiread(info);
    Vfm(Vfm>0) = 1;
    Vfm(Vfm<1) = 0;
    Vm = Vm .* single(Vfm); % and condition.

    % make cube ROI type atlas
    for atlasSize = 12:-1:2
        cubename = [name 'Cube' num2str(atlasSize)];
        atlas = ['data/' cubename 'atlasCal.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            aV = niftiread(atlasinfo);
        else
            c = 1;
            sz = size(Vm);
            aV = zeros(sz(1),sz(2),sz(3),'single');
            for z=1:floor(sz(3)/atlasSize)
                zi = 1+(z-1)*atlasSize:z*atlasSize;
                for y=1:floor(sz(2)/atlasSize)
                    yi = 1+(y-1)*atlasSize:y*atlasSize;
                    for x=1:floor(sz(1)/atlasSize)
                        xi = 1+(x-1)*atlasSize:x*atlasSize;
                        Vi = Vm(xi,yi,zi);
                        if sum(Vi,'all') > 0
                            aV(xi,yi,zi) = c; c = c + 1;
                        end
                    end
                end
            end
            aV = aV .* Vm;

            % use info as template for cubeAtlas, but direction is opposite.
            if max(aV(:)) > 65535
                aV = int32(aV);
                info.Datatype = 'int32';
                info.BitsPerPixel = 32;
            else
                aV = uint16(aV);
                info.Datatype = 'uint16';
                info.BitsPerPixel = 16;
            end

            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'Cube ROI';
            % output nii file
            niftiwrite(aV,atlas,info,'Compressed',true);
        end

        disp([atlas ' ROI count=' num2str(max(aV(:)))]);
    end

    % make piece ROI type atlas (take one voxel per 2,3,4... cube voxel)
    for atlasSize = [12 8 4 3 2]
        piecename = [name 'Piece' num2str(atlasSize)];
        atlas = ['data/' piecename 'atlasCal.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            aV = niftiread(atlasinfo);
        else
            c = 1;
            sz = size(Vm);
            aV = zeros(sz(1),sz(2),sz(3),'single');
            for z=1:floor(sz(3)/atlasSize)
                zi = 1+(z-1)*atlasSize;
                for y=1:floor(sz(2)/atlasSize)
                    yi = 1+(y-1)*atlasSize;
                    for x=1:floor(sz(1)/atlasSize)
                        xi = 1+(x-1)*atlasSize;
                        Vi = Vm(xi,yi,zi);
                        if Vi > 0
                            aV(xi,yi,zi) = c; c = c + 1;
                        end
                    end
                end
            end
            aV = aV .* Vm;

            % use info as template for cubeAtlas, but direction is opposite.
            if max(aV(:)) > 65535
                aV = int32(aV);
                info.Datatype = 'int32';
                info.BitsPerPixel = 32;
            else
                aV = uint16(aV);
                info.Datatype = 'uint16';
                info.BitsPerPixel = 16;
            end

            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'Piece ROI';
            % output nii file
            niftiwrite(aV,atlas,info,'Compressed',true);
        end

        disp([atlas ' ROI count=' num2str(max(aV(:)))]);
    end

    % make branson supervoxel based atlas
    atlas = ['data/' name 'Branson7065atlasCal.nii' ];
    if exist([atlas '.gz'],'file')
        atlasinfo = niftiinfo([atlas '.gz']);
        V = niftiread(atlasinfo);
    else
        mV = niftiread('data/jrc2018f_flyemhemibrainCal_invFDACal.nii.gz');
        mV(mV>0) = 1; mV(mV<1) = 0;

        info = niftiinfo('data/jrc2018f_branson7065Cal_invFDACal.nii.gz');
        aV = niftiread(info); sz = size(aV);
        amV = aV .* mV;
        ids = unique(amV(:)); ids(ids==0) = [];
        V = zeros(sz(1),sz(2),sz(3),'uint16');
        roiid = 1;
        for i=1:length(ids)
            id = ids(i);
            aidx = find(aV==id);
            bidx = find(amV==id);
            if length(aidx) == length(bidx)  % hemiem ROI should contain whole supervoxel.
                V(bidx) = roiid;
                roiid = roiid + 1;
            end
        end
        info.Datatype = 'uint16';
        info.BitsPerPixel = 16;
        info.Description = 'branson 7065 supervoxels atlas';

        % output nii file
        niftiwrite(V,atlas,info,'Compressed',true);
    end
    disp([atlas ' ROI count=' num2str(max(V(:)))]);

    % make neuropil specific voxel type atlas
%    roiids = {[101],[57],[57,51],[51,62,20,111,100]}; % FB, EB, EB-bL(L),bL-b'L-aL-a'L-BU(L)
    roiids = {1	5	7	27	30	32	43	52	54	57	59	63	65	67	78	82	89	93	95	100	101	106	113};
    for i=1:length(roiids)
        idstr = num2str(roiids{i}(1));
        for j=2:length(roiids{i}), idstr=[idstr '-' num2str(roiids{i}(j))]; end
        roiname = [name 'Roi' idstr];
        atlas = ['data/' roiname 'atlasCal.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            aV = niftiread(atlasinfo);
        else
            info = niftiinfo(['data/flyemroi/roi' num2str(roiids{i}(1)) '.nii.gz']);
            aV = niftiread(info); % ROI mask should have same transform with 4D nifti data
            for j=2:length(roiids{i})
                bV = niftiread(['data/flyemroi/roi' num2str(roiids{i}(j)) '.nii.gz']);
                bidx = find(bV>0);
                aV(bidx) = bV(bidx);
            end
            idx = find(aV>0);
            aV(idx) = 1:length(idx);

            if max(aV(:)) > 65535
                aV = int32(aV);
                info.Datatype = 'int32';
                info.BitsPerPixel = 32;
            else
                aV = uint16(aV);
                info.Datatype = 'uint16';
                info.BitsPerPixel = 16;
            end

            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'neuropil specific atlas';
            % output nii file
            niftiwrite(aV,atlas,info,'Compressed',true);
        end
        disp([atlas ' ROI count=' num2str(max(aV(:)))]);
    end

    % make whole flyem ROI voxel atlas (except fibers)
    primaryId = [1	2	4	5	7	8	10	15	16	18	19 20	22	24	27	28	30	31 32	33	34	38	41	42 43	45	47	49	50	51	52	54	56	57	58	59	62	63	65	66	67	68	75	76 78	80	82	87	89	91	93	95	97	98	100	101	102	103	106	107	111	112	113];
    atlas = ['data/' name 'RoiWholeatlasCal.nii' ];
    if exist([atlas '.gz'],'file')
        atlasinfo = niftiinfo([atlas '.gz']);
        aV = niftiread(atlasinfo);
    else
        info = niftiinfo(['data/flyemroi/roi' num2str(primaryId(1)) '.nii.gz']);
        aV = niftiread(info); % ROI mask should have same transform with 4D nifti data
        for i=1:length(primaryId)
            bV = niftiread(['data/flyemroi/roi' num2str(primaryId(i)) '.nii.gz']);
            bidx = find(bV>0);
            aV(bidx) = primaryId(i);
        end
        % remove fibers
        fV = niftiread('data/jrc2018f_IBN_fiber_bundle_mirror_maskCal_invFDACal.nii.gz');
        aV(fV>127) = 0; 
        idx = find(aV>0);
        aV(idx) = 1:length(idx); % if comment out this line, ROI atlas can be saved

        if max(aV(:)) > 65535
            aV = int32(aV);
            info.Datatype = 'int32';
            info.BitsPerPixel = 32;
        else
            aV = uint16(aV);
            info.Datatype = 'uint16';
            info.BitsPerPixel = 16;
        end

        % set info. info.raw is not necessary to set (niftiwrite() does it)
        info.Description = 'neuropil specific atlas';
        % output nii file
        niftiwrite(aV,atlas,info,'Compressed',true);
    end
    disp([atlas ' ROI count=' num2str(max(aV(:)))]);


    % make ROI atlas based on k-means clustering of hemi branson7065 SC.
    % this requires hemibranson7065_connectlist.mat file. so need to run
    % makeStructConnectivity.m (Branson 7065) first.
    % this needs Statistics and Machine Learning Toolbox.
    for k=[20 30 50 100 200 300]
        atlas = ['data/' name 'Branson7065km' num2str(k) 'atlasCal.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            aV = niftiread(atlasinfo);
        else
            load('data/hemibranson7065_connectlist.mat');
            lcm2 = log10(countMat2(:,:,2));
            lcm2(lcm2<0) = 0;
            % to avoid zero for k-means
            Xnorm = sqrt(sum(lcm2.^2, 2));
            lcm2(Xnorm==0,1) = eps(max(Xnorm)) * 2;

            idx = kmeans(lcm2,k,'Distance','cosine');
%{
            % check sorted group result
            [S,I] = sort(idx,'ascend');
            figure; imagesc(log10(countMat2(:,:,2))); colorbar;
            figure; imagesc(log10(countMat2(I,I,2))); colorbar;

            M = zeros(k,k,'single');
            for i=1:k
                a = sum(countMat2(idx==i,:,2),1);
                for j=1:k
                    M(i,j) = sum(a(:,idx==j));
                end
            end
            figure; imagesc(log10(M)); colorbar;
%}
            info = niftiinfo('data/hemiBranson7065atlasCal.nii.gz');
            aV = niftiread(info);
            roimax = max(aV(:));
            for i=1:roimax
                ridx = find(aV==i);
                aV(ridx) = idx(i);
            end
    
            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'k-means clustering based atlas';
            % output nii file
            niftiwrite(aV,atlas,info,'Compressed',true);
        end
        disp([atlas ' ROI count=' num2str(max(aV(:)))]);
    end

    % make ROI atlas based on k-means clustering of hemiem SC.
    % this requires hemiroiwhole_connectlist_cm.mat file. so need to run
    % makeStructConnectivity.m (whole flyem ROI) first.
    % this needs Statistics and Machine Learning Toolbox.
    for k=[20 30 50 100 200 300]
        atlas = ['data/' name 'Cmkm' num2str(k) 'atlasCal.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            aV = niftiread(atlasinfo);
        else
            load('data/hemiroiwhole_connectlist_cm.mat');
            cm = single(full(countMat));
            cm = log10(cm);
            cm(cm<0) = 0;
            % to avoid zero for k-means
            Xnorm = sqrt(sum(cm.^2, 2));
            cm(Xnorm==0,1) = eps(max(Xnorm)) * 2;
            clear Xnorm;

            idx = kmeans(cm,k,'Distance','cosine');
            clear cm;
            
            info = niftiinfo('data/hemiRoiWholeatlasCal.nii.gz');
            aV = niftiread(info);
            aidx = zeros(length(idx),1,'single');
            parfor i=1:length(aidx)
                aidx(i) = find(aV==i);
            end
            aV(aidx) = idx;
    
            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'k-means clustering based atlas';
            % output nii file
            niftiwrite(aV,atlas,info,'Compressed',true);
        end
        disp([atlas ' ROI count=' num2str(max(aV(:)))]);
    end
end
