% make voxel ROI atlas. Two types of atlas, such as cube ROI type, and
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
    for atlasSize = 6:-1:3
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

    % make neuropil specific voxel type atlas
    roiids = {[101],[57],[57,51]}; % FB, EB, EB-bL(L)
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
end
