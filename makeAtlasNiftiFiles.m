% make Atlas nifti files. This script outputs specific ROI atlas with mask image.
% need to run makeStructConnectivity.m first.

function makeAtlasNiftiFiles
    name = 'hemi';

    % read FDA mask
    minfo = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
    Vm = niftiread(minfo);
    Vm(Vm>0) = 1;
    Vm(Vm<1) = 0;

    % output atlas nifti file
    fname = 'hemiFlyem52atlasCal';
    info = niftiinfo(['atlas/' fname '.nii.gz']);
    V = niftiread(info);
    idx = find(V>0);
    aV = Vm;
    aV(aV>0) = max(V(:)) + 10;
    aV(idx) = V(idx);

    % output nii file with mask ROI. to see by ITK-SNAP
    niftiwrite(aV,['atlas/' fname 'wm.nii'],info,'Compressed',true);

    roitypes = {{'hemiBranson7065km',''}, {'hemiCmkm',''}, {'hemiCmkm','r1w1'}, {'hemiDistKm',''}, {'hemiRand',''}};
    for j=1:length(roitypes)
        for i=[50]
            fname = [roitypes{j}{1} num2str(i) roitypes{j}{2} 'atlasCal'];
            info = niftiinfo(['atlas/' fname '.nii.gz']);
            V = niftiread(info);
            idx = find(V>0);
            aV = Vm;
            aV(aV>0) = max(V(:)) + 10;
            aV(idx) = V(idx);
    
            % output nii file
            niftiwrite(aV,['atlas/' fname 'wm.nii'],minfo,'Compressed',true);
        end
    end
end
