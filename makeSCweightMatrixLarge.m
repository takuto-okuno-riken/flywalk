% make SC neuron weight matrix for huge matrix by hemibrain FlyEM structure data.

function [weightMat] = makeSCweightMatrixLarge(sycountMat, type)
    roimax = size(sycountMat,1);
    Nout = {}; Nin = {};
    nfile = ['results/cache/' type '_L_Nin_Nout.mat'];
    load(nfile);

    delete(gcp('nocreate')); % shutdown pools
    parpool(32);

    % calculate weight matrix (full, neurons, others)
    weightMat = sparse(roimax,roimax);
    parfor i=1:roimax
        if isempty(Nout{i}), continue; end
        disp(['ROI ' num2str(i) ' / ' num2str(roimax)]);
        outnids = Nout{i};
        S = sparse(1,roimax);
        for j=1:roimax
            if isempty(Nin{j}), continue; end
            innids = Nin{j};
            % find input neuron rate from ROI(i)
            logi = ismember(innids,outnids);
            S(j) = double(sum(logi)) / length(innids); % in-weight (from i to j)
        end
        weightMat(i,:) = S;
    end
    delete(gcp('nocreate')); % shutdown pools
    weightMat = weightMat .* sycountMat;
end
