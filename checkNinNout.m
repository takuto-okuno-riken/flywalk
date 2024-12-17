function checkNinNout
    % check hb0 and old file
    t1 = load('results/cache-hemiroi_hb0sr80_Nin_Nout.mat');
    t2 = load('results/cache-hemiroi_Nin_Nout.mat');

    for i=1:length(t1.Nin)
        for p=1:3
            nids1 = t1.Nin{i}{p};
            nids2 = t2.Nin{i}{p};
            if length(nids1) ~= length(nids2)
                disp(['Nin bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(nids1-nids2)) > 0
                disp(['Nin bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    for i=1:length(t1.Nout)
        for p=1:3
            nids1 = t1.Nout{i}{p};
            nids2 = t2.Nout{i}{p};
            if length(nids1) ~= length(nids2)
                disp(['Nout bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(nids1-nids2)) > 0
                disp(['Nout bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    for i=1:length(t1.Sin)
        for p=1:3
            sids1 = t1.Sin{i}{p};
            sids2 = t2.Sin{i}{p};
            if length(sids1) ~= length(sids2)
                disp(['Sin bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(sids1-sids2)) > 0
                disp(['Sin bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    t1 = load('data/hemiroi_hb0sr80_connectlist.mat');
    t2 = load('data/hemiroi_connectlist.mat');
    for p=1:2
        mat1 = t1.ncountMat;
        mat2 = t2.ncountMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['ncountMat bad p=' num2str(p)]);
        end
    end

    for p=1:2
        mat1 = t1.sycountMat;
        mat2 = t2.sycountMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['sycountMat bad p=' num2str(p)]);
        end
    end

    for p=1:2
        mat1 = t1.nweightMat;
        mat2 = t2.nweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['nweightMat bad p=' num2str(p)]);
        end
    end

    for p=1:2
        mat1 = t1.outweightMat;
        mat2 = t2.outweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['outweightMat bad p=' num2str(p)]);
        end
    end

    for p=1:2
        mat1 = t1.syweightMat;
        mat2 = t2.syweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['syweightMat bad p=' num2str(p)]);
        end
    end

    % check fw0 and old file
    t1 = load('results/cache-hemiroi_fw0sr50_Nin_Nout.mat');
    t2 = load('results/cache-hemiroi_fw_Nin_Nout.mat');

    for i=1:length(t1.Nin)
        for p=2
            nids1 = t1.Nin{i}{p};
            nids2 = t2.Nin{i}{p};
            if length(nids1) ~= length(nids2)
                disp(['Nin bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(nids1-nids2)) > 0
                disp(['Nin bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    for i=1:length(t1.Nout)
        for p=2
            nids1 = t1.Nout{i}{p};
            nids2 = t2.Nout{i}{p};
            if length(nids1) ~= length(nids2)
                disp(['Nout bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(nids1-nids2)) > 0
                disp(['Nout bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    for i=1:length(t1.Sin)
        for p=2
            sids1 = t1.Sin{i}{p};
            sids2 = t2.Sin{i}{p};
            if length(sids1) ~= length(sids2)
                disp(['Sin bad length i=' num2str(i) ', p=' num2str(p)]);
            end
            if sum(abs(sids1-sids2)) > 0
                disp(['Sin bad nids i=' num2str(i) ', p=' num2str(p)]);
            end
        end
    end

    t1 = load('data/hemiroi_fw0sr50_connectlist.mat');
    t2 = load('data/hemiroi_fw_connectlist.mat');
    for p=2
        mat1 = t1.ncountMat;
        mat2 = t2.ncountMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['ncountMat bad p=' num2str(p)]);
        end
    end

    for p=2
        mat1 = t1.sycountMat;
        mat2 = t2.sycountMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['sycountMat bad p=' num2str(p)]);
        end
    end

    for p=2
        mat1 = t1.nweightMat;
        mat2 = t2.nweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['nweightMat bad p=' num2str(p)]);
        end
    end

    for p=2
        mat1 = t1.outweightMat;
        mat2 = t2.outweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['outweightMat bad p=' num2str(p)]);
        end
    end

    for p=2
        mat1 = t1.syweightMat;
        mat2 = t2.syweightMat;
        if sum(abs(mat1-mat2),'all') > 0
            disp(['syweightMat bad p=' num2str(p)]);
        end
    end
end
