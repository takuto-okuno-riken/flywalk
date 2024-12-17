% analyze Struct Connectivity.

function analyzeStructConnectivity
    % check SC post synapse cloud
    % roitype: hemiroi primary
    checkSCpostSynapse();

    % check neural transmitter type in each neuron
    % roitype: hemiroi primary
    checkNeuralTransmitter();

end

function checkSCpostSynapse()
    % make post-synapse could of FlyEM hemibrain
    rateThs = [50 60 70 80 90];
    synThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(rateThs)
        rateTh = rateThs(r);
        for j=1:length(synThs)
            synTh = synThs(j);

            niifile = ['data/hemibrain_hb' num2str(synTh) 'sr' num2str(rateTh) '_postsynFDACal.nii'];
            if exist([niifile '.gz'],'file')
                hbinfo = niftiinfo([niifile '.gz']);
                hbV = niftiread(hbinfo);
            else
                % read FlyEM hemibrain synapse info & location in FDA
                load('data/hemibrain_v1_2_neurons.mat');
                clear Nconn; clear Ncrop; clear Nsize; 
        
                load('data/hemibrain_v1_2_synapses.mat');
                clear Sloc; clear StoS;
                srate = (Srate >= (rateTh / 100)); % use only accurate synapse more than 'rate'
                straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
                spost = (Sdir == 2); % Find post synapse
        
                load('data/synapseloc_fdacal.mat');
                SpostlocFc = SlocFc(srate & straced & spost,:);
                clear SlocFc;
        
                hbinfo = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
                hbV = niftiread(hbinfo); % mask should have same transform with 4D nifti data
                hbV(:) = 0;
                sz = size(hbV);
        
                for i=1:size(SpostlocFc,1)
                    t = ceil(SpostlocFc(i,:));
                    if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                        hbV(t(1),t(2),t(3)) = hbV(t(1),t(2),t(3)) + 1;
                    else
                        disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
                    end
                end
                clear SpostlocFc;
        
                niftiwrite(hbV,niifile,hbinfo,'Compressed',true);
            end
        end
    end

    % make post-synapse could of FlyWire (hemibrain)
    rateThs = [50 70 100 130 140 150];
    synThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(rateThs)
        rateTh = rateThs(r);
        for j=1:length(synThs)
            synTh = synThs(j);

            niifile = ['data/hemibrain_fw' num2str(synTh) 'sr' num2str(rateTh) '_postsynFDACal.nii'];
            if exist([niifile '.gz'],'file')
                fwinfo = niftiinfo([niifile '.gz']);
                fwV = niftiread(fwinfo);
            else
                % read synapse info
                load('data/flywire783_synapse.mat');
                score = (cleftScore >= rateTh);
                valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.
            
                % read synapse location in FDA
                load('data/flywire783i_sypostloc_fdacal.mat');
                SpostlocFc = SpostlocFc(valid & score,:);
        
                fwinfo = niftiinfo('template/thresholded_FDACal_mask.nii.gz');
                fwV = niftiread(fwinfo); % mask should have same transform with 4D nifti data
                fwV(:) = 0;
                sz = size(fwV);
        
                mV = imgaussfilt3(hbV, [10 10 10] / sqrt(8*log(2))); mV(mV>0)=1; % make mask
        
                for i=1:size(SpostlocFc,1)
                    t = ceil(SpostlocFc(i,:));
                    if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                        if mV(t(1),t(2),t(3)) > 0
                            fwV(t(1),t(2),t(3)) = fwV(t(1),t(2),t(3)) + 1;
                        end
                    else
                        disp(['out of bounds ' num2str(i) ') ' num2str(t)]);
                    end
                end
                clear SpostlocFc;
        
                niftiwrite(fwV,niifile,hbinfo,'Compressed',true);
            end
        end
    end

    % load SC & atlas
    roitypes = {'hemiroi'};

    for i = 1:length(roitypes)
        roitype = roitypes{i};

        fname = ['data/' roitype '_postsyncount.mat'];
        if exist(fname,'file')
            load(fname);
        else
            roiIdxs = {};
            switch(roitype)
            case 'hemiroi'
                primaryIds = [107	16	59	68	65	78	4	49	106	87	100	27	43	5	57	89	101	97	50	58	113	10	32	66	30	67	19	76	31	82	93	54	52	8	7	42	1	63	95	112	98	33	18	103	15	20	111	34	51	62	47	24	38	22	75	41	2	45	80	102	56	28	91];
                listing = dir(['atlas/' roitype '/*.nii.gz']);
                for j=1:length(listing)
                    V = niftiread(['atlas/' roitype '/roi' num2str(j) '.nii.gz']); % ROI mask should have same transform with 4D nifti data
                    roiIdxs{j} = find(V>0);
                end
            end

            hbrateThs = [50 60 70 80 90];
            hbsynThs = [0];
            fwrateThs = [50 70 100 130 140 150];
            fwsynThs = [0];
            hbS = nan(length(roiIdxs),length(hbrateThs),length(hbsynThs),'single');
            fwS = nan(length(roiIdxs),length(fwrateThs),length(fwsynThs),'single');

            for k=1:length(roiIdxs)
                for r=1:length(hbrateThs)
                    rateTh = hbrateThs(r);
                    for c=1:length(hbsynThs)
                        synTh = hbsynThs(c);
                        hbV = niftiread(['data/hemibrain_hb' num2str(synTh) 'sr' num2str(rateTh) '_postsynFDACal.nii.gz']);
                        hbSc = hbV(roiIdxs{k});
                        hbS(k,r,c) = sum(hbSc);
                    end
                end
                for r=1:length(fwrateThs)
                    rateTh = fwrateThs(r);
                    for c=1:length(fwsynThs)
                        synTh = fwsynThs(c);
                        fwV = niftiread(['data/hemibrain_fw' num2str(synTh) 'sr' num2str(rateTh) '_postsynFDACal.nii.gz']);
                        fwSc = fwV(roiIdxs{k});
                        fwS(k,r,c) = sum(fwSc);
                    end
                end
            end
            save(fname,'hbrateThs','hbsynThs','fwrateThs','fwsynThs','hbS','fwS','primaryIds','-v7.3');
        end
    end
end

function checkNeuralTransmitter()
    % read neuron info (id, type)
    load('data/flywire783_neuron.mat'); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % read synapse info
    load('data/flywire783_synapse.mat');

    % make transmitter type of FlyWire (hemibrain)
    rateThs = [50 70 100 130 140 150];
    synThs = [0];
    roitypes = {'hemiroi'};

    for n = 1:length(roitypes)
        for r=1:length(rateThs)
            rth = rateThs(r);
            for j=1:length(synThs)
                sth = synThs(j);
                idstr = [roitypes{n} '_fw' num2str(sth) 'sr' num2str(rth)];
                fname = ['data/' idstr '_transmitter.mat'];
                if ~exist(fname,'file')
                    nfile = ['results/cache-' idstr '_Nin_Nout.mat'];
                    load(nfile);
    
                    roiNum = length(Nout);
                    outNTypes = zeros(roiNum,7,'single');
                    inNTypes = zeros(roiNum,7,'single');
                    inSyTypes = zeros(roiNum,7,'single');
                    for i=1:roiNum
                        nidx = Nout{i}{2};
                        types = Ntype(nidx);
                        numtype = groupcounts(types); % number of neurons in each transmitter
                        tidx = unique(types);
                        outNTypes(i,tidx+1) = numtype;

                        nidx = Nin{i}{2};
                        types = Ntype(nidx);
                        numtype = groupcounts(types); % number of neurons in each transmitter
                        tidx = unique(types);
                        inNTypes(i,tidx+1) = numtype;

                        sidx = Sin{i}{2}; % pre-synapse in this ROI
                        presnidx = preNidx(sidx); % get neuron index of pre-synapse
                        [nidx2,ia,ic] = unique(presnidx);
                        % check neuron index nidx should be same as nidx2
                        if sum(nidx-nidx2) ~= 0
                            disp([idstr ' nidx ~= nidx2 i=' num2str(i)])
                        end
                        prestypes = types(ic); % get pre-synapse transmitter
                        numtype = groupcounts(prestypes); % number of pre-synapses in each transmitter
                        tidx = unique(prestypes);
                        inSyTypes(i,tidx+1) = numtype;
                    end
                    save(fname,'inNTypes','outNTypes','inSyTypes','-v7.3');
                end
            end
        end
    end
end
