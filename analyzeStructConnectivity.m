% analyze Struct Connectivity.

function analyzeStructConnectivity
    % check SC post synapse cloud FlyEM, FlyWire
    % roitype: hemiroi primary
    checkSCpostSynapse();

    % check neural transmitter type in each neuron (FlyWire)
    % roitype: hemiroi primary
    checkNeuralTransmitterFw();

    % check mushroom body huge GABA neuron (APL-R) FlyEM, FlyWire
    checkAPLneuron();
    checkAPLneuronFw();
end

function checkAPLneuron()
    confTh = 0.8; % FlyEM hemibrain synapse confidence threshold

    % FlyEM read neuron info (id, connection number, size)
    load('data/hemibrain_v1_2_neurons.mat');
    clear Nconn; clear Ncrop; clear Nsize; 

    % FlyEM read synapse info
    Sdir = []; StoN = []; Srate = []; StoS = [];
    load('data/hemibrain_v1_2_synapses.mat');
    clear Sloc;
    Sid = uint32(1:length(StoN));
    srate = (Srate >= confTh); % use only accurate synapse more than confidence threshold
    straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
    s1rate = ismember(StoS(:,1),Sid(srate));
    s2rate = ismember(StoS(:,2),Sid(srate));
    ssrate = (s1rate & s2rate);
    s1traced = ismember(StoS(:,1),Sid(straced));
    s2traced = ismember(StoS(:,2),Sid(straced));
    sstraced = (s1traced & s2traced);
    clear straced; clear s1rate; clear s2rate; clear s1traced; clear s2traced;

    Sdir(Srate < confTh) = 0;  % use only accurate synapse more than confidence threshold
    clear Srate;

    % FlyEM read synapse location in FDA
    load('data/hemibrain_v1_2_synapseloc_fdacal.mat');

    % count pre (to post) and post synapse count
    APLnid = 425790257;
    logi = ismember(StoN,APLnid); % find synapses which belong to APL neuron
    presids = Sid(logi & Sdir==1);
    postsids = Sid(logi & Sdir==2);
    logi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
    cpostsids = StoS(logi & ssrate & sstraced,2);
    [cpostsids, ia] = unique(cpostsids);

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0; V = Vt;
    sz = size(Vt);

    % get connected post-synapse counts from APL pre-synapses (output)
    fname = 'results/nifti/hemiAplOutputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        conSlocFc = SlocFc(cpostsids,:); % get 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    else
        V = niftiread([fname '.gz']);
    end
    fname = 'results/nifti/hemiAplOutputSynapsesMask.nii';
    if ~exist([fname '.gz'],'file')
        V(V>0) = 1;
        niftiwrite(V,fname,info,'Compressed',true);
    end

    % get post-synapse count of APL neuron (input)
    fname = 'results/nifti/hemiAplInputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        conSlocFc = SlocFc(postsids,:); V = Vt; % get 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    end

    % get connected post-synapse counts from all traced neurons
    fname = 'results/nifti/hemiTracedNeuronOutputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        logi = ismember(StoN,Nid(Nstatus==1)); % find synapses which belong to traced neuron
        presids = Sid(logi & Sdir==1);
        logi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
        cpostsids = StoS(logi & ssrate & sstraced,2);
        [cpostsids, ia] = unique(cpostsids);
    
        conSlocFc = SlocFc(cpostsids,:); V = Vt; % get 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    end
    clear SlocFc; clear Nid; clear Nstatus; clear Sdir; clear Sid;
end

function checkAPLneuronFw()
    scoreTh = 140; % FlyWire synapse score threshold
    conf = getSCconfig('wire', 0, scoreTh);

    % FlyWire read neuron info
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % FlyWire read synapse info
    load(conf.synapseFile);
    score = (cleftScore >= scoreTh);
    Sidx = int32(1:length(Sid))';
    valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.

    % read synapse location in FDA
    load(conf.syprelocFdaFile);
    load(conf.sypostlocFdaFile);

    APLnid = int64(720575940613583001);
    APLidx = find(Nid==APLnid);
    logi = ismember(preNidx,APLidx); % find pre-synapses which belong to APL neurons
    presidx = Sidx(logi & valid & score);
    logi = ismember(postNidx,APLidx); % find post-synapses which belong to APL neurons
    postsidx = Sidx(logi & valid & score);

    info = niftiinfo('template/thresholded_FDACal.nii.gz');
    Vt = niftiread(info); Vt(:) = 0; V = Vt;
    sz = size(Vt);

    % get connected post-synapse counts from APL pre-synapses (output)
    fname = 'results/nifti/wireAplOutputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        conSlocFc = SpostlocFc(presidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    else
        V = niftiread([fname '.gz']);
    end
    fname = 'results/nifti/wireAplOutputSynapsesMask.nii';
    if ~exist([fname '.gz'],'file')
        V(V>0) = 1;
        niftiwrite(V,fname,info,'Compressed',true);
    end

    fname = 'results/nifti/wireAplInputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        conSlocFc = SpostlocFc(postsidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    end

    % get connected post-synapse counts from all traced neurons
    fname = 'results/nifti/wireTracedNeuronOutputSynapses.nii';
    if ~exist([fname '.gz'],'file')
        presidx = Sidx((preNidx>0) & valid & score);

        conSlocFc = SpostlocFc(presidx,:); V = Vt; % get (pre-post) 3D location in FDA Cal template.
        for j=1:size(conSlocFc,1)
            t = ceil(conSlocFc(j,:));
            if t(1)>0 && t(2)>0 && t(3)>0 && t(1)<sz(1) && t(2)<sz(2) && t(3)<sz(3) 
                V(t(1),t(2),t(3)) = V(t(1),t(2),t(3)) + 1;
            else
                disp(['out of bounds ) ' num2str(t)]);
            end
        end
        niftiwrite(V,fname,info,'Compressed',true);
    end
    clear SpostlocFc;
end

function checkSCpostSynapse()
    if ~exist('results/nifti','dir'), mkdir('results/nifti'); end

    % make post-synapse could of FlyEM hemibrain
    confThs = [0 10 20 30 40 50 60 70 80 90];
    hbsynThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(confThs)
        confTh = confThs(r);
        for j=1:length(hbsynThs)
            synTh = hbsynThs(j);

            niifile = ['results/nifti/hemibrain_hb' num2str(synTh) 'sr' num2str(confTh) '_postsynFDACal.nii'];
            if exist([niifile '.gz'],'file')
                hbinfo = niftiinfo([niifile '.gz']);
                hbV = niftiread(hbinfo);
            else
                % read FlyEM hemibrain synapse info & location in FDA
                load('data/hemibrain_v1_2_neurons.mat');
                clear Nconn; clear Ncrop; clear Nsize; 
        
                load('data/hemibrain_v1_2_synapses.mat');
                clear Sloc;
                Sid = uint32(1:length(StoN))';
                srate = (Srate >= (confTh / 100)); % use only accurate synapse more than confidence threshold
                straced = ismember(StoN,Nid(Nstatus==1)); % Find synapses belong to Traced neuron.
                s1rate = ismember(StoS(:,1),Sid(srate));
                s2rate = ismember(StoS(:,2),Sid(srate));
                ssrate = (s1rate & s2rate);
                s1traced = ismember(StoS(:,1),Sid(straced));
                s2traced = ismember(StoS(:,2),Sid(straced));
                sstraced = (s1traced & s2traced);

                presids = Sid(straced & Sdir==1);
                logi = ismember(StoS(:,1),presids); % get pre-synapse to connected post-synapse
                cpostsids = StoS(logi & ssrate & sstraced,2);
                [cpostsids, ia] = unique(cpostsids);

                load('data/hemibrain_v1_2_synapseloc_fdacal.mat');
                SpostlocFc = SlocFc(cpostsids,:); % get 3D location in FDA Cal template.
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
%            if confTh == 80
%                hbV(1,1,1) = 800; % for 80 (figure.3)
%                niftiwrite(hbV,niifile,hbinfo,'Compressed',true);
%            end
            disp(['hb' num2str(synTh) 'sr' num2str(confTh) ' maxV=' num2str(max(hbV(:)))]);
        end
    end

    % make post-synapse could of FlyWire (hemibrain)
    scoreThs = [50 70 100 130 140 150];
    fwsynThs = [0]; % 5 10 20 30 50 100];
    for r=1:length(scoreThs)
        scoreTh = scoreThs(r);
        for j=1:length(fwsynThs)
            synTh = fwsynThs(j);

            conf = getSCconfig('wire', synTh, scoreTh);
            niifile = ['results/nifti/hemibrain_fw' num2str(synTh) 'sr' num2str(scoreTh) '_postsynFDACal.nii'];
            if exist([niifile '.gz'],'file')
                fwinfo = niftiinfo([niifile '.gz']);
                fwV = niftiread(fwinfo);
            else
                % read synapse info
                load(conf.synapseFile);
                score = (cleftScore >= scoreTh);
                valid = (postNidx>0 & preNidx>0); % Find synapses belong to Traced neuron.
            
                % read synapse location in FDA
                load(conf.sypostlocFdaFile);
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
%            if scoreTh == 140
%                fwV(1,1,1) = 800; % for 140 (figure.3)
%                niftiwrite(fwV,niifile,hbinfo,'Compressed',true);
%            end
            disp(['fw' num2str(synTh) 'sr' num2str(scoreTh) ' maxV=' num2str(max(fwV(:)))]);
        end
    end

    % load SC & atlas
    roitypes = {'hemiroi', 'hemidistkm500'};

    for i = 1:length(roitypes)
        roitype = roitypes{i};

        fname = ['results/sc/' roitype '_postsyncount.mat'];
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
            otherwise
                V = niftiread(['atlas/' roitype 'atlasCal.nii.gz']); % ROI mask should have same transform with 4D nifti data
                roitype = lower(roitype);
                roimax = max(V(:));
                for i=1:roimax
                    roiIdxs{i} = find(V==i);
                end
                primaryIds = 1:roimax;
            end

            hbS = nan(length(roiIdxs),length(confThs),length(hbsynThs),'single');
            fwS = nan(length(roiIdxs),length(scoreThs),length(fwsynThs),'single');

            for k=1:length(roiIdxs)
                for r=1:length(confThs)
                    scoreTh = confThs(r);
                    for c=1:length(hbsynThs)
                        synTh = hbsynThs(c);
                        hbV = niftiread(['results/nifti/hemibrain_hb' num2str(synTh) 'sr' num2str(scoreTh) '_postsynFDACal.nii.gz']);
                        hbSc = hbV(roiIdxs{k});
                        hbS(k,r,c) = sum(hbSc);
                    end
                end
                for r=1:length(scoreThs)
                    scoreTh = scoreThs(r);
                    for c=1:length(fwsynThs)
                        synTh = fwsynThs(c);
                        fwV = niftiread(['results/nifti/hemibrain_fw' num2str(synTh) 'sr' num2str(scoreTh) '_postsynFDACal.nii.gz']);
                        fwSc = fwV(roiIdxs{k});
                        fwS(k,r,c) = sum(fwSc);
                    end
                end
            end
            save(fname,'confThs','hbsynThs','scoreThs','fwsynThs','hbS','fwS','primaryIds','-v7.3');
        end
    end
end

function checkNeuralTransmitterFw()
    conf = getSCconfig('wire', 0, 140);

    % read neuron info (id, type)
    load(conf.neuronFile); % type, da(1),ser(2),gaba(3),glut(4),ach(5),oct(6)

    % read synapse info
    load(conf.synapseFile);

    % make transmitter type of FlyWire (hemibrain)
    rateThs = [50 70 100 130 140 150];
    synThs = [0];
    roitypes = {'hemiroi', 'hemidistkm500'};

    for n = 1:length(roitypes)
        for r=1:length(rateThs)
            rth = rateThs(r);
            for j=1:length(synThs)
                sth = synThs(j);
                idstr = [roitypes{n} '_fw' num2str(sth) 'sr' num2str(rth)];
                fname = ['results/sc/' idstr '_transmitter.mat'];
                if ~exist(fname,'file')
                    nfile = ['results/cache/' idstr '_Nin_Nout.mat'];
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
