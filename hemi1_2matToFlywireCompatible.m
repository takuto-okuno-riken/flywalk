% make flywire compatible data type from hemibrain 1.2 mat data (neuprint) 

function hemi1_2matToFlywireCompatible
    % transform neuron mat file
    fname = 'data/hemibrain1_2fw_neuron.mat';
    if exist(fname,'file')
        load(fname);
    else
        % FlyEM read neuron info (id, connection number, size)
        h = load('data/hemibrain_v1_2_neurons.mat');
        tlogi = (h.Nstatus==1); % traced neuron logical
        Nid = h.Nid(tlogi);
        Nscore = zeros(length(Nid),1,'single'); % [0 1] basically, unused.
        Ntype = zeros(length(Nid),1,'uint16'); % neural transmitter type. hemi doesn't have this info.
        Ncrop = h.Ncrop(tlogi); % croped or not. flywire doesn't have this info.
        Nsize = h.Nsize(tlogi); % neuron size. flywire doesn't have this info.
        save(fname,'Nid','Nscore','Ntype','Ncrop','Nsize','-v7.3');
    end

    % transform synapse mat file
    fname = 'data/hemibrain1_2fw_synapse.mat';
    if exist(fname,'file')
        load(fname);
    else
        % FlyEM read synapse info
        load('data/hemibrain_v1_2_synapses.mat');
        Sid = uint32(1:length(StoN))';
        straced = ismember(StoN,Nid); % Find synapses belong to Traced neuron.
        s1traced = ismember(StoS(:,1),Sid(straced));
        s2traced = ismember(StoS(:,2),Sid(straced));
        sstraced = (s1traced & s2traced);
        clear straced; clear s1traced; clear s2traced; clear Sid;

        preSid = StoS(sstraced,1); % flywire doesn't have preSid.
        Sid = StoS(sstraced,2);    % post-sid is used as Sid.
        preNid = StoN(preSid);
        postNid = StoN(Sid);
        [prenlogi,preNidx] = ismember(preNid,Nid);
        if sum(prenlogi) ~= length(preNid) % should be same
            disp('prenlogi error!');
        end
        if sum(preNid-Nid(preNidx)) ~= 0 % should be same
            disp('prenlogi error!');
        end
        [postnlogi,postNidx] = ismember(postNid,Nid);
        if sum(postnlogi) ~= length(preNid) % should be same
            disp('prenlogi error!');
        end
        if sum(postNid-Nid(postNidx)) ~= 0 % should be same
            disp('prenlogi error!');
        end
        preSrate = Srate(preSid);
        postSrate = Srate(Sid);
        cleftScore = uint16(min([preSrate,postSrate],[],2) * 100); % take smaller confidence of pre or post. Then, [0 1] -> [0 100].
        save(fname,'Sid','preSid','cleftScore','preNidx','postNidx','-v7.3');
        save('data/hemibrain1_2fw_synapse_nid.mat','preNid','postNid','-v7.3');

        % transform synapse location file
        Spreloc = Sloc(preSid,:);
        Spostloc = Sloc(Sid,:);
        save('data/hemibrain1_2fw_sypreloc.mat','Spreloc','-v7.3');
        save('data/hemibrain1_2fw_sypostloc.mat','Spostloc','-v7.3');
        % check pre-post distance
        D = sqrt(sum((double(Spreloc)-double(Spostloc)).^2,2));
        if any(D>100) % voxel (*8nm)
            disp('pre-post distance error!');
        end

        % FlyEM read synapse location in FDA
        load('data/hemibrain_v1_2_synapseloc_fdacal.mat');
        SprelocFc = SlocFc(preSid,:);
        SpostlocFc = SlocFc(Sid,:);
        save('data/hemibrain1_2fw_sypreloc_fdacal.mat','SprelocFc','-v7.3');
        save('data/hemibrain1_2fw_sypostloc_fdacal.mat','SpostlocFc','-v7.3');
        D = sqrt(sum((SprelocFc-SpostlocFc).^2,2));
        if any(D>1.0) % voxel
            disp('pre-post in FDA distance error!');
        end
    end

    % transform synapse separation index mat file (fw type to original)
    synTh = 0; scoreTh = 80; epsilon = 3000; minpts = 1;
    fname = ['data/hemibrain_v1_2_synapses_sepidx' num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat'];
    if exist(fname,'file')
        load(fname);
    else
        % FlyEM read synapse info
        load('data/hemibrain_v1_2_synapses.mat');

        % read separation index (fw type)
        load('data/hemibrain1_2fw_synapse.mat');
        load(['data/hemibrain1_2fw_synapse_sepidx' num2str(synTh) 'sr' num2str(scoreTh) '_' num2str(epsilon) 'mi' num2str(minpts) '.mat']);

        Spidx = ones(length(StoN),1,'int16') * -1;
        Spidx(Sid) = postSpidx;
        Spidx(preSid) = preSpidx;

        save(fname,'Spidx','-v7.3');
    end

    % transform synapse separation index mat file (fw type to original)
    synTh = 0; scoreTh = 80;
    fname = ['data/hemibrain_v1_2_synapses_reci' num2str(synTh) 'sr' num2str(scoreTh) '.mat'];
    if exist(fname,'file')
        load(fname);
    else
        % FlyEM read synapse info
        load('data/hemibrain_v1_2_synapses.mat');

        % read separation index (fw type)
        load('data/hemibrain1_2fw_synapse.mat');
        load(['data/hemibrain1_2fw_synapse_reci' num2str(synTh) 'sr' num2str(scoreTh) '.mat']);

        SrcCloseDist = nan(length(StoN),1,'single');
        SrcCloseDist(Sid) = SrcpostCloseDist;
        SrcCloseDist(preSid) = SrcpreCloseDist;

        SrcCloseSid = zeros(length(StoN),1,'int32'); % TODO:
%        SrcCloseSid(Sid) = preSid(SrcpostCloseSidx); % set close pre sid
%        SrcCloseSid(preSid) = Sid(SrcpreCloseSidx);  % set close post sid

        save(fname,'SrcCloseDist','SrcCloseSid','-v7.3');
    end
end
