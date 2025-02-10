% get structure config file (FlyEM, FlyWire)

function conf = getSCconfig(scname, synTh, scoreTh)
    conf.scname = scname;
    conf.synTh = synTh;        % connected synapse number at one neuron threshold
    conf.scoreTh = scoreTh;    % score threshold
    switch(scname)
    case 'hemi'
        conf.neuronFile = 'data/hemibrain1_2fw_neuron.mat';
        conf.synapseFile = 'data/hemibrain1_2fw_synapse.mat';
        conf.syprelocFile = 'data/hemibrain1_2fw_sypreloc.mat';
        conf.sypostlocFile = 'data/hemibrain1_2fw_sypostloc.mat';
        conf.syprelocFdaFile = 'data/hemibrain1_2fw_sypreloc_fdacal.mat';
        conf.sypostlocFdaFile = 'data/hemibrain1_2fw_sypostloc_fdacal.mat';
        conf.neuSepidxFile = 'data/hemibrain1_2fw_neuron_sepidx';
        conf.neuReciFile = 'data/hemibrain1_2fw_neuron_reci';
        conf.sySepidxFile = 'data/hemibrain1_2fw_synapse_sepidx';
        conf.syReciFile = 'data/hemibrain1_2fw_synapse_reci';
        conf.sypostCellFile = 'results/hemibrain1_2fw_sypostCell.mat';
        conf.swcPath = 'swc/hemibrain_v1_2';
        conf.swcSize = [1 1 1];   % swc unit is voxel
        conf.voxelSize = [8 8 8]; % nano meter
        conf.brainMeshFile = 'data/hemibrain_mesh.mat';
        conf.brainMeshView = [0 180];
    case 'wire'
        conf.neuronFile = 'data/flywire783_neuron.mat';
        conf.synapseFile = 'data/flywire783_synapse.mat';
        conf.syprelocFile = 'data/flywire783_sypreloc.mat';
        conf.sypostlocFile = 'data/flywire783_sypostloc.mat';
        conf.syprelocFdaFile = 'data/flywire783i_sypostloc_fdacal.mat';
        conf.sypostlocFdaFile = 'data/flywire783i_sypostloc_fdacal.mat';
        conf.neuSepidxFile = 'data/flywire783_neuron_sepidx';
        conf.neuReciFile = 'data/flywire783_neuron_reci';
        conf.sySepidxFile = 'data/flywire783_synapse_sepidx';
        conf.syReciFile = 'data/flywire783_synapse_reci';
        conf.sypostCellFile = 'results/flywire783i_sypostCell.mat';
        conf.swcPath = 'swc/flywire783';
        conf.swcSize = [4 4 40];   % swc unit is nano meter
        conf.voxelSize = [4 4 40]; % nano meter
        conf.brainMeshFile = 'data/flywire_mesh.mat';
        conf.brainMeshView = [0 -90];
    case 'wire630'
        conf.neuronFile = 'data/flywire630_neuron.mat';
        conf.synapseFile = 'data/flywire630_synapse.mat';
        conf.syprelocFile = '';
        conf.sypostlocFile = '';
        conf.syprelocFdaFile = '';
        conf.sypostlocFdaFile = '';
        conf.swcPath = 'swc/flywire630';
        conf.swcSize = [4 4 40];   % swc unit is nano meter
        conf.voxelSize = [4 4 40]; % nano meter
        conf.brainMeshFile = 'data/flywire_mesh.mat';
        conf.brainMeshView = [0 -90];
    end
    conf.voxelUnit = 'nm';
    conf.voxelSizeFda = [2.45, 2.28, 3.715]; % micro meter
    conf.voxelUnitFda = 'um';
end
