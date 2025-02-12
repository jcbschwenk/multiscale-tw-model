function eeg = pvn_eegGetModel(eeg)
% Get forward projection of mean-field output as cortical dipoles to EEG.

%% parse area input:
srcPos = cat(1, eeg.source_pos{:});

%% Get standard headmodel:
headmodel = ft_read_headmodel('standard_bem.mat');

%% Get electrode positions:

elec = ft_read_sens('standard_1020.elc');
elec = ft_convert_units(elec, 'm'); 

% realign to headmodel:
cfg = [];
cfg.method = 'project'; 
cfg.headshape = headmodel.bnd(1); 

elec = ft_electroderealign(cfg, elec);

%% Get source lead fields for the sources (for plotting only):
cfg             = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = srcPos;
cfg.headmodel   = headmodel; 
cfg.inwardshift = 0.005;
srcModel = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.sourcemodel = srcModel;
cfg.headmodel   = headmodel;
cfg.elec        = elec;
srcLeadField = ft_prepare_leadfield(cfg);

% Get a coarse grid of possible noise-source positions:
cfg             = [];
cfg.headmodel   = headmodel;
cfg.resolution = 5; 
cfg.inwardshift = 0.005;
nsModel = ft_prepare_sourcemodel(cfg);
nsPosGrid = nsModel.pos(nsModel.inside,:);

%%
eeg.ns_pos = nsPosGrid;
eeg.src_leadfield = srcLeadField;
eeg.headmodel = headmodel;
eeg.elec = elec;

end