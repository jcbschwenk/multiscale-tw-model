function param = pvn_param()
% generate param struct for the model
param.N = 3;
param.Areas = 0:param.N+1;
param.fiCurve = @pvn_mfFICurve;
param.spkIConv = @(I) 4.*I;

% Architecture:
param.NodeLabels = {...
    'ILSt', 'ILPr',... % Input layers: stimulus and prior
    'ILStIN', 'ILPrIN',...
    'L4X', 'L4IN',... % Layer 4
    'SGX', 'SGIN', 'SGE',... % Supragranular
    'IGX', 'IGIN',... % Infragranular
    'Pul', 'PulIN'}; % Pulvinar

param.Nodes = fxNode.empty;


param.igxSpkNeur = fxSpkNeurIzh(...
    'a', 1/150,...
    'b', 0.2,...
    'c', -50,...
    'd', 2,...
    'e', [6 4]...
    );

igx = fxNode();
igx.SpkChild = repmat(param.igxSpkNeur, 1, 350);

% NODE SPECIFICATIONS:
for iType = param.NodeLabels
    
    ndType = char(iType);
    
    switch ndType
        
        case 'IGX'

            param.Nodes(end+1) = igx;
            
        otherwise % STANDARD NODE
            
            param.Nodes(end+1) = fxNode(...
                'lm', 270,...
                'b', 108,...
                'c', 0.154,...
                'Ib', 0.33,...
                'e', 0.025); 
    end
    
end
% -> each object in param.Nodes is a separate handle. This means that all
% nodes of the same type in the final graph will share parameters and
% changes in one will affect all others



% define which nodes exist in which area:
nodeInArea = cell(1, numel(param.Nodes));
% Default is all areas (not input layers), so we preset this for all node-types:
for iType = 1:numel(param.NodeLabels)
   nodeInArea{iType} = 1:param.N;
end

nodeInArea(strcmpi(param.NodeLabels, 'SGE')) = {1:param.N-1};
nodeInArea(strcmpi(param.NodeLabels, 'SGIN')) = {1:param.N-1};

nodeInArea(strcmpi(param.NodeLabels, 'Pul')) = {1:param.N};
nodeInArea(strcmpi(param.NodeLabels, 'PulIN')) = {1:param.N};

nodeInArea(strcmpi(param.NodeLabels, 'ILSt')) = {0};
nodeInArea(strcmpi(param.NodeLabels, 'ILStIN')) = {0};
nodeInArea(strcmpi(param.NodeLabels, 'ILPr')) = {param.N+1};
nodeInArea(strcmpi(param.NodeLabels, 'ILPrIN')) = {param.N+1};

param.nodeInArea = nodeInArea;




%% EEG forward model:

eeg.midline = {'Iz' 'Oz' 'POz' 'Pz' 'CPz' 'Cz' 'FCz'};

eeg.source_labels = {'Occ_L' 'Par_L' 'Front_L'                      'Occ_R' 'Par_R' 'Front_R'};
eeg.source_pos =    {[-8 -76 10] [-24 -61 58] [-5 48 30]            [8 -76 10] [24 -61 58] [5 48 30]};

eeg.num_noise_src = 5;
eeg.snr_range = [0.4 1.6];

% get forward model for EEG projection based on a standard headshape, and
% generate a grid of possible noise-source positions.
eeg = pvn_eegGetModel(eeg);  

eeg.sr = 100;

lowCutOff = 20;
[bpFilt.b,bpFilt.a]=butter(3,lowCutOff/(1e3/2),'low'); % applied before resampling, so sr = 1e3
eeg.bpFilt = bpFilt;

eeg.isfw = @(x) abs(x+pi/2) < 0.5; % determine fw state from fitted wavDir
eeg.isbw = @(x) abs(x-pi/2) < 0.5; % 

param.eeg = eeg;


end

