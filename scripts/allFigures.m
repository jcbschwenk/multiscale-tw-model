% Main script to generate all results figures for the final paper 

clearvars

figFolder = pwd; % set 

rand_seed = 536234551;



%% Basic response behavior (cx-only)

[panelFig, phaseTopoFig] = pvn_singleCxStreamWaveDynamics();

savefig(panelFig, [figFolder 'Fig_SingleDynPanel'])
savefig(phaseTopoFig, [figFolder 'Fig_EEGTopoPhasePlotsFWBW'])

%% Dynamic state threshold estimation (cx-only):
rng(rand_seed)
estArgs = {...
    'MaxCurr', 1.8,...
    'PriorScale', 0.3,...
    'NumSteps', 15,...
    'NumTrialsPerStep', 20,...
    'Figure', true,...
    'SaveOutput', true};

[g, param] = pvn_mfmodel();
[dynStateThresh, dstFig] = pvn_dynStateThresh(g, param, estArgs{:});
savefig(dstFig, [figFolder 'Fig_DynStateThresh'])


%% Cross-Correlation Stim/Prior <-> Wave State
% uses 1/sqrt(f) noise for both
rng(rand_seed)
figXCorr = pvn_waveStateXCorr();
savefig(figXCorr, [figFolder 'Fig_XCorrWaveStateStimPrior'])

%% Comparison with real data:

% Experimental data:
rng(rand_seed)
loadFromDisk = true;
[figExpEEG, wavEEG] = pvn_expEEGWaveAnalysis(loadFromDisk); % Pang et al. 2020 (5sec dynamic and static stimulus)

% Simulation:
rng(rand_seed)
[figExpEEGSim, figExpEEGSubjParamVar] = pvn_expEEGSimulate(dynStateThresh.thresh);
% takes stimulus scale as input arg

savefig(figExpEEG, [figFolder 'Fig_expEEG'])
savefig(figExpEEGSim, [figFolder 'Fig_expEEGSimulation'])
savefig(figExpEEGSubjParamVar, [figFolder 'Fig_expEEGSimulation_SubjVar'])


%% Temp Frequency Spectra:
rng(rand_seed)
freqSpecFig = pvn_freqSpectra();
savefig(freqSpecFig, [figFolder 'Fig_freqSpec'])


%% Pathway isolation:
rng(rand_seed)
[isoFigs, varDel] = pvn_isoIGFwLoop();
savefig(isoFigs(1), [figFolder 'Fig_sgxIsoFreq'])
savefig(isoFigs(2), [figFolder '/suppl/Fig_sgxIsoFreq_waveStateProb'])


%% Pulvinar Feed-Up Weight Parameter Search:
rng(rand_seed)
[pulFigSingleStream, pulFigSingleStreamBounds] = pvn_fdupWeightSingleStream();
savefig(pulFigSingleStream, [figFolder 'Fig_PulFUpSingleStreamDST'])
savefig(pulFigSingleStreamBounds, [figFolder 'Fig_PulFUpSingleStreamBounds'])

%% Node-level comparison of evoked and spontaneous waves:
N = 3; % 5
[figEpochs, figPowDist] = pvn_fwWaveComparisonStimVSSpont(N);
savefig(figEpochs, [figFolder 'Fig_spontWavesEpochMaps_N' num2str(N)])
savefig(figPowDist, [figFolder 'Fig_spontWavesPowDist_N'  num2str(N)])


%% Dual-Stream Model:
rng(rand_seed)
dualStreamFig = pvn_dualStreamWaveDynamics();
savefig(dualStreamFig, [figFolder 'Fig_dualStreamLateralization'])


%% Pulvinar Feed-Up Weight Parameter Search for Dual Model:
rng(rand_seed)
[panelPulDualStream, boundsDSTPulDualStream] =  pvn_fdupWeightDualStreamParallel();
savefig(boundsDSTPulDualStream, [figFolder 'Fig_PulFUpDualStream'])
savefig(panelPulDualStream, [figFolder '/suppl/Fig_PulFUpDualStreamPanels'])

