# multiscale-tw-model

MATLAB code used to generate the multiscale traveling wave model (Schwenk & Alamia, 2025)  
Requires any recent version of the fieldtrip toolbox to run


The basic version of the model is constructed by calling:

[g, param, stim, prior] = pvn_mfmodel();  
% returns a handle for the graph object (g), which can be used to alter the network  

% the parameter struct 'param' also contains handles to the nodes in the graph (identical between areas for nodes of the same type)  

nTr = 1; % simulate a single trial  
[stim, ~, isOnStim] = pvn_defaultStim('dcPulse_wnPrior', nTr, stim, prior);   
% set stimulus using default parameters (DC pulse & white noise top-down input to the infragranular node of the last area)  

mf = g.run(); % run the simulation  
% returns mf, a 1 x nTr array of structs, each containing activity for individual nodes in field 'r'  


scripts/allFigures.m contains the code to recreate the main figures from the paper  
