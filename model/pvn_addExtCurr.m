function [g, param, stim, prior] = pvn_addExtCurr(g, param, varargin)
% add external current (stim/prior) to the network
% inputs are g: graph object; param: parameter struct; 
% optional arguments are: 'Stim': Stimulus object, 'Prior': Prior object


isdual = any(contains(param.NodeLabels, '_R'));

p = inputParser();
p.addOptional('Stim', fxCurr('Name', 'Stim'));
p.addOptional('Prior', fxCurr('Name', 'Prior'));
p.parse(varargin{:});


stim = p.Results.Stim;
prior = p.Results.Prior;

% dual networks have hemisphere-specific currents (L/R) and shared (M) for
% both stim and prior

if isdual && numel(stim) == 1
    stim = [stim stim.copy stim.copy];
    stim(1).Name = [stim(1).Name '_L'];
    stim(2).Name = [stim(2).Name '_R'];
    stim(3).Name = [stim(3).Name '_M'];
end
if isdual && numel(prior) == 1
    prior = [prior prior.copy prior.copy];
    prior(1).Name = [prior(1).Name '_L'];
    prior(2).Name = [prior(2).Name '_R'];
    prior(3).Name = [prior(3).Name '_M'];
end

assert(numel(stim) == numel(prior));

if ~isdual
    g.addCurr(stim, 'Stim')
    g.addCurr(prior, 'Prior')
else
    g.addCurr(stim(1), 'Stim_L')
    g.addCurr(stim(2), 'Stim_R')
    g.addCurr(stim(3), 'Stim_M')
    g.addCurr(prior(1), 'Prior_L')
    g.addCurr(prior(2), 'Prior_R')
    g.addCurr(prior(3), 'Prior_M')
end

param.Curr = [stim prior];
param.CurrLabels = {param.Curr.Name};

end