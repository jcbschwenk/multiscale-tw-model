function [out, figHandle, mf] = pvn_dynStateThresh(g, param, varargin)
% estimates the dynamic state stim-curr threshold at which the network
% switches from BW to FW. Output contains mean state probabilities +
% parameters for a sigmoid fit of the form: y = a/(1+exp(-b*(x-c))) + d,
% with bounds [0 1] for a & d and b > 0.
%
% Uses a 2sec off- 3sec on stimulus pattern scaled in NumSteps steps
% between 0 and MaxCurr, running NumTrialsPerStep for each curr intensity

isdual = any(contains(param.NodeLabels, '_R'));

if ~isdual
    defStim = 'Stim';
else
    defStim = 'Stim_M';
end

p = inputParser();
p.addParameter('MaxCurr', 1);
p.addParameter('CurrLabel', defStim);
p.addParameter('PriorScale', 0.3);
p.addParameter('NumSteps', 10);
p.addParameter('NumTrialsPerStep', 10);
p.addParameter('Figure', false);
p.addParameter('Axes', []);
p.addParameter('SaveOutput', false);

p.parse(varargin{:});


if ~isdual
    stim = g.Curr(g.getCurrIdx(p.Results.CurrLabel));
    prior = g.Curr(g.getCurrIdx('Prior'));
    roi = {{[-0.2 0.2] [-Inf 0.25]}};
    adjustArg = false; % use A-P axis to classify wave direction
else
    stim = g.Curr(g.getCurrIdx(p.Results.CurrLabel));
    prior = g.Curr(g.getCurrIdx('Prior_M'));
    roi = {{[-0.2 0] [-Inf 0.25]}  {[0 0.2] [-Inf 0.25]}};
    adjustArg = true; % use median axis to classify wave direction
end

% get copies of stim and prior to reset to after:
stim0 = stim.copy;
prior0 = prior.copy;

% define stimulus as incr scaled versions of a DC pulse:

scale = linspace(0,p.Results.MaxCurr,p.Results.NumSteps);
idx = sort(repmat(1:p.Results.NumSteps, 1, p.Results.NumTrialsPerStep));

preDur = 2e3;
stimOn = 1:3e3;
t = (1:(preDur+max(stimOn))) - preDur;

nTr = p.Results.NumTrialsPerStep * p.Results.NumSteps;

stim.I = zeros(1,numel(t),nTr);
stim.I(1,:,:) = ismember(t, stimOn)'.*scale(idx);
prior.I = rand(1,numel(t),nTr).*p.Results.PriorScale;
[stim.t, prior.t] = deal(t);

mf = g.run();
eeg = pvn_eegProject(mf, g, param);

out = [];
out.roi = roi; % roi is last dim for all output fields

for iROI = 1:numel(roi)
    
    wav = pvn_fitPlaneEEG(eeg, 'Frequency', [7 13], 'ROI', roi{iROI}, 'AdjustDirAxis', adjustArg);
    
    tWin = wav.t >= stimOn(1)./1e3 & wav.t < stimOn(end)./1e3;
    
    nBoot = 1e2;
    
    [pfwOn, pbwOn] = deal(nan(p.Results.NumSteps, nBoot));
    for iStep = 1:p.Results.NumSteps
        for iBoot = 1:nBoot
            trIdx = find(idx == iStep);
            rsample = randsample(trIdx,numel(trIdx),true);
            pfwOn(iStep, iBoot) = mean(mean(wav.fw(tWin, rsample) & wav.sig(tWin, rsample)));
            pbwOn(iStep, iBoot) = mean(mean(wav.bw(tWin, rsample) & wav.sig(tWin, rsample)));
        end
    end
    
    % fit a sigmoid:
    ft = fittype( 'a/(1+exp(-b*(x-c)))+d', 'Independent', 'x', 'Dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares');
    opts.Lower = [0 0.01 0 0];
    opts.Upper = [1 Inf p.Results.MaxCurr 1];
    
    [ftRez, gof] = fit(scale', mean(pfwOn,2), ft, opts);
    thresh = ftRez.c;
    
    out.thresh(iROI) = thresh;
    out.slope(iROI) = ftRez.b;
    out.lower(iROI) = ftRez.d;
    out.upper(iROI) = ftRez.a + ftRez.d;
    out.rsq(iROI) = gof.rsquare;
    out.pfw(:,iROI) = mean(pfwOn,2);
    out.pfwCI(:,:,iROI) = quantile(pfwOn, [0.025 0.975],2);
    out.pbw(:,iROI) = mean(pbwOn,2);
    out.pbwCI(:,:,iROI) = quantile(pbwOn, [0.025 0.975],2);    
    
    if p.Results.Figure
        if ~isempty(p.Results.Axes)
            axes(p.Results.Axes)
        else
            figHandle(iROI) = figure;
        end
        
        pvn_shplot(scale,out.pfw(:,iROI),out.pfwCI(:,:,iROI))
        
        hold on
        plot(ftRez)
        
        xline(thresh, 'k--')
        xlabel('Input Current [nA]')
        ylabel('P(FW)')
        text(thresh, 0.3, 'Threshold', 'FontWeight', 'Bold', 'FontName', 'Arial')
        legend off
        set(gcf, 'Color', 'w')
        ylim([-0.1 1.1])
        yticks([0 0.5 1])
    else
        figHandle = [];
    end
    
    
end

if p.Results.SaveOutput
    basefolder = pwd; % set manually
    tmpfolder = [basefolder 'stateThreshEstimation/']; mkdir(tmpfolder)
    save([tmpfolder 'workspace_' datestr(now, 'MMDDhhmm')])
end

out.curr = scale;


% reset stim and prior to original state:
stim.I = stim0.I;
stim.t = stim0.t;

prior.I = prior0.I;
prior.t = prior0.t;

end