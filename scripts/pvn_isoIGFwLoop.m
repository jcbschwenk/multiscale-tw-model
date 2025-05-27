% analysis to isolate the respective influence of the IG pacemaker and the
% recurrent FFw pathway on the dominant temporal frequency of the wave.
function [f, varDel] = pvn_isoIGFwLoop

%%

basefolder = pwd; % set manually
global tmpfolder

tmpfolder = [basefolder 'isoIGFwLoop/'];
mkdir(tmpfolder);


%%
dT = 3:40; % vary delay in integers

varDel(1) = varDelay_fixedIGPacemaker(dT);

varDel(2) = varDelay_nonRhythmicIG(dT);

%%

%% FIGURES:
col.IGX = [0.2 0.8 0.4];
col.SGX = [0.2 0.2 0.8];


f(1) = figure;
tiledlayout(1,3)

f(2) = figure;
tiledlayout(1,2)


for iCond = 1:2 % rhythmic / non-rhythmic
    
    figure(f(1))
    c = varDel(iCond);
    
    nexttile(iCond)
    
    % first dimension here and for pks is {'IGX' 'SGX'}
    
    freq = c.specOn(1).f;
    
    meanSpecOn = squeeze(mean(cat(4,c.specOn.spec),3));
    contourf(dT, freq, squeeze(meanSpecOn(2,:,:)), 64, 'EdgeColor', 'none');
    ylim([2 15])
    if iCond == 1
        ylabel('Frequency [Hz]')
    else
        yticks([])
    end
    xticks([10 20 30])
    if iCond == 2
    xlabel('Inter-Areal Delay [ms]')
    end
    colormap(brewermap(32, '*RdYlBu'))
    
    nexttile(3)
    hold on; box on
    pkOn = cat(3, c.specOn.pk);
    % -> dimord is [IGX/SGX x nTr x nDelay];
    
    if iCond == 1 % rhythmic IGX
        pvn_shplot(dT, squeeze(mean(pkOn(1,:,:),2)),quantile(squeeze(pkOn(1,:,:)), [0.025 0.975]), '-', 'LineWidth', 2, 'Color', col.IGX)
        pvn_shplot(dT, squeeze(mean(pkOn(2,:,:),2)),quantile(squeeze(pkOn(2,:,:)), [0.025 0.975]), '-', 'LineWidth', 2, 'Color', col.SGX + 0.2)
        text(dT(end), squeeze(mean(pkOn(1,:,end),2)), 'IGX', 'Color', col.IGX, 'FontWeight', 'Bold', 'HorizontalAlignment', 'right');
        text(dT(end), squeeze(mean(pkOn(2,:,end),2)), 'SGX', 'Color', col.SGX, 'FontWeight', 'Bold', 'HorizontalAlignment', 'right');
    else % non-rhythmic IGX
        pvn_shplot(dT, squeeze(mean(pkOn(1,:,:),2)),quantile(squeeze(pkOn(1,:,:)), [0.025 0.975]), '--', 'LineWidth', 2, 'Color', col.SGX - 0.2)
        text(dT(end), squeeze(mean(pkOn(1,:,end),2)), 'SGX [isolated]', 'Color', col.SGX - 0.2, 'FontWeight', 'Bold','HorizontalAlignment', 'right');
        
    end
    ylabel('Peak Freq. [Hz]')
    xticks([10 20 30])
    ylim([2 15]);
    xlim([min(dT) max(dT)])
    
    figure(f(2))
    nexttile(iCond)
    hold on; box on
    plot(dT, c.pfwOn, 'k-', 'LineWidth', 2)
    plot(dT, c.pfwOff, 'k--', 'LineWidth', 2)
    legend({'STIM ON' 'STIM OFF'}, 'Location', 'EastOutside')
    ylabel('P(FW)')
    xlabel('Inter-Areal Delay [ms]')
    
    
end

for i = 1:2
    figure(f(i))
nexttile(1)
title('SGX')
nexttile(2)
title('SGX [isolated]')
set(gcf, 'Color', 'w')
end

end
%%
function [g, param, m, stim, prior, nTr] = init()

[g, param, stim, prior] = pvn_mfmodel();
nTr = 50;
pvn_defaultStim('dcPulse_wnPrior_wnStim', nTr, stim, prior);
m = fxMapper(g);

end
%%
function getStack(m, param)
% Define pipeline for parameter searches:

% add frequency estimation to the pipeline:
tWin = {[-1e3 0] [0 1e3]};
tWinLabel = {'Off' 'On'};
for iWin = 1:2
    m.addToStack(['spec' tWinLabel{iWin}],...
        @(x, args)...
        pvn_getWTSpec(x, args{:}),...
        {m.Graph, 'TimeWindow', tWin{iWin}, 'Nodes', {'IGX', 'SGX'}, 'FreqBandPeak', [2 12]}, [1;1]);
end

% add second branch w/ eeg projection and wave estimation:
m.addToStack('eeg', @(x, args) pvn_eegProject(x, args{:}), {m.Graph, param}, [1; 1]); % eeg projection
[~, stackPosWav] = m.addToStack('wav', @(x, args) pvn_fitPlaneEEG(x, args{:}), {'Frequency', [7 13]}); % wave fit

for iWin = 1:2
    m.addToStack(['pfw' tWinLabel{iWin}],...
        @(x, args)...
        mean(x.pFitFW(x.t >= tWin{iWin}(1) & x.t < tWin{iWin}(2))),...
        {}, [stackPosWav; 1]);
end

end
%%
function out = varDelay_fixedIGPacemaker(dT)
global tmpfolder

%% 1) VARYING INTERAREAL DELAY

[g, param, m, stim, prior, nTr] = init();
getStack(m, param);

% get all synapses with inter-areal delay (excluding connections to the input-layers):
delayedSyn = find(cellfun(@(x) contains(x, {'FFw' 'FBk'}) && ~contains(x, 'IL'), param.SynLabels));

for iSyn = delayedSyn
    m.add(param.Syn(iSyn), 'T', dT, 1);
end

%
out = m.run('struct');
out.run = [];
save([tmpfolder 'varDelay_fixedIGPacemaker'], 'out', 'm', 'param');

end
%%
function out = varDelay_nonRhythmicIG(dT)
%% 2) VARYING INTERAREAL DELAY, WITHOUT IGX PACEMAKER
global tmpfolder

% get model and mapper:
[g, param, m, stim, prior, nTr] = init();
getStack(m, param);


% get all synapses with inter-areal delay (excluding connections to the input-layers):
delayedSyn = find(cellfun(@(x) contains(x, {'FFw' 'FBk'}) && ~contains(x, 'IL'), param.SynLabels));

for iSyn = delayedSyn
    m.add(param.Syn(iSyn), 'T', dT, 1);
end

% change IGX dynamics to non-rhythmic:
igxHandle = param.Nodes(find(strcmpi(param.NodeLabels, 'IGX'))); % all neurons in all IGX nodes are identical
igxHandle.SpkChild = fxSpkNeurIzh.empty; % make non-spiking
param.Syn(strcmpi(param.SynLabels, 'FBk_IGX_IGX')).W = 0; % disconnect IGX Feedback
param.Syn(strcmpi(param.SynLabels, 'FBk_ILPr_IGX')).W = 0; % disconnect prior

out = m.run('struct');
out.run = [];

save([tmpfolder 'varDelay_nonRhythmicIG'], 'out', 'm', 'param');

end