% runs the simulation used to compare model output to real EEG data (Pang et
% al.); input arg is the stimulus amplitude to use 
% (threshold value is used in the paper for maximum dynamic range)
% variable basefolder needs to be set manually

function [f1, f2] = pvn_expEEGSimulate(stimScale)
basefolder = pwd; % set manually

global tmpfolder 
tmpfolder = [basefolder 'eegDataAnalysis/'];

dur = 11.5e3;
preDur = 1.5e3;
stimDur = 5e3;
priorScale = 0.2 * stimScale;

fitArgs = {...
    'Frequency', [7 13],...
    'ROI', {[-0.2 0.2] [-Inf 0.25]},...
    'WindowSize', 100,...
    'NumStepsSpatFreq', 30,...
    'NumStepsWaveDir', 60,...
    'RandShuffleIter', 10,...
    'RandShuffleTimeStep', 200,...
    'DirClassWindowSize', 0.5,...
    };

[g, param, stim, prior] = pvn_mfmodel();
nTr = 1e2;


% STATIC:
stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', [0 stimDur], 'Scale', stimScale, 'Trial', 1:nTr);
prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);

mf{1} = g.run();

eeg(1) = pvn_eegProject(mf{1}, g, param);
wav(1) = pvn_fitPlaneEEG(eeg(1), fitArgs{:});


%%
conds = {'Static'};

f1 = figure();

tiledlayout(1,3);
condcols = {[0.7 0.2 1] [1 0.4 0]};

for iCond = 1
    nexttile(iCond)
    bwln = pvn_shplot(wav(iCond).t, wav(iCond).pFitBW, wav(iCond).pFitBW_CI, 'Color', pvn_figCols('bw'));
    hold on
    fwln = pvn_shplot(wav(iCond).t, wav(iCond).pFitFW, wav(iCond).pFitFW_CI, 'Color', pvn_figCols('fw'));
    title(conds{iCond})
    xlim([-0.4 9.4])
    yline(0, 'k:')
    yline(1, 'k:')
    patch([0 5 5 0], [-0.1 -0.1 1.1 1.1], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
%     ylim([0.15 0.5])
    ylabel('P(State)')
    box on
    
    nexttile(3)
    hold on
    % bootstrap CI for the log ratio FW/BW:
    for iter = 1:1e2
        smpl = randsample(1:size(wav(iCond).fw,2), size(wav(iCond).fw,2), true);
        pFitFW = mean(wav(iCond).fw(:,smpl) & wav(iCond).sig(:,smpl),2);
        pFitBW = mean(wav(iCond).bw(:,smpl) & wav(iCond).sig(:,smpl),2);
        dd(:,iter) = log(pFitFW ./ pFitBW);
    end
    ci = quantile(dd, [0.025 0.975], 2);
    ci(ci == -Inf) = -10;
    ci(ci == Inf) = 10;
    ln(iCond) = pvn_shplot(wav(iCond).t, mean(dd,2), ci, 'Color', condcols{iCond});
    xlim([-0.4 9.4])
    ylim([-2.5 1.5])
    yline(0, 'k:')
    if iCond == 1
        patch([0 5 5 0], [-2 -2 2 2], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    ylabel('log[P_{FW}/P_{BW}]')
    box on
end

nexttile(1)
legend([fwln bwln], {'FW' 'BW'}, 'Location', 'Southwest');
legend('boxoff')
xlabel('Time [secs]')

nexttile(2)
legend(ln, conds, 'Location', 'Northeast');
legend('boxoff')

set(gcf, 'Color', 'w')



%% load wave fits for experimental data from disk:

expData = load([tmpfolder 'wavFits'], 'wav');
t = expData.wav{1}(1).t;

for iCond = 1:2
    fw(:,:,iCond) = cat(2,expData.wav{iCond}.pFitFW);
    bw(:,:,iCond) = cat(2,expData.wav{iCond}.pFitBW);
end
logp = log(fw(:,:,2)./bw(:,:,2));


% do the exponential fit on the time-course of the log-ratio:
isExp = true;
for iSubj = 1:size(logp,2)

    [ft(iSubj), pltFt(:,iSubj)] = pvn_expEEGFitWavResp(t, logp(:,iSubj), isExp);

end


% do the same for the simulated data:
% (fit on the original model param values)
logpSim = log(wav.pFitFW./wav.pFitBW);

isExp = false;
[ftSim, pltFtSim] = pvn_expEEGFitWavResp(wav.t, logpSim, isExp);


% small-range parameter search to map variability in fit parameters:
fitParamModelMap(g, param); % saves output to disk to save memory

%% load the results:
load([tmpfolder 'out_l4w'], 'out_l4w');
load([tmpfolder 'out_l4tau'], 'out_l4tau');
load([tmpfolder 'out_fbkw'], 'out_fbkw');

%% PARAM SEARCH PLOTS:

f2 = figure;

subplot(1,4,2)
hold on
cols = jet(numel(out_l4w.wav));
for i = 1:numel(out_l4w.wav) 
    plot(out_l4w.wav(i).t, smooth(log(out_l4w.wav(i).pFitFW./out_l4w.wav(i).pFitBW)',25), 'LineWidth', 2, 'Color', cols(i,:));
end
title('L4 Input Adaptation: Weight')


subplot(1,4,3)
hold on
cols = jet(numel(out_l4tau.wav));
for i = 1:numel(out_l4tau.wav) 
    plot(out_l4tau.wav(i).t, smooth(log(out_l4tau.wav(i).pFitFW./out_l4tau.wav(i).pFitBW)',25), 'LineWidth', 2, 'Color', cols(i,:));
end
title('L4 Input Adaptation: Tau')


subplot(1,4,4)
hold on
cols = jet(numel(out_fbkw.wav));
for i = 1:numel(out_fbkw.wav) 
    plot(out_fbkw.wav(i).t, smooth(log(out_fbkw.wav(i).pFitFW./out_fbkw.wav(i).pFitBW)',25), 'LineWidth', 2, 'Color', cols(i,:));
end
title('IG_{IB} Feedback: Weight')


subplot(1,4,1)
exSubs = [6 11];
for sub = exSubs
    hold on
    plot(t, smooth(logp(:,sub), 100), 'LineWidth', 2);
end
legend({'S1' 'S2'},'AutoUpdate','off')
title('Experimental Data')

% EXPERIMENT DATA: 2 Subs 
for sp = 1:4
    subplot(1,4,sp)
    xline(0, 'k:')
    yline(0, 'k-')
    patch([0 5 5 0], [-4 -4 4 4], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    ylim([-4 4])
%     xlim([-0.4 9.4])
    box on
end
end
%%
function fitParamModelMap(g, param)

global tmpfolder

m = fxMapper(g);

m.addToStack('eeg', @(x, args) pvn_eegProject(x, args{:}), {m.Graph, param}); % eeg projection
[~, stackPosWav] = m.addToStack('wav', @(x, args) pvn_fitPlaneEEG(x, args{:}), {'Frequency', [7 13]}); % wave fit
m.addToStack('fit', @(x, args) pvn_expEEGFitWavResp(x, args{:}), {false}); % eeg projection

L4AdaptSyn = param.Syn(strcmpi(param.SynLabels, 'Loc_L4IN_L4X'));
ILAdaptSyn = param.Syn(strcmpi(param.SynLabels, 'Loc_ILStIN_ILSt'));

% LAYER 4 input adaptation: weight
L4Adapt_w = -1.1:-0.05:-1.4;
m.add(L4AdaptSyn, 'W', L4Adapt_w);
m.add(ILAdaptSyn, 'W', L4Adapt_w, 1);
out_l4w = m.run('struct');
out_l4w = rmfield(out_l4w, 'run');
save([tmpfolder 'out_l4w'], 'out_l4w');
clear out_l4w


% time constant
m.clearGrid;
L4Adapt_tau = 100:5:140;
m.add(L4AdaptSyn, 'Tau', L4Adapt_tau);
m.add(ILAdaptSyn, 'Tau', L4Adapt_tau, 1);
out_l4tau = m.run('struct');
out_l4tau = rmfield(out_l4tau, 'run');
save([tmpfolder 'out_l4tau'], 'out_l4tau');
clear out_l4tau


% INFRAGRANULAR FEEDBACK weight
m.clearGrid;
IGIBFBkSyn = param.Syn(strcmpi(param.SynLabels, 'FBk_IGX_IGX'));
FBk_w = 0.8:0.05:1.2;
m.add(IGIBFBkSyn, 'W', FBk_w);

out_fbkw = m.run('struct');
out_fbkw = rmfield(out_fbkw, 'run');
save([tmpfolder 'out_fbkw'], 'out_fbkw');
clear out_fbkw

end