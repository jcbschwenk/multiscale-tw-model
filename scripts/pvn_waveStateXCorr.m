function f1 = pvn_waveStateXCorr()

basefolder = pwd; % set manually
tmpfolder = [basefolder 'dynXCorrs/']; mkdir(tmpfolder)

[g, param, stim, prior] = pvn_mfmodel();
nTr = 50;

dur = 1e4;
stimScale = 1.6;
priorScale = 0.3;
stim_tau = 1e2;

for iTr = 1:nTr
    x = interp1(stim_tau:stim_tau:dur, round(rand(1,dur/stim_tau)), 1:dur, 'previous');
    x(isnan(x)) = round(rand);
    x = x.*stimScale;
    stt(1,:,iTr) = x;
    prt(1,:,iTr) = rand(1,dur).*priorScale;
end
stim.I = stt;
stim.t = 1:dur;
stim.Name = 'Stim';

prior.I = prt;
prior.t = 1:dur;
prior.Name = 'Prior';
%%

mf = g.run();
eeg = pvn_eegProject(mf, g, param);
wav = pvn_fitPlaneEEG(eeg, 'Frequency', [7 13]);

cx = pvn_getAreaAvg(mf, g, param, 'Cx');



% extract curr amplitudes and wave state and cross correlate:
for iTr = 1:nTr
    st(:,iTr) = resample(squeeze(stim.I(1,:,iTr)), param.eeg.sr, 1e3);
    pr(:,iTr) = resample(squeeze(prior.I(1,:,iTr)), param.eeg.sr, 1e3);
end
st = stimScale.*(abs(st) > 0.5.*stimScale); % eliminate resampling artifacts

%%
nShuffle = 10;
maxlag = 1.5; % in sec
for iTr = 1:nTr
    for wdir = {'fw' 'bw'}
        wavState = wav.(char(wdir))(:,iTr) & wav.sig(:,iTr);
        [xc_stim.(char(wdir))(:,iTr), lags] = ptBiSerialXCorr(wavState,st(:,iTr), maxlag*param.eeg.sr);
        
        for iShuffle = 1:nShuffle
            trIdx = randsample(setdiff(1:nTr, iTr), 1);
            xc_shuffle.(char(wdir))(:,iTr,iShuffle) = ptBiSerialXCorr(wavState,st(:,trIdx), maxlag*param.eeg.sr);
        end
    end
end

lags = 1e3.*(lags./param.eeg.sr);

%%
f1 = figure;
tiledlayout(3,1)
wdirs = {'fw' 'bw'};

nexttile([1 1])
plot(st(:,3), 'Color', 0.15.*[1 0.8 0.8], 'LineWidth', 2)
yticks([])
axis off
title('Stimulus')

nexttile([2 1])
hold on
leg = [];
for wdir = 1:2
    shmean = squeeze(mean(xc_shuffle.(wdirs{wdir}),2));
    shErr = quantile(shmean, [0.05 0.95], 2);
    pst.(wdirs{wdir}) = plot(lags, mean(xc_stim.(wdirs{wdir}),2), 'Color', pvn_figCols(wdirs{wdir}), 'LineWidth', 2);
    pvn_shplot(lags,mean(shmean,2),shErr,'Color', pvn_figCols(wdirs{wdir}), 'LineWidth', 1.5);
end

legend([pst.fw pst.bw], {'FW STIM' 'BW STIM' 'Stim Autocorr.'})
xlabel('Lag [ms]')
xlim([-500 1500])
xticks([-300 0 300 600 900 1200])
ylabel('Corr.)')
xline(0, 'k--')
set(gca, 'FontName', 'Arial')
set(gcf, 'Color', 'w')
title('Cross-correlation Input Current Wave State')

%%
save([tmpfolder 'xcorr_workspace.mat'])
%%

end
%%
function [xc, lags] = ptBiSerialXCorr(x,y, maxlag)
% computes the cross-correlation of x and y using the point-biserial
% correlation coefficient. x is categorical, y is continuous. Not optimized
% for speed (subfunction is called for every value)

lags = -maxlag:maxlag;
xc = nan(1,numel(lags));
for ilag = 1:numel(lags)
    t = lags(ilag);
    if t < 0
        xc(ilag) = pvn_ptBiSerialCorr(x(-t+1:end), y(1:end+t));
    elseif t >= 0
        xc(ilag) = pvn_ptBiSerialCorr(x(1:end-t), y(t+1:end));
    end
end

xc = flip(xc);

end
