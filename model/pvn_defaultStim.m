function [stim, prior, isOnStim, isOnPrior] = pvn_defaultStim(type, nTr, stim, prior)
% defines default stimulation parameters

isdual = numel(stim) > 1 && numel(prior) > 1;

switch type
    case 'stimOff5sec_wnPrior'  
        
        dur = 5e3;
        preDur = 0;
        priorScale = 0.3;
        assert(~isdual, 'Stim not defined for dual network');
        
        stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', [-Inf Inf], 'Scale', 0, 'Trial', 1:nTr);
        prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
     case 'stimOn5sec_wnPrior'  
        
        dur = 5e3;
        preDur = 0;
        stimScale = 1;
        priorScale = 0.3;
        assert(~isdual, 'Stim not defined for dual network');
        
        stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', [-Inf Inf], 'Scale', stimScale, 'Trial', 1:nTr);
        prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
     case 'dcPulse_wnPrior'
        
        [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef();
        
        assert(~isdual, 'Stim not defined for dual network');
        
        stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
    case 'dcPulse_wnPrior_wnStim'
        
        [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef();
        
        assert(~isdual, 'Stim not defined for dual network');
        
        stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        stimWnScale = 0.5*priorScale;
        stim.I = stim.I + stimWnScale.*rand(size(stim.I));
        prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
    case 'shortDCPulse_wnPrior'
        
        dur = 3e3;
        preDur = 4e3;
        stimScale = 1;
        priorScale = 0.3;
        
        stim.set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', [0 50], 'Scale', stimScale, 'Trial', 1:nTr);
        prior.set('WhiteNoise', 'Time', stim.t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
    case 'sharedDCPulse_sharedWnPrior'
        
        [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef();
        
        assert(isdual, 'Stim not defined for single network');
        
        stimIdx = strcmpi({stim.Name}, 'Stim_M');
        priorIdx = strcmpi({prior.Name}, 'Prior_M');
        
        stim(stimIdx).set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        
        prior(priorIdx).set('WhiteNoise', 'Time', stim(stimIdx).t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
        
    case 'dcPulseLeft_sharedWnPrior_sharedWnStim'
        
        % dc pulse to the left hemisphere, shared white-noise prior and
        % stim (stim-side noise is half scale)
        
        [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef();
        
        assert(isdual, 'Stim not defined for single network');
        
        stimIdxL = strcmpi({stim.Name}, 'Stim_L');
        stimIdxM = strcmpi({stim.Name}, 'Stim_M');
        priorIdxM = strcmpi({prior.Name}, 'Prior_M');
        
        stimWnScale = 0.5*priorScale;
        
        stim(stimIdxL).set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        stim(stimIdxM).set('WhiteNoise', 'Time', stim(stimIdxL).t, 'Scale', stimWnScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
        prior(priorIdxM).set('WhiteNoise', 'Time', stim(stimIdxL).t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
    case 'sharedDCPulse_sharedWnPrior_sharedWnStim'
        
        % dc pulse to both hemispheres, shared white-noise prior and
        % stim (stim-side noise is half scale)
        
        [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef();
        
        assert(isdual, 'Stim not defined for single network');
        
        stimIdxL = strcmpi({stim.Name}, 'Stim_L');
        stimIdxR = strcmpi({stim.Name}, 'Stim_R');
        stimIdxM = strcmpi({stim.Name}, 'Stim_M');
        priorIdxM = strcmpi({prior.Name}, 'Prior_M');
        
        stimWnScale = 0.5*priorScale;
        
        stim(stimIdxL).set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        stim(stimIdxR).set('Pulse', 'Dur', dur, 'Pre', preDur, 'On', isOnStim, 'Scale', stimScale, 'Trial', 1:nTr);
        stim(stimIdxM).set('WhiteNoise', 'Time', stim(stimIdxL).t, 'Scale', stimWnScale, 'Trial', 1:nTr, 'FixedSeed', false);
        
        prior(priorIdxM).set('WhiteNoise', 'Time', stim(stimIdxL).t, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);

    case '10sec_sqfNoiseStim_sqfNoisePrior'
        dur = 1e4;
        stimScale = 1.2;
        priorScale = 1.2;
        
        assert(~isdual, 'Stim not defined for dual network');
        
        [stim.I, prior.I] = deal(zeros(1,dur,nTr));
        for iTr = 1:nTr
            stim.I(1,:,iTr) = getSQFNoise(dur, stimScale);
            prior.I(1,:,iTr) = getSQFNoise(dur, priorScale);
        end
        stim.t = 1:dur;
        prior.t = 1:dur;
        
    case '10sec_wnStim_wnPrior'
        dur = 1e4;
        stimScale = 1.2;
        priorScale = 1.2;
        
        stim.set('WhiteNoise', 'Pre', 0, 'Dur', dur, 'Scale', stimScale, 'Trial', 1:nTr, 'FixedSeed', false);
        prior.set('WhiteNoise', 'Pre', 0, 'Dur', dur, 'Scale', priorScale, 'Trial', 1:nTr, 'FixedSeed', false);

    case '10sec_sqfNoiseStim'
        dur = 1e4;
        stimScale = 0.8;
        
        assert(~isdual, 'Stim not defined for dual network');
        
        stim.I = zeros(1,dur,nTr);
        for iTr = 1:nTr
            stim.I(1,:,iTr) = getSQFNoise(dur, stimScale);
        end
        stim.t = 1:dur;
        
    otherwise
        
        error('Not defined')
end

end
%%
function [dur, isOnStim, isOnPrior, stimScale, priorScale, preDur] = dcPulseDef()
dur = 2e3;
isOnStim = {[0 1e3]};
isOnPrior = {[-Inf Inf]};
stimScale = 1;
priorScale = 0.3;
preDur = 2e3;
end
%%
function out = getSQFNoise(dur, scale)
[wt, f] = cwt(rand(1,dur), 1e3);
wt = wt./abs(wt);
out = icwt(wt./sqrt(f));
out = scale*(0.2.*zscore(out) + 0.5);
end