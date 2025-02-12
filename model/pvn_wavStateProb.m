function [fw, bw] = pvn_wavStateProb(wav, param, tWin)
% returns probability of fw/bw state. if time window is passed, average is
% returned.

if nargin > 2
    tIdx = wav.t >= tWin(1) & wav.t < tWin(2);
    fw = mean(mean(param.eeg.isfw(wav.wavDir(tIdx,:)),2));
    bw = mean(mean(param.eeg.isbw(wav.wavDir(tIdx,:)),2));
else
    fw = mean(param.eeg.isfw(wav.wavDir),2);
    bw = mean(param.eeg.isbw(wav.wavDir),2);
end

end