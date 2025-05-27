% runs the wave analysis on the human EEG data from Pang et al. (2020).
% returns figure of state probabilities and wave fit output

function [f1, wav] = pvn_expEEGWaveAnalysis(loadFromDisk)


basefolder = pwd; % set manually 
tmpfolder = [basefolder 'eegDataAnalysis/'];

conds = {'dynamic' 'static'};

if nargin && loadFromDisk
    
    load([tmpfolder 'wavFits'], 'wav');
    
else
    
    mkdir(tmpfolder);
    
    load('dataA3.mat', 'dataA3');
    
    nSubs = numel(dataA3);
    
    
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
    
    wav = cell(1,2);
    
    for iSub = 1:nSubs
        d = dataA3(iSub);
        fprintf('\nSubject #%d ', iSub);
        
        % create data structure produced by pvn_eegProject:
        eeg = [];
        eeg.t = d.timepoints; % in secs
        eeg.label = d.channelInfo;
        
        for iCond = 1:2
            
            fprintf('Condition %d/2 ... ', iCond);
            
            eeg.eeg = double(d.([conds{iCond} 'EEG']));
            % dimord is [nChan x nTm x nTr];
            
            wav{iCond}(iSub) = pvn_fitPlaneEEG(eeg, fitArgs{:});
        end
        
    end
    
    save([tmpfolder 'wavFits'], 'wav');
    
end

%%
t = wav{1}(1).t;

for iCond = 1:2
    fw(:,:,iCond) = cat(2,wav{iCond}.pFitFW);
    bw(:,:,iCond) = cat(2,wav{iCond}.pFitBW);
end

f1 = figure;
tiledlayout(1,2);
condcols = {[0.7 0.2 1] [1 0.4 0]};

for iCond = 2
    nexttile(1)
    hold on
    
    fwln = pvn_shplot(t, mean(fw(:,:,iCond),2), std(fw(:,:,iCond), [], 2), 'Color', pvn_figCols('fw'));
    bwln = pvn_shplot(t, mean(bw(:,:,iCond),2), std(bw(:,:,iCond), [], 2), 'Color', pvn_figCols('bw'));
    
    title(conds{iCond})
    xlim([-0.4 9.4])
    yline(0, 'k:')
    yline(1, 'k:')
    patch([0 5 5 0], [-0.1 -0.1 1.1 1.1], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    ylim([0.15 0.5])
    ylabel('P(State)')
    box on
    
    nexttile(2)
    dd = log(fw(:,:,iCond)./bw(:,:,iCond));
    ln(1) = pvn_shplot(t, mean(dd,2), std(dd, [], 2), 'Color', condcols{iCond});
    xlim([-0.4 9.4])
    ylim([-1.5 0.5])
    yline(0, 'k:')
    patch([0 5 5 0], [-2 -2 2 2], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    ylabel('log[P_{FW}/P_{BW}]')
    box on
    
end
nexttile(1)
legend([fwln bwln], {'FW' 'BW'}, 'Location', 'Southwest')
legend('boxoff')
xlabel('Time [secs]')

nexttile(2)
legend(ln, conds{2}, 'Location', 'Southwest')
legend('boxoff')

set(gcf, 'Color', 'w')

%%
end
%%
