function f = pvn_freqSpectra

%%

basefolder = pwd; % set manually
global tmpfolder

tmpfolder = [basefolder 'freqSpectra/'];
mkdir(tmpfolder);


%%
numRep = 10;
numSteps = 20;
estArgs = {...
    'MaxCurr', 1.5,...
    'PriorScale', 0.3,...
    'NumSteps', numSteps,...
    'NumTrialsPerStep', numRep,...
    'Figure', false,...
    'SaveOutput', false};

[g, param] = pvn_mfmodel();
[dynStateThresh, ~, mf] = pvn_dynStateThresh(g, param, estArgs{:});

%%
tWin = {[-1e3 0] [0 1e3]};

for iWin = 1:2
    spec(iWin) = pvn_getWTSpec(mf, g, 'TimeWindow', tWin{iWin}, 'Nodes', {'IGX', 'SGX'}, 'FreqBandPeak', [2 12]);
end

%%
curr = dynStateThresh.curr;

for iNode = 1:2 % IGX/SGX
    for iWin = 1:2 % OFF/ON
        for iStep = 1:numel(curr)
            spec_avg(:,iStep,1) = squeeze(mean(spec(iWin).spec(iNode,:,numRep*(iStep-1)+1:numRep*iStep),3));
            spec_avg(:,iStep,2:3) = squeeze(quantile(spec(iWin).spec(iNode,:,numRep*(iStep-1)+1:numRep*iStep),[0.025 0.975],3));
            pk_avg(1,iStep) = squeeze(mean(spec(iWin).pk(iNode,numRep*(iStep-1)+1:numRep*iStep),2));
            pk_avg(2:3,iStep) = squeeze(quantile(spec(iWin).pk(iNode,numRep*(iStep-1)+1:numRep*iStep),[0.025 0.975]));
        end
        if iWin == 2
            pk_slope(iNode, 1) = mean(pk_avg(1,(curr-dynStateThresh.thresh) > 0 & (curr-dynStateThresh.thresh) < 0.1));
            pk_slope(iNode, 2) = mean(pk_avg(1,(max(curr)-curr) < 0.1));
        end
        
        figure(1)
        subplot(2,2,(iNode-1)*2 + iWin)
        contourf(curr, spec(iWin).f, spec_avg(:,:,1), 'EdgeColor', 'none');
        ylim([1 20])
        
        figure(2)
        subplot(1,3,iNode)
        qnt = squeeze(mean(spec_avg(:,curr > dynStateThresh.thresh,[2 3]),2));
        qnt(isnan(qnt)) = mean(qnt(~isnan(qnt)));
        ln(iNode,iWin) = pvn_shplot(spec(iWin).f, mean(spec_avg(:,curr > dynStateThresh.thresh,1),2), qnt);
        hold on;
        
        if iWin == 2 % STIM ON
            subplot(1,3,3)
            hold on
            pvn_shplot(curr, pk_avg(1,:), pk_avg(2:3,:))
        end
    end
    
end

figure(2)
nodes = {'IGX' 'SGX'};
for i = 1:2
    subplot(1,3,i)
    xline(9, 'k-')
    xlim([1 25])
    title(nodes{i})
    legend(ln(iNode,:), {'OFF' 'ON'})
    xlabel('Frequency [Hz]')
ylabel('Mean Wavelet Power [a.u.]')

end

subplot(1,3,3)
xline(dynStateThresh.thresh, 'k--')
ylim([3 12])
xlabel('Input Current [nA]')
ylabel('Peak Freq. [Hz]')

save([tmpfolder 'freqSpecs_workspace.mat'])

f = figure(2);
end
