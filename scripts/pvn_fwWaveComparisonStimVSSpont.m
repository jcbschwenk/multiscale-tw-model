function [fig1, fig2] = pvn_fwWaveComparisonStimVSSpont(N)

basefolder = pwd; % set manually
tmpfolder = [basefolder 'fwAnalysisSpont/']; mkdir(tmpfolder)
%%

conds = {'Stim ON' 'Stim OFF' 'Pul (Stim OFF)'};
comp = {'SGX' 'IGX'}; % laminar compartments

[g, param, stim, prior] = pvn_mfmodel(pvn_param(N));
bp = designfilt('bandpassfir',...
    'StopbandFrequency1', 6, 'PassbandFrequency1', 7,...
    'PassbandFrequency2', 13, 'StopbandFrequency2', 14,...
    'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 2e2);

fig1 = figure;
fig2 = figure;

for icond = 1:3
    
    nTr = 25;
    if icond == 1
        pvn_defaultStim('stimOn5sec_wnPrior', nTr, stim, prior);
    else
        pvn_defaultStim('stimOff5sec_wnPrior', nTr, stim, prior);
    end
    
    w = double(icond == 3)*2; % w = 2 if Pul
    
    fUpSyn = param.Syn(strcmpi(param.SynLabels, 'FUp_Pul_L4X'));
    fUpSyn.W = w;
    
    mf = g.run();
    cx = pvn_getAreaAvg(mf, g, param, 'Cx');
    
    
    
    
    cxs{icond} = [];
    apow{icond} = [];
    cut_edge = 40; % in samples at 200 Hz to cut from each edge
    for icomp = 1:2 % SGX/IGX
        
        for i = 1:param.N
            for iTr = 1:nTr
                sig = mf(iTr).r(g.getNodeIdx([comp{icomp} num2str(i)]), :);
                sig = resample(sig, 2e2, 1e3);
                sig_filt = filtfilt(bp, [flip(sig) sig flip(sig)]);
                cxs{icond}.filt.(comp{icomp})(i,:, iTr) = sig_filt;
                cxs{icond}.raw.(comp{icomp})(i,:, iTr) = sig;
                
                apow{icond}.(comp{icomp})(i,:,iTr) = abs(hilbert(sig_filt));
                
                
            end
        end
        N = size(cxs{icond}.filt.(comp{icomp}),2)/3;
        cxs{icond}.filt.(comp{icomp}) = cxs{icond}.filt.(comp{icomp})(:,N+1+cut_edge:2*N-cut_edge,:);
        cxs{icond}.raw.(comp{icomp}) = cxs{icond}.raw.(comp{icomp})(:,cut_edge+1:end-cut_edge,:);
        
        apow{icond}.(comp{icomp}) = apow{icond}.(comp{icomp})(:,N+1+cut_edge:2*N-cut_edge,:);        
    end
    
    % Cortex average:
    for i = 1:param.N
        for iTr = 1:nTr
            sig = cx(iTr).avg(i, :);
            sig = resample(sig, 2e2, 1e3);
            sig_filt = filtfilt(bp, [flip(sig) sig flip(sig)]);
            apow{icond}.cx(i,:,iTr) = abs(hilbert(sig_filt));
        end
    end
    N = size(apow{icond}.cx,2)/3;
    apow{icond}.cx = apow{icond}.cx(:,N+1+cut_edge:2*N-cut_edge,:);
    
    t = resample(mf(1).t, 2e2,1e3);
    t = t(cut_edge+1:end-cut_edge);
    
    % align to IGX2 peaks:
    wnsz = 60; % at 200 Hz
    
    for iTr = 1:nTr
        [pks,locs] = findpeaks(cxs{icond}.filt.IGX(2,:, iTr));
        
        for icomp = 1:2 % SGX/IGX
            
            epochs{icond}.filt.(comp{icomp}) = [];
            epochs{icond}.raw.(comp{icomp}) = [];
            
            for i = 1:param.N
                
                for ipk = 1:numel(pks)
                    if locs(ipk) <= wnsz/2 || locs(ipk) >= numel(t) - wnsz/2
                        continue
                    end
                    for ifld = {'raw' 'filt'}
                        epochs{icond}.(char(ifld)).(comp{icomp})(:,:,end+1) = cxs{icond}.(char(ifld)).(comp{icomp})(:,locs(ipk)-0.5*wnsz:locs(ipk)+0.5*wnsz,iTr);
                    end
                end
                
            end
        end
    end
    
    epochtime = (-wnsz/2:wnsz/2)./2e2;
    scfac = 2e2;
    figure(fig1)
    
    if icond ~= 2
        subplot(2,3,icond)
        plot(epochtime, mean(epochs{icond}.raw.SGX,3) - (1:param.N)'.*scfac, 'k')
        imagesc(epochtime, 1:param.N, mean(epochs{icond}.raw.SGX,3))
        for i = 1:param.N
            [~, pks] = findpeaks(mean(epochs{icond}.filt.SGX(i,:,:),3), epochtime);
            pktm(i) = pks(find(pks > -0.08 & pks < 0.08,1,'last'));
        end
        hold on
        scatter(pktm, 1:param.N, 'wx');
        title(sprintf('SGX %s [%02d %02d] Hz', conds{icond}, round(min(caxis)), round(max(caxis))));
        
    end
    
    subplot(2,3,3+icond)
    imagesc(epochtime, 1:param.N, mean(epochs{icond}.raw.IGX,3))
    
    for i = 1:param.N
        [~, pks] = findpeaks(mean(epochs{icond}.filt.IGX(i,:,:),3), epochtime);
        pktm(i) = pks(find(pks > -0.08 & pks < 0.08,1,'last'));
    end
    hold on
    scatter(pktm, 1:param.N, 'wx');
    title(sprintf('IGX %s [%02d %02d] Hz', conds{icond}, round(min(caxis)), round(max(caxis))));
    
    
    if icond ~= 2
        figure(fig2)
        
        subplot(1,2,double(icond > 2) + 1)
        hold on
        errorbar(1:param.N, mean(mean(apow{icond}.IGX,3),2)', std(mean(apow{icond}.IGX, 3),[], 2), 'k')
        errorbar(1:param.N, mean(mean(apow{icond}.SGX,3),2)', std(mean(apow{icond}.SGX, 3),[], 2), 'r')
    end
end

figure(fig1)
for isp = [1 3 4 5 6]
    subplot(2,3,isp)
    xline(0, 'w--', 'LineWidth', 2)
    colormap(brewermap(64,'YlOrBr'))
    axis off
end

figure(fig2)
for isp = 1:2
    subplot(1,2,isp)
    xlim([-0.5 param.N+1.5])
    xticks(1:param.N)
    ylim([-10 120])
    errorbar(1:param.N, mean(mean(apow{2}.IGX,3),2)', std(mean(apow{2}.IGX, 3),[], 2), 'LineStyle', '--', 'Color', 'k'); % add stimulus off IGX
    ylabel('Power [a.u.]')
    set(gca, 'FontName', 'Arial')
    legend({'IGIB' 'SGX' 'IGIB (BW)'}, 'AutoUpdate', 'off')
    yline(0, 'k-')
    
end


%%
end

