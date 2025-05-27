function pvn_plotEEGTopo(eeg, varargin)
% topography plot for EEG data

p = inputParser();
p.addParameter('Parameter', 'RMS'); % 'Phase' for phase gradient plot
p.addParameter('PhaseReference', 'Iz');
p.addParameter('TimeWindow', [-Inf Inf]);
p.addParameter('PlotType', '3D');
p.addParameter('ROIMask', {[-0.2 0.2] [-0.35 0.25]});
p.addParameter('FigHandle', []);

p.parse(varargin{:});

if isempty(p.Results.FigHandle)
    fh = gcf;
else
    fh = p.Results.FigHandle;
end

if isstruct(eeg) && nargin > 1
    tIdx = eeg.t >= p.Results.TimeWindow(1) & eeg.t < p.Results.TimeWindow(2);
    
    st = [];
    
    switch p.Results.Parameter
        case 'RMS'
            avg = mean(squeeze(rms(eeg.eeg(:,tIdx,:), 2)),2);
        case {'Phase', 'Power'}
            freq = [7 13];
            ref = find(strcmpi(eeg.label, p.Results.PhaseReference));
            bpFilt = designfilt('bandpassfir', 'StopbandFrequency1', 0.85*freq(1), 'PassbandFrequency1', freq(1),...
                'PassbandFrequency2', freq(2), 'StopbandFrequency2', 1.25*freq(2),...
                'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 1/diff(eeg.t(1:2)));
            [phi, pow] = getPhiPow(eeg.eeg, bpFilt);
            phi = phi./phi(ref,:,:);
            powmask = mean(mean(pow(:,tIdx,:),3),2);
            
            if strcmp(p.Results.Parameter, 'Phase')
                avg = wrapTo2Pi(angle(squeeze(mean(mean(phi(:,tIdx,:), 2),3))));
                st.mask = (powmask./max(powmask)).^0.5;
            else
                avg = powmask;
                st.mask = ones(size(powmask));
            end
            
    end
    
    
    st.dimord = 'chan_time';
    st.avg = avg;
    st.label = eeg.label;
    st.time = 1;
    
    sens = ft_read_sens('standard_1005.elc');
    [~, idx] = ismember(sens.label, eeg.label);
    idx = idx(idx > 0);
    
    idx2 = ismember(eeg.label, sens.label);
    avg(~idx2) = [];
    
    switch p.Results.PlotType
        case '3D'
            ft_plot_topo3d(sens.chanpos(idx,:), avg);
            alpha 0.8
            hold on
            ft_plot_sens(sens);
            midLineIdx = find(ismember(sens.label(idx), {'Oz' 'POz' 'Pz' 'CPz' 'Cz' 'FCz' 'Fz'}));
            plot3(sens.chanpos(midLineIdx,1), sens.chanpos(midLineIdx,2), sens.chanpos(midLineIdx,3), 'w');
            view([0 30])
            
        case '2D'
            cfg = [];
            cfg.layout = 'biosemi64.lay';
           
            cfg.parameter = 'avg';
            cfg.style = 'straight';
            cfg.interactive = 'no';
            cfg.comment = 'no';
            cfg.figure = fh;
            ft_topoplotER(cfg, st);
           
            if strcmpi(p.Results.Parameter, 'Phase')
                colormap(hsv)
                caxis([0 2*pi])
                ax = gca;
                xWin = mean(ax.Children(end).XData(:)) + p.Results.ROIMask{1};
                yWin = mean(ax.Children(end).YData(:)) + p.Results.ROIMask{2};
                mask = ax.Children(end).XData > xWin(1) & ax.Children(end).XData < xWin(2);
                mask = mask & ax.Children(end).YData > yWin(1) & ax.Children(end).YData < yWin(2);
                ax.Children(end).CData(~mask) = nan;
            end
                  
    end
    
end

end
%%
function [phi, pow] = getPhiPow(raw, bpFilt)
[nElec, nTm, nTr] = size(raw);
phi = nan(size(raw));

try
    for itrial = 1:nTr
        for iElec = 1:nElec
            eegFilt(iElec,:,itrial) = filtfilt(bpFilt, raw(iElec,:,itrial));
        end
    end
catch errmsg
    if strcmpi(errmsg.identifier, 'signal:filtfilt:InvalidDimensionsDataShortForFiltOrder')
        rawPad = cat(2,raw,raw,raw);
        eegFilt = nan(nElec, 3*nTm, nTr);
        for itrial = 1:nTr
            for iElec = 1:nElec
                eegFilt(iElec,:,itrial) = filtfilt(bpFilt, rawPad(iElec,:,itrial));
            end
        end
        eegFilt = eegFilt(:,nTm+1:2*nTm,:);
    else
        throw(errmsg)
    end
end

for iTr = 1:nTr
    for iElec = 1:nElec
        h = hilbert(eegFilt(iElec,:,iTr));
        phi(iElec,:,iTr) = exp(1i.*angle(h));
        pow(iElec,:,iTr) = abs(h);
    end
end
end