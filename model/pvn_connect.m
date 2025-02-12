function [g, param] = pvn_connect(g, param)

isdual = any(contains(param.NodeLabels, '_R'));

conn.synTypes = {'E_fast'   'E_slow'    'I_fast',  'I_slow',    'IGself'};
conn.synTau =   [20         120         3           120          1];
conn.synGamma = [0.8        0.5         0.8         0.1         0.8];

conn.interAreaDelay = 12;
conn.FDnDelay = ceil(conn.interAreaDelay/2);
conn.FUpDelay = floor(conn.interAreaDelay/2);

conn.connTypes = {'Loc' 'FFw' 'FBk' 'FDn' 'FUp' 'FBkUp' 'PulSync'};

% LOCAL CONNECTIONS
conn.Loc.rule = 0;
conn.Loc.delay = 0;
conn.Loc.def = {...
    {'L4X', 'SGX', 'E_fast', 1},...
    {'L4X', 'L4IN', 'E_fast', 1},...
    {'L4IN', 'L4X', 'I_slow', -1.2},...
    {'ILSt', 'ILStIN', 'E_fast', 1},...
    {'ILStIN', 'ILSt', 'I_slow', -1.2},...
    {'SGX', 'SGE', 'E_fast', 1},...
    {'SGIN', 'SGX', 'I_fast', -2},...
    {'SGX', 'IGX', 'E_fast', 1},...
    {'IGX', 'IGX', 'IGself', 0.7},...
    {'IGIN', 'IGX', 'I_fast', 0},...
    };

% FEEDFORWARD CONNECTIONS
conn.FFw.rule = 1;
conn.FFw.delay = conn.interAreaDelay;
conn.FFw.def = {...
    {'ILSt', 'L4X', 'E_fast', 2},...
    {'SGE', 'L4X', 'E_fast', 2},...
    };

% FEEDBACK CONNECTIONS
conn.FBk.rule = -1;
conn.FBk.delay = conn.interAreaDelay;
conn.FBk.def = {...
    {'ILPr', 'IGX', 'E_fast', 1.5},...
    {'IGX', 'IGX', 'E_slow', 1},...
    {'IGX', 'SGIN', 'E_fast', 1.5},...
    };

% FEED-DOWN CONNECTIONS (TO PUL)
conn.FDn.rule = 0;
conn.FDn.delay = conn.FDnDelay;
conn.FDn.def = {...
    {'IGX', 'Pul', 'E_fast', 1},...
    };

% FEED-UP CONNECTIONS (FROM PUL)
conn.FUp.rule = 1;
conn.FUp.delay = conn.FUpDelay;
conn.FUp.def = {...
    {'Pul', 'L4X', 'E_fast', 0},...
    };

% FEEDBack-UP CONNECTIONS (FROM PUL)
conn.FBkUp.rule = 0;
conn.FBkUp.delay = conn.FUpDelay;
conn.FBkUp.def = {...
    {'Pul', 'IGIN', 'E_fast', 0},... % needs local IGIN->IGX connection! (above)
    };


% IL TO PUL CONNECTIONS (single delay, i.e. for synchronized pulvinar)
conn.PulSync.rule = 1:param.N;
conn.PulSync.delay = conn.interAreaDelay;
conn.PulSync.def = {...
    {'ILSt', 'Pul', 'E_fast', 0},...
    };


% expand all connections to include hemisphere suffix if dual-stream:
if isdual
    hemis = {'L' 'R'};
    for cType = conn.connTypes
            newDef = cellfun(@(h)...
                cellfun(@(d) {[d{1} '_' h] [d{2} '_' h] d{3} d{4}},...
                conn.(char(cType)).def, 'UniformOutput', false),...
                hemis, 'UniformOutput', false);
            conn.(char(cType)).def = [newDef{:}];
    end
end


param.Syn = fxSyn.empty;
param.SynLabels = {};

for cType = conn.connTypes
    
    
    connSet = conn.(char(cType));
    N = numel(connSet.def);
    
    for iConn = 1:N
        
        [node1, node2] = connSet.def{iConn}{1:2};
        
        synType = find(strcmpi(conn.synTypes, connSet.def{iConn}{3}));
        w = connSet.def{iConn}{4};
        
        param.Syn(end+1) = fxSyn(...
            'W', w,...
            'T', connSet.delay,...
            'Tau', conn.synTau(synType),...
            'G', conn.synGamma(synType)...
            );
        % -> this is a handle to the synapse in the final graph!
        % all connections using this synapse definition will be updated together
        
        param.SynLabels{end+1} = [char(cType) '_' node1 '_' node2];
        
        pvn_nodePairLayerConnect(g, param, connSet.rule, node1, node2, param.Syn(end));
        
    end
    
end

conn.stimWeight = 1;
conn.priorWeight = 1;

if ~isdual
    g.connect('Stim', 'ILSt0', 'W', conn.stimWeight);
    g.connect('Prior', ['ILPr' num2str(param.N+1)], 'W', conn.priorWeight);
else
    for h = hemis
        % hemisphere-specific input:
        g.connect(['Stim_' char(h)], ['ILSt_' char(h) '0'], 'W', conn.stimWeight);
        g.connect(['Prior_' char(h)], ['ILPr_' char(h) num2str(param.N+1)], 'W', conn.priorWeight);
        
        % non-specific input (connects to both hemispheres):
        g.connect('Stim_M', ['ILSt_' char(h) '0'], 'W', conn.stimWeight);
        g.connect('Prior_M', ['ILPr_' char(h) num2str(param.N+1)], 'W', conn.priorWeight);
    end
end

param.conn = conn;

end