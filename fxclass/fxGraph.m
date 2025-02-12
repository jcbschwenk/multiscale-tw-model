classdef fxGraph < matlab.mixin.Copyable
    properties
        Nodes fxNode
        Curr fxCurr
        Labels cell = {}
        CurrLabels cell = {}
        Syn fxSyn
        FI function_handle = @(I,lm,b,c)(lm.*I-b)./(1-exp(-c.*(lm.*I-b)));
        spkIConv function_handle = @(I) I; % conversion mean-field current -> spk neuron current
        WCurr double % weights curr -> nodes
        name char = ''
    end
    methods
        function obj = fxGraph(varargin)
            obj.Nodes = fxNode.empty();
            obj.Curr = fxCurr.empty();
        end
        %%
        function obj = addNode(obj, node, name)
            if isa(node, 'fxNode')
                obj.Nodes(end+1) = node;
                if nargin < 3
                    name = ['N' num2str(numel(obj.Nodes))];
                end
                obj.Labels{end+1} = name;
                obj.Syn(1:obj.numNodes,obj.numNodes) = fxSyn; % fill with empty synapse objects
                obj.WCurr(1:obj.numCurr,obj.numNodes) = 0;
            else
                error('Node invalid');
            end
        end
        %%
        function obj = addCurr(obj, curr, name)
            assert(isa(curr, 'fxCurr'), 'Curr object invalid')
            obj.Curr(end+1) = curr;
            obj.CurrLabels{end+1} = name;
            obj.WCurr(obj.numCurr,1:obj.numNodes) = 0;
            
        end
        %%
        function obj = clearCurr(obj, name)
            if nargin < 2 || isempty(name)
               name = obj.CurrLabels;
            elseif ~iscell(name)
                name = {name};
            end
            obj.Curr(obj.getCurrIdx(name)) = [];
            obj.CurrLabels(obj.getCurrIdx(name)) = [];
        end
        %%
        function obj = connect(obj, node1, node2, varargin)
            
            assert(ismember(node1, obj.Labels) || ismember(node1, obj.CurrLabels),...
                'Bad node references')
            assert(ismember(node2, obj.Labels) || ismember(node2, obj.CurrLabels),...
                'Bad node references')
            
            % parse input as property list for fxSyn:
            p = inputParser;
            
            p.addParameter('Update', false); 
            % if true, use existing synapse, only changing requesting parameters
            
            for iProp = cell(properties('fxSyn'))'
                p.addParameter(char(iProp), []);
            end
             
            % check if fxSyn object is passed:
            synObjPassed = isa(varargin{1}, 'fxSyn');
            if synObjPassed
                synObj = varargin{1};
                varargin = varargin(2:end);
            else
                synObj = fxSyn;
            end
            
            p.parse(varargin{:})
            
        
            if p.Results.Update
                assert(~synObjPassed, 'Cannot update if synapse object passed')
                synObj = obj.getSyn(node1, node2);
                assert(~isempty(synObj), 'Cannot update, no existing connection')
                assert(sum(synObj == obj.Syn(:)) == 1, 'Cannot update, target synapse has multiple instances in graph.')
            end
            
            
            
            % parse anything passed after the fxSyn object:
            % (this overwrites any existing param specifications in the obj)
            parseRez = p.Results;
            if ~isempty(varargin)
                for specParam = cell(properties('fxSyn'))'
                    if ~isempty(parseRez.(char(specParam)))
                        synObj.(char(specParam)) = parseRez.(char(specParam));
                    end
                end
            end
            
            if isscalar(node1)
                assert(node1 > 0 && node1 <= obj.numNodes, 'Bad node reference');
                idx1 = node1;
            else
                idx1 = obj.getNodeIdx(node1);
            end
            
            if isscalar(node2)
                assert(node2 > 0 && node2 <= obj.numNodes, 'Bad node reference');
                idx2 = node2;
            else
                idx2 = obj.getNodeIdx(node2);
            end
            
            if ~idx1
                idx1 = obj.getCurrIdx(node1);
                assert(~isempty(idx1), 'Bad node reference')
                assert(~isempty(parseRez.W), 'Must pass weight argument if connecting current to node')
                obj.WCurr(idx1, idx2) = parseRez.W;
            else
                assert(idx2 > 0, 'Bad node reference')
                assert(isa(synObj, 'fxSyn'), 'Bad synapse argument')
                obj.Syn(idx1,idx2) = synObj;
            end
        end
        %% 
        function obj = setWeight(obj, fromNode, toNode, newWeight)
            % Separate function to update weights only. If w is zero
            % (request disconnect), either input node may be left empty to
            % cut all incoming / outgoing connections.
            
            if isempty(fromNode)
                fromNode = char;
            end
            if isempty(toNode)
                toNode = char;
            end
              
            nd1 = obj.getNodeIdx(fromNode);
            nd2 = obj.getNodeIdx(toNode);
            
            assert(~(newWeight ~= 0 && ~all([nd1 nd2])),...
                'Open connection definitions only valid for disconnect (w=0)')
            
            if ~nd1
               nd1 = 1:obj.numNodes;
            elseif ~nd2
               nd2 = 1:obj.numNodes;
            end
            

            [obj.Syn(nd1, nd2).W] = deal(newWeight);
            
        end
        %% 
        function out = W(obj, varargin)
            % return weight matrix
            out = obj.Syn.getParam;
            if ~isempty(varargin)
                out = out(varargin{:});
            end
        end
        %% 
        function out = T(obj, varargin)
            % return weight matrix
            [~, out] = obj.Syn.getParam;
            if ~isempty(varargin)
                out = out(varargin{:});
            end
        end
        %%
        function [outNode, outLabel, outIdx] = getRecipients(obj, node)
            myIdx = obj.getNodeIdx(node);
            W = obj.Syn.getParam;
            outIdx = find(abs(W(myIdx,:))>0);
            outNode = obj.Nodes(outIdx);
            outLabel = obj.Labels(outIdx);
        end
         %%
        function [outNode, outLabel, outIdx] = getSenders(obj, node)
            myIdx = obj.getNodeIdx(node);
            W = obj.Syn.getParam;
            outIdx = find(abs(W(:,myIdx))>0);
            outNode = obj.Nodes(outIdx);
            outLabel = obj.Labels(outIdx);
        end
        %%
        function out = run(obj, varargin)
            
            p = inputParser;
            p.addParameter('Dur', []);
            p.addParameter('NumTrials', []);
            p.addParameter('Trial', []);
            p.parse(varargin{:});
            
            dur = p.Results.Dur;
            nTr = p.Results.NumTrials;
            trID = p.Results.Trial;

            refCurr = obj.Curr(find(~isNull(obj.Curr), 1, 'first'));
            
            if isempty(dur)
                dur = refCurr.dur;
            end
            
            if isempty(nTr) && ~isempty(trID)
                nTr = 1;
            elseif isempty(nTr)
                nTr = refCurr.numTrials;
            end
            
            if nTr > 1
               for iTrial = 1:nTr
                  out(iTrial) = obj.run('Dur', dur, 'NumTrials', 1, 'Trial', iTrial);
               end
               return
            elseif isempty(trID)
                trID = 1;
            end
            
            % External currents:
            if obj.numCurr == 0 || all(isNull(obj.Curr))
                extCurr = zeros(1,dur);
            else
                idx = find(~isNull(obj.Curr));
                tmp = cat(1, obj.Curr(idx).I);
                extCurr = zeros(obj.numCurr, size(tmp,2), size(tmp,3));
                extCurr(idx,:,:) = tmp;
                extCurr = extCurr(:,:,trID);
            end
            
            
            % NODE PARAMETERS:
            lm = [obj.Nodes.lm]; % FI curve gain (lambda)
            b = [obj.Nodes.b]; %
            c = [obj.Nodes.c]; %

            IBase = [obj.Nodes.Ib]; % base current
            
            sigN = [obj.Nodes.e]; % noise amplitude
            tauN = [obj.Nodes.eTD]; % noise tau
                       
            % Synapse-specific parameters:
            [W, T, Tau, G] = obj.Syn.getParam();
            
            % Spiking children:
            [iz_a, iz_b, iz_c, iz_d, iz_e, neur2Node] = getSpkChildParam(obj.Nodes);
            nSpkNeur = numel(neur2Node);
            spkNodeIdx = hasSpk(obj.Nodes);
            if any(spkNodeIdx)
                izN = rand(nSpkNeur, dur) .* iz_e(2,:)' + iz_e(1,:)'; % noise
            end

            % Get unique delays:
            allT = unique(T(:));
            
            % Initialize running variables:
            s = 0 .* W; % synaptic gating variable
            r = zeros(size(W,1), dur); % rate r
            [v, u] = deal(zeros(1, nSpkNeur)); % tmp u,v for spiking neurons
            spkV = nan(nSpkNeur, dur); % output v 
            spkThresh = 30; % in mV
            
            I = zeros(1, size(W,1)); % current I
            n = 0.*I;
            
            for t = 1:dur
                
                % perform spikes:
                if any(spkNodeIdx)
                    spkidx = v >= spkThresh;
                    v(spkidx) = iz_c(spkidx);
                    u(spkidx) = u(spkidx) + iz_d(spkidx);
                end
                
                % Current to rate conversion:
                r(:, t) =  obj.FI(I, lm, b, c);                
                if any(spkNodeIdx)
                    r(spkNodeIdx, t) = obj.getInstantSpkRate(spkidx, neur2Node);
                end
                
                % Get synaptic decay before updating synapses:
                synDecay = s./Tau;
                
                % Update synaptic gating variables, split by delays:
                for iT = 1:numel(allT)
                    if allT(iT) >= t
                        continue
                    end
                    s = s...
                        + (T == allT(iT)) .* G./1e3 .* (1-s) .* r(:, t-allT(iT)); 
                end
                
                % Apply synaptic decay:
                s = s - synDecay;
                
                % Get noise:
                n = n - n./tauN + randn(size(n)).*sqrt(tauN.*sigN.^2); 
                
                % Compute final summed current:
                Iapp = sum(extCurr(:,t) .* obj.WCurr, 1);
                I = sum(W .* s, 1) + Iapp + IBase + n;
                
                % Convert current to input for spiking neurons:
                if any(spkNodeIdx)
                    spkI = obj.spkIConv(I(neur2Node));
                    spkI = spkI + izN(:,t)';
                    
                    % Update spk neurons:
                    v = v + 0.5.*(0.04*v.^2+5*v+140-u+spkI);
                    v = v + 0.5.*(0.04*v.^2+5*v+140-u+spkI);
                    u = u + iz_a.*(iz_b.*v-u);
                    
                    spkV(:,t) = v;
                end
                
            end
            
            
            out.t = refCurr.t;
            out.r = r;
            out.ext = extCurr;
            
            if any(spkNodeIdx)
                spkBool = spkV > spkThresh;
                
                spkBool(:,1:2) = 0; % remove first two timepoints
                out.r(spkNodeIdx,1:2) = 0;
                spkV(spkBool) = 50;
                
                out.v = spkV;
                out.spk = obj.getSparseSpkList(spkBool, neur2Node);
                out.neur2Node = neur2Node;
                out.spkNodeIdx = find(spkNodeIdx);
            else % keep output uniform
                out.v = [];
                out.spk = [];
                out.neur2Node = [];
                out.spkNodeIdx = [];
            end
        end
        %%
        function sOut = getSynCurr(obj, runOut, fromNode, toNode)
            % get synaptic current for a single synapse. runOut should be
            % the output of obj.run, fromNode and toNode are labels. 
            % if fromNode is empty, returns sum of all input to toNode
            
            if isempty(fromNode)
                [~,allSend] = obj.getSenders(toNode);
                for iNd = 1:numel(allSend)
                   sOut(iNd,:) = obj.getSynCurr(runOut, allSend{iNd}, toNode);
                end
                sOut = sum(sOut,1);
                return
            end
            
            nd1 = obj.getNodeIdx(fromNode);
            nd2 = obj.getNodeIdx(toNode);
            
            dur = numel(runOut.t);
            [sW, sT, sTau, sG] = obj.Syn(nd1, nd2);
                      
            v = runOut.v(nd1,:);
            s = 0;
            sOut = nan(1,dur);
            
            for t = 2:dur
                
                decay = s./sTau;
                
                if t > sT
                    s = s + sG./1e3 .* (1-s) .* v(:, t-sT); % saturating
                end
                
                s = s - decay;
                sOut(t) = s;
            end
            sOut = sW .* sOut;
        end
        %%
        function idx = getNodeIdx(obj, label)
            if nargin < 2
               idx = 0;
               return
            end
            [~, idx] = ismember(label, obj.Labels);
            
        end
        %%
        function idx = getCurrIdx(obj, label)
            if nargin < 2
               idx = 0;
               return
            end
            [~, idx] = ismember(label, obj.CurrLabels);
        end
        %%
        function syn = getSyn(obj, fromNode, toNode)
            
            nd1 = obj.getNodeIdx(fromNode);
            nd2 = obj.getNodeIdx(toNode);
            syn = obj.Syn(nd1, nd2);
            
        end
        %%
        function out = numNodes(obj)
            out = numel(obj.Nodes);
        end
        %%
        function out = numCurr(obj)
            out = numel(obj.Curr);
        end
        %%
        function plotFI(obj, nodes, x)
            if nargin < 2 || isempty(nodes)
               nodes = obj.Labels;
            end
            if nargin < 3
                x = -0.2:0.01:4;
            end
            for iNode = 1:numel(nodes)
                nd = obj.Nodes(strcmpi(obj.Labels, nodes{iNode}));
                plot(x, obj.FI(x, nd.lm, nd.rmax));
                hold on
            end
            xlabel('Current')
            ylabel('Firing Rate [Hz]')
            legend(nodes)
            ylim([-30 max(ylim)+50])
            yline(0,'k:')
        end
    end
    %%
    methods (Static)
               function r = getInstantSpkRate(spkBool, id)
            % get instantenous spike rate from the population (no
            % windowing). Output rate is the normalized count per node. 
            
            % TODO: adjust for
            % synchrony w/ transduction function
            
            r = 1e3 .* arrayfun(@(x) mean(double(spkBool(id == x))), unique(id));
            
        end
        %%
        function spk = getSparseSpkList(spkBool, id)
            % returns sparse spike list for each node as cell array.
            % columns of list are: 
            % [spiketime neurID_in_graph neurID_in_node];
            
            nNeur = size(spkBool,1);
            
            spk = [];
            for iNeur = 1:nNeur
                nSpk = find(spkBool(iNeur,:));
                spk = [spk; nSpk' repmat(iNeur,numel(nSpk),1)];
            end
            
            spk = arrayfun(@(x) spk(id(spk(:,2)) == x,:), unique(id), 'UniformOutput', false);
            spk = cellfun(@(x) [x x(:,2)-x(1,2)+1], spk, 'UniformOutput', false);
            
        end
        %% 
    end
    %%
end
