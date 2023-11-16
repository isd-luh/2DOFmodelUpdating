% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef GlobalPattern < Optimizer
    % GLOBALPATTERN is a single-objective deterministic algorithm derived from Optimizer.
    %
    % LITERATURE:
    %   B. Hofmeister, M. Bruns, R. Rolfes (2019) Finite element model updating
    %   using deterministic optimisation: A global pattern search approach,
    %   Engineering Structures (195)
    
    methods
        function this = GlobalPattern()
            
            this.defaultparams = struct('nTrack', 50, ...   % Local minima to track
                'nResolution', 2^24 );                      % Resolution of the grid
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Initialize all states
            nDims = numel(x_lb);
            
            this.states = struct('fBest', Inf, ...                  % Current best objective value
                'viResolution', zeros(nDims), ...  % Allocate space for resolution vector
                'miSampleCache', [], ...           % Cache for sample coordinates
                'vfObjectiveValues', []);          % Cache for objective values
            
            % Initialize resolution per dimension
            % Use given resolution for real dimensions
            this.states.viResolution = ones(1,nDims) * this.params.nResolution;
            % Use actual resolution for integer dimensions
            if ~isempty(vartype)
                for m = 1:nDims
                    if strcmpi(vartype{m}, 'real')
                        this.states.viResolution(m) = this.params.nResolution;
                    elseif strcmpi(vartype{m}, 'int')
                        this.states.viResolution(m) = x_ub(m) - x_lb(m);
                    end
                end
            end
            
            % Initialize step sizes to half the search space
            this.states.viStep = floor((this.states.viResolution+1) / 2);
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % First the algorithm chooses the bases for a variation
            
            % Initialize and use all previous locations initially
            viBases = [this.states.viStep; this.states.miSampleCache];
            
            % If we have enough samples already, only use best ones
            if numel(this.states.miSampleCache) > 0
                if numel(this.states.miSampleCache(:, 1)) > this.params.nTrack
                    % Sort by minima
                    [vfSortedMins, viSortedIndices] = sort(this.states.vfObjectiveValues);
                    viBases = this.states.miSampleCache(viSortedIndices(1:this.params.nTrack),:);
                end
            end
            
            % Now the bases are fixed and the samples are created by
            % Variation along the problem axes
            
            miSamples = []; % stores sampling positions
            
            for viBaseT = viBases'
                viBase = viBaseT';
                
                % Move base back to even grid
                viBase = round(viBase ./ this.states.viStep) .* this.states.viStep;
                
                % Confine bases to be inside the search space
                viBase = max(min(viBase, this.states.viResolution ), zeros(size(this.states.viResolution)));
                
                % Take care of initialization by adding the sample itself
                % This will be filtered out when already in cache
                miSamples = [miSamples ; viBase];
                
                % Wiggle each dimension up and down
                for iDim = [1 : numel(viBase)]
                    viSample = viBase;
                    iBase = viSample(iDim);
                    % Wiggle up and sample
                    viSample(iDim) = min(iBase + this.states.viStep(iDim), this.states.viResolution(iDim));
                    % Save to sample list
                    miSamples = [miSamples ; viSample];
                    % Wiggle down and sample
                    viSample(iDim) = max(iBase - this.states.viStep(iDim), 0);
                    % Save to sample list
                    miSamples = [miSamples ; viSample];
                end
            end
            
            % Deduplicate samples using cache
            if numel(this.states.miSampleCache)
                miDedupSamples = setdiff(miSamples, this.states.miSampleCache, 'rows');
            else
                miDedupSamples = miSamples;
            end
            
            % Filter samples and output
            samples = [];
            for viSampleT = miDedupSamples'
                viSample = viSampleT';
                
                % We did not sample here before, put this sample to the cache
                this.states.miSampleCache = [this.states.miSampleCache; viSample];
                % Transform samples to floating point
                samples = [samples; viSample ./ (this.states.viResolution) .* (x_ub-x_lb) + x_lb];
            end
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Append results to objective value cache
            this.states.vfObjectiveValues = [this.states.vfObjectiveValues; objectiveValues];
            
            % Get minimum objective value
            fMin = min(objectiveValues);
            
            % Check if we have a new minimum
            if fMin < this.states.fBest
                this.states.fBest = fMin;
                terminate = false;
                return
            end
            
            % Find out which dimension has the largest step
            [iMaxStep, iMaxDim] = max(this.states.viStep);
            
            terminate = true;
            
            % If not lowest resolution reached
            if iMaxStep > 1
                terminate = false;
                
                % Reduce largest step by factor 2
                this.states.viStep(iMaxDim) = floor(this.states.viStep(iMaxDim)/ 2);
            end
            
        end
        
    end
end