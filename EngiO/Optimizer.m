% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef Optimizer < handle
    % OPTIMIZER is base class for all optimizers. The optimization problem
    % should be stated as follows:
    %   minimize f(x) for x_lb <= x <= x_ub
    % Supports parallel computing (NumWorkers>1)
    
    properties
        % Struct defined in base class
        
        % Global options of optimizer
        options = struct('maxEvals', 2000,... % Maximum count of objective function evaluations (useful to saturate overall optimization time)
            'maxIters', Inf,... % Maximum total count of iterations
            'maxItersNoChange', Inf,... % Maximum count of iterations in succession with no change of objective
            'numWorkers', 1,... % Number of workers
            'outputStatus', true,... % Select whether status should be displayed (true) in command line or not (false)
            'rndFact', 1.0e-3,... % Factor by which the values are rounded
            'saveStates', false); % Select whether the current states shall be saved to a file in each generation
        
        % Structs to be defined in derived classes
        
        params = struct; % Parameters of optimizer
        states = struct; % States of optimizer
        
        defaultoptions = struct; % Default options of optimizer
        defaultparams = struct; % Default parameters of optimizer
    end
    
    methods(Access = public)
        %% Public methods to be implemented in derived classes
        
        function samples = initialize(this, nObj, vartype, x_lb, x_ub) %#ok<STOUT, INUSD>
            % INITIALIZE must initialize all samples before starting the iteration.
            %
            % INPUT:
            % - nObj (integer scalar)
            %   Number of objectives
            % - vartype (string|cell of strings)
            %   Variable types of design variables; valid inputs are 'real'
            %   or 'int'; if empty all design variables are 'real'
            % - x_lb (vector of doubles)
            %   Lower boundaries of design variables
            % - x_ub (vector of doubles)
            %   Upper boundaries of design variables
            %
            % OUPUT:
            % - samples (m by n matrix of doubles)
            %   Sampled design varibles, where m is number of samples
            %   and n number of design variables per sample
            %
            % See also OPTIMIZE.
        end
                
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub) %#ok<STOUT, INUSD>
            % GENERATESAMPLES must generate all samples for the current iteration.
            %
            % INPUT:
            % - nObj (integer scalar)
            %   Number of objectives
            % - vartype (string|cell of strings)
            %   Variable types of design variables; valid inputs are 'real'
            %   or 'int'; if empty all design variables are 'real'
            % - x_lb (vector of doubles)
            %   Lower boundaries of design variables
            % - x_ub (vector of doubles)
            %   Upper boundaries of design variables
            %
            % OUPUT:
            % - samples - (m by n matrix of doubles)
            %   Sampled design varibles, where m is number of samples
            %   and n number of design variables per sample
            %
            % See also OPTIMIZE.
        end
        
        function processResults(this, samples, objectiveValues) %#ok<INUSD>
            % PROCESSRESULTS processes the results of each iteration.
            %
            % INPUT:
            % - samples (m by n matrix of doubles)
            %   Sampled design varibles, where m is number of samples
            %   and n number of design variables per sample
            % - objectiveValues (m by n matrix of doubles)
            %   Objective function values of current samples, where m is
            %   number of samples and n number of objectives
            %
            % See also OPTIMIZE.
        end
                
        %% Public methods to be used in derived classes
        
        function varargout = optimize(this, objective, vartype, x_lb, x_ub, options, params, states_init)
            % OPTIMIZE performs the optimization procedure.
            %
            % SYNTAX:
            % varargout = OPTIMIZE(this, objective, vartype, x_lb, x_ub)
            % varargout = OPTIMIZE(this, objective, vartype, x_lb, x_ub, options, params, states_init)
            % states_temp = OPTIMIZE(___)
            % [x_pareto, f_pareto]= OPTIMIZE(___)
            % [x_pareto, f_pareto, numEvals, mSamples, vResults]= OPTIMIZE(___)
            % 
            % INPUT:
            % - objective (function handle)
            %   Objective function
            % - vartype (string|cell of strings)
            %   Variable types of design variables; valid inputs are 'real'
            %   or 'int'; if empty all design variables are 'real'
            % - x_lb (vector of doubles)
            %   Lower boundaries of design variables
            % - x_ub (vector of doubles)
            %   Upper boundaries of design variables
            % - options (struct)
            %   Specific options of optimizer (optional)
            % - params (struct)
            %   Specific parameters of optimizer (optional)
            % - states_init (struct)
            %   Initial optimization states (optional)
            %
            % OUTPUT:
            % - varargout (cell|vectors|multiple matrices)
            %   Output arguments;
            %   if nargout==1, all states for all generations;
            %   if nargout==2, optimal variables and objective values of last iteration.
            %   if nargout==5, optimal variables and objective values of 
            %   last iteration, numEvals, all samples and corresponding 
            %   objective values.

            % Check number of inputs and outputs
            narginchk(5, 8);
            nargoutchk(1, 5);
            
            if nargout==3 || nargout==4
                error('Number of outputs must be 1, 2, or 5.')
            end
            
            % Some simple checks on validity of inputs
            if iscell(objective)
                objective = objective{1};
            end
            
            if ~isa(objective, 'function_handle')
                error('The objective must be a function handle.')
            end
            
            % In case of vartype dimension of one, transform it to cell array
            if ~isempty(vartype)
                if numel(vartype)==1
                    if ~iscell(vartype)
                        vartype ={vartype};
                    end
                end
            end
            
            if ~isempty(vartype)
                if ~iscell(vartype) || ~isvector(vartype) || length(vartype)~= length(x_lb)
                    error('Input vartype must be empty or a cell array comprising length(x_lb) elements.')
                else
                    for m = 1:length(vartype)
                        if ~strcmpi(vartype{m}, 'real') && ~strcmpi(vartype{m}, 'int')
                            error('All vartypes must be either real or int.')
                        end
                    end
                end
            end
            
            if ~isnumeric(x_lb) || ~isvector(x_lb)
                error('Input x_lb must be a numeric vector.')
            end
            
            if ~isnumeric(x_ub) || ~isvector(x_ub)
                error('Input x_ub must be a numeric vector.')
            end
            
            if length(x_lb)~= length(x_ub)
                error('x_lb and x_ub must have the same length.')
            end
            
            if any(x_ub<= x_lb)
                error('x_lb must be smaller than x_ub.')
            end
            
            % Find out the number of objective functions
            if nargin(objective)==1
                nObj = objective([]);
            else
                nObj = objective([], 1);
            end
            
            % Set options if desirable
            this.options = this.setOptions(this.options, this.defaultoptions);
            if exist('options', 'var')
                this.options = this.setOptions(this.options, options);
            end
            
            % Set parameters if desirable
            this.params = this.defaultparams;
            if exist('params', 'var')
                this.params = this.setOptions(this.params, params);
            end
            
            % Initialize sample and result history
            mSamples = [];
            vResults = [];
            
            % Initialize states
            if nargin<8
                initialize(this, nObj, vartype, x_lb, x_ub);
            elseif nargin==8
                if ~isstruct(states_init)
                    error('Input states_init must be either a struct or a string.')
                end
                
                % Restore states from data
                this.states = states_init;
                mSamples = states_init.mSamples;
                vResults = states_init.vResults;
            end
            
            % Allocate parallel pool if necessary
            if this.options.numWorkers>1
                if ~exist ('OCTAVE_VERSION', 'builtin')
                    pool = gcp;
                    if pool.NumWorkers<this.options.numWorkers % Reserve appropriate number of workers if necessary
                        delete(pool);
                        parpool('local', this.options.numWorkers); % Profile 'local' is used by default
                    end
                end
       
            end
            
            % Start timer
            tic;
            
            % Allocate local variables
            numEvals = 0;
            numIters = 0;
            numItersNoChange = 0;
            objectiveBest = inf;
            objectiveBestChange = 0;
            states_temp ={};
            termFlag = 0;
            
            % Output initial status to command window
            if this.options.outputStatus
                outputStatus(this, numIters, numEvals, objectiveBest, objectiveBestChange, termFlag,[]);
            end
            
            % Repeat iteration until convergence is reached (termFlag ~=0)
            while ~termFlag
                
                % Iteration counter
                numIters = numIters+1;
                
                % Generate sample set
                samples = generateSamples(this, nObj, vartype, x_lb, x_ub);
                
                % Test if more samples than allowed were requested
                if numEvals+size(samples, 1) > this.options.maxEvals
                    nSamplesAllowed = this.options.maxEvals - numEvals;
                    
                    % Truncate number of samples to fill up to maximum
                    samples = samples(1:nSamplesAllowed, :);
                end
                
                % Round integer design variables
                if ~isempty(vartype) && ~isempty(samples)
                    for m = 1:length(vartype)
                        if strcmpi(vartype{m}, 'int')
                            samples(:, m) = round(samples(:, m));
                        end
                    end
                end
                
                mSamples = [mSamples; samples]; %#ok<AGROW>
                
                % Evaluate objectives
                objectiveValues = this.evalObjective(objective, nObj, samples);
                numEvals = numEvals+size(samples, 1);
                vResults = [vResults; objectiveValues]; %#ok<AGROW>
                
                % Get pareto frontier
                frontier = this.paretoSort(vResults);
                x_pareto = mSamples(frontier,:);
                f_pareto = vResults(frontier,:);
                
                if nObj == 1
                    minObj = min(vResults);
                else
                    % determine reference point for hypervolume metric
                    indices = isfinite(sum(vResults, 2));
                    
                    if any(indices) > 0
                        vMax = max(vResults(indices,:), [], 1);
                        minObj = - this.hyperVolumeMetric(f_pareto, vMax);
                    else
                        minObj = inf;
                    end
                end
                
                objectiveBestChange = objectiveBest-minObj;
                objectiveBest = minObj;
                
                % Set objectiveBest if necessary
                if minObj<objectiveBest
                    numItersNoChange = 0;
                else
                    numItersNoChange = numItersNoChange+1;
                end
                
                % Output current status to command window
                if this.options.outputStatus
                    outputStatus(this, numIters, numEvals, objectiveBest, objectiveBestChange, termFlag,[]);
                end
                
                % Check termination conditions
                if numItersNoChange == this.options.maxItersNoChange % No change of objective (convergence)
                    termFlag = 1;
                end
                
                if numIters == this.options.maxIters % Maximum number of iterations
                    termFlag = 2;
                end
                
                if numEvals>= this.options.maxEvals % Maximum number of objective evaluations
                    termFlag = 3;
                end
                
                % Process results
                if termFlag == 0
                    terminate = this.processResults(samples, objectiveValues);
                    
                    if terminate
                        termFlag = 4;
                    end
                end
                
                % Store states to temporary struct
                states_temp{numIters}= this.states;
                
                % Store sampling history
                states_temp{numIters}.samples = samples;
                states_temp{numIters}.objectiveValues = objectiveValues;
                
                % Save temporary states to file if desirable
                if this.options.saveStates == true
                    save('states.mat', 'states_temp');
                elseif ischar(this.options.saveStates) && strfind(this.options.saveStates, '.mat')
                    save(this.options.saveStates, 'states_temp');
                end
                
            end
            
            % Output last status message to command window
            if this.options.outputStatus
                this.outputStatus(numIters, numEvals, objectiveBest, objectiveBestChange, termFlag, toc);
            end
            
            % Check for number of outputs
            if nargout==1
                varargout{1}= states_temp;
            elseif nargout==2
                varargout{1}= x_pareto;
                varargout{2}= f_pareto;
            elseif nargout==5
                varargout{1}= x_pareto;
                varargout{2}= f_pareto;
                varargout{3}= numEvals;
                varargout{4}= mSamples;
                varargout{5}= vResults;
            end
            
        end
        
        function indices = paretoSort(this, f)
            % PARETOSORT sorts result according to their dominance.
            %
            % INPUT:
            % - f (vector|matrix)
            %   Objective function values
            %
            % OUTPUT:
            % - indices (scalar|vector)
            %   Indices of Pareto optimal solutions
            %
            % See also OPTIMIZE, HYPERVOLUMEMETRIC.
            
            % Choose sorting strategy according to number of objectives
            switch size(f, 2)
                case 1
                    [~, indices] = min(f);
                    indices = indices(1);
                    return
                case 2
                    % Sort by first objective value
                    [~, sortIndex] = sort(f(:, 1));
                    
                    % Indicates the min. value of second objective
                    currMin = Inf;
                    
                    indices = sortIndex([1; diff(cummin(f(sortIndex,2)))<0]>0);
                otherwise
                    [~, sorted] = sort(f(:, 1));
                    indices = this.paretoSortRecursive(f, sorted);
            end
            
        end
        
        function minIndices = paretoSortRecursive(this, f, indices)
            % PARETOSORTRECURSIVE sorts results with more than 2 objectives.
            %
            % INPUT:
            % - f (vector|matrix)
            %   Objective function values
            % - indices (vector)
            %   Indices sorted by first objective
            %
            % OUTPUT:
            % - minIndices (vector)
            %   Indices of Pareto optimal solutions
            %
            % See also PARETOSORT.
            
            nDim = size(f, 2);
            nSize = numel(indices);
            if nSize == 1
                minIndices = indices;
                return
            end
            
            % Divide and conquer
            nHalf = round(nSize/2);
            indices1 = this.paretoSortRecursive(f, indices(1:nHalf));
            indices2 = this.paretoSortRecursive(f, indices(nHalf+1:end));
            
            minIndices = [];
            
            % Check both sets against each other
            for iPass = 1:2
                indicesA = indices1;
                indicesB = indices2;
                
                if iPass == 2
                    indicesA = indices2;
                    indicesB = indices1;
                end
                
                % For each point in A, check against all points in B
                for iPoint = indicesA
                    dominated = ones(numel(indicesB), 1);
                    
                    for iDim = 1:nDim
                        dominated = and(dominated, f(indicesB, iDim) < f(iPoint, iDim));
                    end
                    if ~ any(dominated)
                        minIndices = [minIndices, iPoint];
                    end
                end
            end
            
        end
        
        function V= hyperVolumeMetric(this, points, ref)
            % HYPERVOLUMEMETRIC recursively calculates hypervolume for multi-objective optimization.
            %
            % INPUT:
            % - points (m by n matrix of doubles)(f_pareto, max(vResults, [], 1));
            %   Objective values, where m is number of solutions and n is 
            %   number of objective functions
            % - ref (m by n matrix of doubles)
            %   Reference solution to compare with, where m is number of 
            %   reference solutions and n is number of objective functions
            %
            % OUTPUT:
            % - V (scalar)
            %   Value of hypervolume metric
            %
            % See also OPTIMIZE.
            
            % Pareto sort the points
            paretoIndices = this.paretoSort(points);
            
            paretoSorted = points(paretoIndices, :);
            
            % Sort by first coordinate
            [~, i] = sort(paretoSorted(:, 1));
            
            sorted = paretoSorted(i,:);
            lengths = diff([sorted(:, 1); ref(1)]);
            
            if size(points, 2) == 2
                V = sum(lengths .* (ref(2) - sorted(:, 2)));
            else
                V = 0;
                % March along line and integrate
                for i = 1:size(sorted, 1)
                    V_i = this.hyperVolumeMetric(sorted(1:i, 2:end), ref(2:end));
                    V = V + V_i * lengths(i);
                end
            end
            
        end
        
        function options_temp = setOptions(this, defaultoptions, options) %#ok<INUSL>
            % SETOPTIONS sets options and params before starting the iteration.
            %
            % INPUT:
            % - defaultoptions (struct)
            %   Default options of optimizer or default parameters of optimizer
            % - options (struct)
            %   Specific options of optimizer or specific parameters of optimizer
            %
            % OUPUT:
            % - options_temp (struct)
            %   Update properties (options/params)
            %
            % See also OPTIMIZE.
            
            options_temp = defaultoptions;
            options_names = fieldnames(options);
            for num_opt = 1:length(options_names)
                eval(['options_temp.', options_names{num_opt}, '= options.', options_names{num_opt}, ';']);
            end
            
        end
                
    end
    
    methods(Access = private)
        %% Private methods called by Optimizer
        
        function objectiveValues = evalObjective(this, objective, nObj, samples)
            % EVALOBJECTIVE evaluates objective function for samples.
            %
            % - objective (function handle)
            %   Objective function of optimization problem 
            % - nObj (integer scalar)
            %   Number of objectives
            % - samples (m by n matrix of doubles)
            %   Sampled design varibles, where m is number of samples 
            %   and n number of design variables per sample
            %
            % OUPUT:
            % - objectiveValues (m by n matrix of doubles)
            %   Objective functions values corresponding to samples, where
            %   m is number of evaluated samples and n is number of
            %   objectives
            %
            % See also OPTIMIZE.
            
            % Initialize
            objectiveValues = zeros(size(samples, 1), nObj);
            
            if this.options.numWorkers==1 % Single-core optimization
                
                for j = 1:size(samples, 1)
                    if nargin(objective) == 1
                        objectiveValues(j,:) = objective(samples(j,:));
                    else % If worker number is input
                        objectiveValues(j,:) = objective(samples(j,:), j); 
                    end
                end
                
            else % Multi-core optimization (REQUIRES Parallel Computing Toolbox)
                
                parfor j = 1:size(samples, 1)
                    if nargin(objective)==1
                        objectiveValues(j,:) = objective(samples(j,:));
                    else % If worker number is input
                        objectiveValues(j,:) = objective(samples(j,:), j); 
                    end
                end
                
            end
            
        end
        
        function outputStatus(this, numIters, numEvals, objectiveValue, objectiveChange, termFlag, time)
            % OUTPUTSTATUS writes optimization status to command line.
            %
            % INPUT:
            % - numIters (integer scalar)
            %   Number of iterations 
            % - numEvals (integer scalar)
            %   Number of objective function evaluations
            % - objectiveValue (double scalar)
            %   Best objective value or hypervolume in multi-objective case
            % - objectiveChange (double scalar)
            %   Change in best objective value (compared to previous
            %   iteration)
            % - termFlag (integer scalar)
            %   Termination condition
            % - time (double scalar)
            %   End of timer
            %
            % See also OPTIMIZE.
            
            % Choose output depending on termination flag
            switch termFlag
                case 0 % No termination
                    if numIters==0 % Initialization status
                        disp(['Starting optimization with ', class(this)])
                        disp('Iteration Objective-evals Objective-value Objective-change')
                        disp(['Init',...
                            repmat(char(32), 1, 16-length(num2str(numEvals))), num2str(numEvals),...
                            repmat(char(32), 1, 16-length(num2str(objectiveValue))), num2str(objectiveValue)])
                    else % Intermediate status
                        disp([repmat(char(32), 1, 4-length(num2str(numIters))), num2str(numIters),...
                            repmat(char(32), 1, 16-length(num2str(numEvals))), num2str(numEvals),...
                            repmat(char(32), 1, 16-length(num2str(objectiveValue))), num2str(objectiveValue),...
                            repmat(char(32), 1, 17-length(num2str(objectiveChange))), num2str(objectiveChange)])
                    end
                    
                case 1 % Termination due to maximum number of iterations without objective change
                    disp(['Maximum count of iterations without change of objective reached. Aborted optimization after ', num2str(numIters), ' iterations!'])
                    disp(['An objective function value of ', num2str((objectiveValue/this.options.rndFact)*this.options.rndFact), ' has been reached in ', num2str(time), ' seconds'])
                case 2 % Termination due to maximum number of iterations reached
                    disp(['Termination of optimization due to limit of ', num2str(this.options.maxIters), ' iterations!'])
                    disp(['An objective function value of ', num2str((objectiveValue/this.options.rndFact)*this.options.rndFact), ' has been reached in ', num2str(time), ' seconds'])
                case 3 % Termination due to maximum number of objective function evaluations reached
                    disp(['Termination of optimization after ', num2str(numIters), ' iterations due to limit of ', num2str(this.options.maxEvals), ' objective function evaluations!'])
                    disp(['An objective function value of ', num2str((objectiveValue/this.options.rndFact)*this.options.rndFact), ' has been reached in ', num2str(time), ' seconds'])
                case 4 % Termination due to termination by optimizer
                    disp(['Termination of optimization after ', num2str(numIters), ' iterations due to termination by algorithm.'])
                    disp(['An objective function value of ', num2str((objectiveValue/this.options.rndFact)*this.options.rndFact), ' has been reached in ', num2str(time), ' seconds.'])
            end
            
        end

    end
    
end