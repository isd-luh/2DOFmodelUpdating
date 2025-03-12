% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef Genetic < Optimizer
    % GENETIC is a single-objective Genetic Algorithm derived from Optimizer.
    %
    % LITERATURE:
    %   D.E. Goldberg (1989) Genetic Algorithms in Search,
    %   Optimization, and Machine Learning, Addison-Wesley
    
    methods
        
        function this = Genetic()
            
            this.defaultparams = struct('crossPercentage', 0.8,...     % Crossover percentage
                'crossRange', 0.4,...          % Crossover range factor
                'mutationRange', 0.2,...       % Mutation range factor
                'popSize', 50,...              % Population size
                'selectionType', 'default',... % Valid selection types: 'Random', 'RouletteWheel' (default), 'Tournament'
                'tournamentSize', 3);          % Tournament size (only used in case of Tournament Selection)
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % Initialize all states
            this.states = struct('f', inf(this.params.popSize, 1),... % Objective function values for current samples
                'f_g', inf,... % Objective function values for global best sample
                'f_temp', inf(this.params.popSize, 1),... % Objective function values from last iteration loop
                'numLoops', 0,... % Current generation
                'x', zeros(this.params.popSize, length(x_lb)),... % Current samples
                'x_c', zeros(floor(this.params.popSize*this.params.crossPercentage), length(x_lb)),... % Crossover samples
                'x_g', 0,... % Global best sample
                'x_m', zeros(this.params.popSize-ceil(this.params.popSize*this.params.crossPercentage), length(x_lb)),... % Mutated samples
                'x_temp', zeros(this.params.popSize, length(x_lb))); % Samples from last iteration loop
            
            % Initialize particle positions
            for j = 1:this.params.popSize
                if isempty(vartype)
                    this.states.x(j,:) = x_lb+(x_ub-x_lb).*rand(size(this.states.x(j,:)));
                else
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'real')
                            this.states.x(j, m) = x_lb(m)+(x_ub(m)-x_lb(m)).*rand(1);
                        elseif strcmpi(vartype{m}, 'int')
                            this.states.x(j, m) = randi([x_lb(m) x_ub(m)]);
                        end
                    end
                end
            end
            
            % Initialize roulette wheel selection array
            this.states.indArray = [];
            
            for k = 1:size(this.states.x, 1)
                this.states.indArray = [this.states.indArray; repmat(k, size(this.states.x, 1)-k+1, 1)]; %#ok<AGROW>
            end
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % Increment loop counter
            this.states.numLoops = this.states.numLoops+1;
            if this.states.numLoops == 1
                samples = this.states.x;
                return
            end
            
            % Allocate states
            this.states.x_temp = this.states.x;
            this.states.f_temp = this.states.f;
            this.states.x_c = [];
            this.states.x_m = [];
            
            % Loop over crossover particles
            for j = 1:max(1, floor(this.params.popSize*this.params.crossPercentage/2))
                
                % Select parents
                if strcmpi(this.params.selectionType, 'Random')
                    
                    ind1 = randi([1 this.params.popSize]);
                    ind2 = randi([1 this.params.popSize]);
                    
                elseif strcmpi(this.params.selectionType, 'default') || strcmpi('this.params.selectionType', 'RouletteWheel')
                    ind1 = this.states.indArray(randi([1 length(this.states.indArray)]));
                    ind2 = this.states.indArray(randi([1 length(this.states.indArray)]));
                    
                elseif strcmpi(this.params.selectionType, 'Tournament')
                    
                    indDuelist1 = randi([1 this.params.popSize]);
                    indDuelist2 = randi([1 this.params.popSize]);
                    
                    for numTournaments = 1:this.params.tournamentSize
                        
                        if this.states.f(indDuelist1)>this.states.f(indDuelist2)
                            indDuelist1 = randi([1 this.params.popSize]);
                            ind1 = indDuelist2;
                        else
                            indDuelist2 = randi([1 this.params.popSize]);
                            ind1 = indDuelist1;
                        end
                        
                    end
                    
                    indDuelist1 = randi([1 this.params.popSize]);
                    indDuelist2 = randi([1 this.params.popSize]);
                    
                    for numTournaments = 1:this.params.tournamentSize
                        
                        if this.states.f(indDuelist1)>this.states.f(indDuelist2)
                            indDuelist1 = randi([1 this.params.popSize]);
                            ind2 = indDuelist2;
                        else
                            indDuelist2 = randi([1 this.params.popSize]);
                            ind2 = indDuelist1;
                        end
                        
                    end
                    
                else
                    
                    error('Invalid selectionType provided.')
                    
                end
                
                parent1 = this.states.x(ind1,:); % First parent
                parent2 = this.states.x(ind2,:); % Second parent
                
                alpha = rand(size(parent1)) * (1+2*this.params.crossRange) - this.params.crossRange; % Index alpha
                
                offspring1 = alpha.*parent1+(1-alpha).*parent2; % First offspring position
                offspring2 = alpha.*parent2+(1-alpha).*parent1; % Seond offspring position
                
                % Saturate offsprings at boundaries
                offspring1 = min(max(offspring1, x_lb), x_ub);
                offspring2 = min(max(offspring2, x_lb), x_ub);
                
                % Append offspring positions to crossover position state matrix
                this.states.x_c = [this.states.x_c;offspring1;offspring2];
                
            end
            
            % Loop over mutation particles
            for j = 1:this.params.popSize-ceil(this.params.popSize*this.params.crossPercentage)
                
                % Select random individuum
                xrand = this.states.x(randi([1 this.params.popSize]),:);
                
                % Mutate it
                mutant = xrand+this.params.mutationRange*(x_ub-x_lb).*randn(size(xrand));
                
                % Saturate mutant
                mutant = min(max(mutant, x_lb), x_ub);
                
                this.states.x_m = [this.states.x_m;mutant];
            end
            
            this.states.x = [this.states.x_c;this.states.x_m];
            
            % Integer variables are rounded
            if ~isempty(vartype)
                for j = 1:this.params.popSize
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'int')
                            this.states.x(j, m) = round(this.states.x(j, m));
                        end
                    end
                end
            end
            
            % Store particle positions as samples
            samples = this.states.x;
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Append values from last generations to current one
            this.states.f = [this.states.f_temp;objectiveValues];
            this.states.x = [this.states.x_temp;samples];
            
            % Sort Lagrangian values
            [this.states.f, sortOrder]= sort(this.states.f);
            this.states.x = this.states.x(sortOrder,:);
            
            % Truncate vectors
            this.states.f = this.states.f(1:this.params.popSize);
            this.states.x = this.states.x(1:this.params.popSize,:);
            
            % Loop over all particles
            for j = 1:this.params.popSize
                
                % Check whether best solution was found
                if this.states.f(j)<this.states.f_g
                    this.states.x_g = this.states.x(j,:);
                    this.states.f_g = this.states.f(j);
                end
                
            end
            
            terminate = false;
            
        end
        
    end
    
end