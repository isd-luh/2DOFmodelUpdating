% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef EvolutionStrategy < Optimizer
    % EVOLUTIONSTRATEGY is a single-objective global Algorithm derived from
    % Optimizer.
    %
    % LITERATURE:
    %   H.-P. Schwefel (1981) Numerical optimization of computer models
    %   H.-G. Beyer and H.-P. Schwefel (2002) Evolution strategies - A
    %   Comprehensive Introduction, Natural Compution 1, 3-52
    
    methods
        
        function this = EvolutionStrategy()
            
            this.defaultparams = struct('maxLoops', 250,... % maximum total count of generations
                'mu', 5,... % number of parental samples per generation
                'rho', 2,... % number of mating participants (rho <= mu)
                'lambda', 10,... % number of offspring samples per generation (lambda >= mu)
                'kappa', 50,... % maximum generation number for survival of samples
                'xi', 10,... % participants in tournament selection (xi <= lambda)
                'mutlim', 0.2,... % limit for mutation factor (default = 0.2 --> Rechenberg's 1/5 success rule)
                'mutfactor', 0.85,... % mutation factor for changes (default = 0.85, 'mutfactor' < 1)
                'sigma_scale', 0.1,... % scale factor for mutation factor sigma
                'sigma_maxscale', 0.1,... % maximum scale factor for mutation factor sigma
                'rho_radius', 999,... % visibility of mating partners (to do...)
                'recombination', 'intermediary',... % RECOMBINATION ('intermediary'|'dominant'|'fitprop')
                'mutation', 'individual',... % MUTATION ('uniform'|'individual')
                'selection', 'plus',... % SELECTION ('comma'|'plus'|'kappa'|'tournament'|'fitprop')
                'paracontrol', 'mutationrate'); % PARAMETER CONTROL ('no'|'deterministic'|'mutationrate'|'selfadaptive')
            
            % RECOMBINATION -------------------
            % r1: 'intermediary' > intermediary recombination for each dimension
            % r2: 'dominant' > recombination for each dimension as random dominant
            % r3: 'fitprop' > fitness-proportional selection for mating partners (to do...)
            
            % GAUSSIAN MUTATION ---------------
            % m1: 'uniform' > Gaussian mutation with same value for all dimensions
            % m2: 'individual' > Gaussian mutation with individual values
            
            % SELECTION -----------------------
            % s1: 'comma' > only from offsprings
            % s2: 'plus' > from direct parents + offsprings
            % s3: 'kappa' > from all parental generations (with max. generation age of 'kappa') + offsprings
            % s4: 'tournament' > random tournaments of 'xi' competitors (to do...)
            % s5: 'fitprop' > fitness-proportional selection (to do...)
            
            % PARAMETER CONTROL ---------------
            % c1: 'constant' > constant stepsize of mutation parameter
            % c2: 'deterministic' > deterministic control of mutation stepsize
            % c3: 'mutationrate' > mutation rate control by RECHENBERG
            % c4: 'selfadaptive' > self-adaption of mutation stepsize by KRAMER (to do...)
        end
        
        function initialize(this, nObj,~, x_lb, x_ub)
            
            % Determination of problem dimension
            dim = numel(x_lb);
            
            % Initialize all states
            this.states = struct('x_mu', zeros(this.params.mu, dim),... % initial samples / samples of parental generaltion
                'f_mu', inf(this.params.mu, 1),... % parental fitness
                'x_lambda', zeros(this.params.lambda, dim),... % offspring samples
                'f_lambda', zeros(this.params.lambda, 1),... % offspring fitness
                'x_history', [],... % history of all samples (design parameters)
                'f_history', [],... % history of all fitness values (objective parameters)
                'sigma_ini', this.params.sigma_scale*(max(x_ub)-min(x_lb)),... % initial stepsize of mutation parameter
                'sigma', this.params.sigma_scale*(max(x_ub)-min(x_lb)),... % stepsize of mutation parameter
                'sigma_max', this.params.sigma_maxscale*(max(x_ub)-min(x_lb)),... % upper bound of mutation factor sigma
                'gen_no', 0); % counter for number of generations
            
            % Initialize sampling positions
            for i = 1:this.params.mu
                for j = 1:dim
                    this.states.x_mu(i, j) = x_lb(j)+(x_ub(j)-x_lb(j))*rand;
                end
            end
            
        end
        
        function samples = generateSamples(this, nObj,~, x_lb, x_ub)
            
            % Determination of problem dimension
            dim = numel(x_lb);
            
            % Recombination + mutation
            for i = 1:this.params.lambda % lambda offsprings
                
                % RECOMBINATION
                
                
                % Draw randomly 'rho' individuals
                [~, ind] = sort(randperm(this.params.mu));
                x_mu_mixed = this.states.x_mu(ind,:);
                x_mating = x_mu_mixed(1:this.params.rho,:);
                
                % 'intermediary' recombination (r1)
                if isequal(this.params.recombination, 'intermediary')
                    x_recomb = mean(x_mating);
                    
                    % 'dominant' recombination (r2)
                elseif isequal(this.params.recombination, 'dominant')
                    x_recomb = zeros(1, dim);
                    for j = 1:dim
                        x_recomb(:, j) = x_mating(randi(this.params.rho, 1, 1), j);
                    end
                    
                    % 'fitprop' recombination (r3)
                elseif isequal(this.params.recombination, 'fitprop')
                    % ... (to do)
                    
                    % Error interception
                else
                    error('select = ''intermediary'' or ''dominant'' or ''fitprop''!')
                end
                
                % GAUSSIAN MUTATION
                
                % With 'uniform' parameter (m1)
                if isequal(this.params.mutation, 'uniform')
                    mutation_value = this.states.sigma * randn;
                    mutation_vector = repmat(mutation_value, 1, dim);
                    
                    % With 'individual' parameters (m2)
                elseif isequal(this.params.mutation, 'individual')
                    mutation_vector = this.states.sigma * randn(1, dim);
                    
                    % Error interception
                else
                    error('select = ''uniform'' or ''individual''!')
                end
                
                % Mutate recombined individuals
                this.states.x_lambda(i,:) = x_recomb + mutation_vector;
                
                % Check if staying inside search area (correction if required)
                for j = 1:dim
                    if this.states.x_lambda(i, j)<x_lb(j)
                        this.states.x_lambda(i, j) = x_lb(j);
                    elseif this.states.x_lambda(i, j)>x_ub(j)
                        this.states.x_lambda(i, j) = x_ub(j);
                    end
                end
                
            end
            
            % Schedule of sampling points (design variables)
            samples = this.states.x_lambda;
            
            % History of all samples ('lambda' offsprings over all generations)
            this.states.x_history = [this.states.x_history
                this.states.x_lambda];
        end
        
        function terminate = processResults(this,~, objectiveValues)
            
            % Evaluation of fitness values for offspring individuals (new generation)
            this.states.f_lambda = objectiveValues;
            
            % Increase generation counter
            this.states.gen_no = this.states.gen_no+1;
            
            % History of all fitness values
            this.states.f_history = [this.states.f_history
                this.states.f_lambda];
            
            % SELECTION
            
            % 'comma'-selection (s1)
            if isequal(this.params.selection, 'comma')
                [f_sort, ind] = sort(this.states.f_lambda);
                x_sort = this.states.x_lambda(ind,:);
                this.states.x_mu = x_sort(1:this.params.mu,:);
                this.states.f_mu = f_sort(1:this.params.mu,:);
                
                % 'plus'-selection (s2)
            elseif isequal(this.params.selection, 'plus')
                p_plus_o_samples = [this.states.x_mu ; this.states.x_lambda];
                p_plus_o_fitness = [this.states.f_mu ; this.states.f_lambda];
                [f_sort, ind] = sort(p_plus_o_fitness);
                x_sort = p_plus_o_samples(ind,:);
                this.states.x_mu = x_sort(1:this.params.mu,:);
                this.states.f_mu = f_sort(1:this.params.mu,:);
                
                % 'kappa'-selection (s3)
            elseif isequal(this.params.selection, 'kappa')
                rememberN = this.params.kappa*this.params.lambda;
                if size(this.states.x_history, 1) < rememberN
                    p_plus_o_samples = [this.states.x_history ; this.states.x_lambda];
                    p_plus_o_fitness = [this.states.f_history ; this.states.f_lambda];
                else
                    p_plus_o_samples = [this.states.x_history(end-rememberN+1:end,:) ; this.states.x_lambda];
                    p_plus_o_fitness = [this.states.f_history(end-rememberN+1:end,:) ; this.states.f_lambda];
                end
                [f_sort, ind] = sort(p_plus_o_fitness);
                x_sort = p_plus_o_samples(ind,:);
                this.states.x_mu = x_sort(1:this.params.mu,:);
                this.states.f_mu = f_sort(1:this.params.mu,:);
                
                % 'tournament'-selection (s4)
            elseif isequal(this.params.selection, 'tournament')
                % ... (to do)
                
                % 'fitprop'-selection (s5)
            elseif isequal(this.params.selection, 'fitprop')
                % ... (to do)
                
                % Error interception
            else
                error('select = ''comma'' or ''plus'' or ''kappa'' or ''tournament'' or ''fitprop''!')
            end
            
            % CONTROL OF MUTATION PARAMETER
            
            % 'constant' > constant stepsize of mutation parameter (c1)
            if isequal(this.params.paracontrol, 'constant')
                this.states.sigma = this.states.sigma_ini;
                
                % 'deterministic' > deterministic parameter control (c2)
            elseif isequal(this.params.paracontrol, 'deterministic')
                this.states.sigma = this.states.sigma_ini*(1-0.9*this.states.gen_no/this.params.maxLoops);
                
                % 'mutationrate' > mutation rate control by RECHENBERG (c3)
            elseif isequal(this.params.paracontrol, 'mutationrate')
                % Fitness of 10 preceding generations
                % ---------------------------------------------------------
                % "After every n mutations, check how many successes have
                % Occured over the preceding 10*n mutations. If this number
                % Is less than 2*n, multiply the step lengths by the factor
                % 0.85; divide them by 0.85 if more than 2*n successes
                % Occured." (formulation by SCHWEFEL, 1995)
                % ---------------------------------------------------------
                anc_gen = this.params.kappa;
                if size(this.states.f_history, 1) >= anc_gen*this.params.lambda
                    % Adjustment of mutation parameter sigma according to 'success'
                    f_pgen = this.states.f_history(end-anc_gen*this.params.lambda+1:end,:);
                    success = sum(f_pgen <= min(this.states.f_lambda))/size(f_pgen, 1);
                    if success < this.params.mutlim
                        this.states.sigma = this.states.sigma*this.params.mutfactor;
                    elseif success > this.params.mutlim
                        this.states.sigma = min(this.states.sigma/this.params.mutfactor, this.states.sigma_max);
                    else % success == this.params.mutlim
                        % 'this.states.sigma' remains the same
                    end
                else % to little data basis --> constant stepsize of mutation parameter
                    this.states.sigma = this.states.sigma_ini;
                end
                
                % 'selfadaptive' > self-adaption by KRAMER (c4)
            elseif isequal(this.params.paracontrol, 'selfadaptive')
                % ... (to do)
                
                % Error interception
            else
                error('select = ''constant'' or ''deterministic'' or ''mutationrate'' or ''selfadaptive''!')
            end
            
            terminate = false;
            
        end
        
    end
    
end
