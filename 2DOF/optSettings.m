function [optimizerRDMU, optParamsRDMU] = optSettings(strOptimizer)

% Settings for optimization runs (dependent on choice of optimization algorithm)
if isequal(strOptimizer,'GPS')
    optimizerRDMU = GlobalPattern();
    optParamsRDMU = struct('nTrack', 20);
elseif isequal(strOptimizer,'GA')
    optimizerRDMU = Genetic();
    optParamsRDMU = struct('popSize', 50);
elseif isequal(strOptimizer,'ES')
    optimizerRDMU = EvolutionStrategy();
    optParamsRDMU = struct('mu', 5, 'maxLoops', 250);
end

end