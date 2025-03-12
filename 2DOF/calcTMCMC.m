function [samplesTMCMC, resultsTMCMC] = calcTMCMC(LogL, lowerBound, upperBound, nSamplesTMCMC, pj_0)
% This function calculates the transitional Markov chain Monte Carlo (TMCMC) method 
% 
% References: 
% Ching and Chen (2007), doi: 10.1061/(ASCE)0733-9399(2007)133:7(816)
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% INPUT: 
% - Problem function handle: LogL
% - Lower and upper bounds: lowerBound, upperBound
% - Number of samples utilized: nSamplesTMCMC
% - Initial tempering parameter: pj_0

% OUTPUT: 
% - Generated samples: samplesTMCMC
% - Calculated objective function values: resultsTMCMC

%%
% Generate Monte Carlo samples
nDim = size(lowerBound,2);
iIteration = 1;
pj(iIteration) = pj_0;
samples_0 = zeros(nSamplesTMCMC, nDim);
for iDim = 1:nDim
    samples_0(:,iDim) = unifrnd(lowerBound(iDim), upperBound(iDim), nSamplesTMCMC, 1);
end

% Calculate first distribution in whole parameter space
parfor iSample = 1:nSamplesTMCMC 
    currentResultsLogL = LogL(samples_0(iSample,:), iSample);
    % "Death Penalty"
    if isequal(currentResultsLogL{1}, Inf)
        currentResultsLogL{1} = -realmax;
    end
    results_0(iSample) = currentResultsLogL{1};
end
samplesTMCMC{iIteration} = samples_0;
resultsTMCMC{iIteration} = results_0;
fprintf('\t \t Iteration:\t%d\t\t pji:\t%f\n', 0, 0);

iCount = 0;
while pj(iIteration) < 1
    iCount = iCount + 1;

    % Evaluation stepsize
    wj = @(e) exp(abs(e) * resultsTMCMC{iIteration});
    fmin = @(e) std(wj(e)) - 1 * mean(wj(e)) + realmin;
    e = abs(fzero(fmin, 0));
    pji = min(1, e + pj(iIteration));
    fprintf('\t \t Iteration:\t%d\t\t pji:\t%f\n', iIteration, pji);

    mu = zeros(1, nDim);
    a = (pji - pj(iIteration)) * (resultsTMCMC{iIteration});
    wji = exp(a);
    wj_norm = wji./sum(wji);
    for iSample = 1:nSamplesTMCMC
        mu = mu + wj_norm(iSample) * samplesTMCMC{iIteration}(iSample,:);
    end

    % Calculate covariance matrix
    covGauss = zeros(nDim);
    for iSample = 1:nSamplesTMCMC
        tk_mu = samplesTMCMC{iIteration}(iSample,:) - mu;
        covGauss = covGauss + wj_norm(iSample)*(tk_mu'*tk_mu);
    end
    beta = 0.2; 
    covGauss = beta^2 * covGauss;

    % Sample generation
    samples_i_index = randsample(nSamplesTMCMC, nSamplesTMCMC, true, wj_norm);
    samples_i = samplesTMCMC{iIteration}(samples_i_index,:);
    results_i = resultsTMCMC{iIteration}(samples_i_index);

    % MCMC with one step for new samples
    parfor iSample = 1:nSamplesTMCMC 
        thetaLead = samples_i(iSample,:);
        logLLead = results_i(iSample);

        % Candidate sample generation
        while true
            thetaCand = mvnrnd(thetaLead, covGauss);
            if all([(thetaCand > lowerBound),  (thetaCand < upperBound)])
                break;
            end
        end

        % Sample calculation
        currentCand = LogL(thetaCand, iSample);
        % "Death Penalty"
        if isequal(currentCand{1}, Inf)
            currentCand{1} = -realmax; 
        end
        logLCand = currentCand{1};

        % Acceptance / rejection step
        alpha = exp((logLCand - logLLead));
        if rand <= min(1, alpha)
            theta_j1(iSample, :) = thetaCand;
            logL_j1(iSample) = logLCand;
            thetaLead = thetaCand;
            logLLead = logLCand;
        else
            theta_j1(iSample, :) = thetaLead;
            logL_j1(iSample) = logLLead;
        end
    end
    iIteration = iIteration + 1;
    pj(iIteration) = pji;
    samplesTMCMC{iIteration} = theta_j1;
    resultsTMCMC{iIteration} = logL_j1;
end

end