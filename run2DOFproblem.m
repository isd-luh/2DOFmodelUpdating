clear
close all
clc

% Add part of EngiO framework to path
% Reference:
% Berger et al. (2021), doi: 10.1016/j.advengsoft.2020.102959
addpath(genpath('EngiO'), '-frozen');

% Add problem folder to path
addpath(genpath('2DOF'), '-frozen');

%% Problem definition 2DOF system
% Reference:
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% String design variables
strDVs = {'k_1', 'k_2'};
nDVs = numel(strDVs);

% Lower and upper bounds
lowerBound = [0.01, 0.01];
upperBound = [2, 2];

% Choose distribution function (choices are "Gaussian" or "Weibull")
strInputDist = 'Gaussian'; 
if isequal(strInputDist, 'Gaussian')
    correct_k = [0.5, 1.5];
    meanEigenvalues = calcEigenvalues(correct_k);
    sigmaEigenvalues = [0.1, 0.05];
    scaleWB = [];
    shapeWB = [];
elseif isequal(strInputDist, 'Weibull')
    scaleWB = [3.4,0.22];
    shapeWB = [60,40];
    example_k = importdata(fullfile('2DOF','example_k_Weibull.mat'));
    exampleEigenvalues = calcEigenvalues(example_k);
    meanEigenvalues = mean(exampleEigenvalues);
    sigmaEigenvalues = std(exampleEigenvalues);
end

% Create folder for saving results
strResultFolder = ['ResultFolder2DOF_',strInputDist];
mkdir(strResultFolder)


%% Bayesian model updating (BMU) using the transitional Markov chain Monte Carlo (TMCMC) sampling technique
fprintf('\nBMU calculation using TMCMC\n')

% Set number of samples for TMCMC
nSamplesTMCMC = 500;
save(fullfile(strResultFolder,'nSamplesTMCMC'), 'nSamplesTMCMC')

nRunsTMCMC = 3;
fprintf('\n  - Calculate %d runs\n', nRunsTMCMC);
optimalDVsBMU = {};
optimalEVsBMU = {};
for iRunTMCMC = 1:nRunsTMCMC
    fprintf('\n  - Run No. %d \n', iRunTMCMC);

    % Set model updating problem
    problemTMCMC = class2DOF(meanEigenvalues, sigmaEigenvalues, [], 'BMU', []);
    LogL = @(x, index)problemTMCMC.calcObjValue(x, index);
    
    % Run TMCMC
    [DVsamplesBMU, ~] = calcTMCMC(LogL, lowerBound, upperBound, nSamplesTMCMC, 0);
    optimalDVsBMU{1,iRunTMCMC} = DVsamplesBMU{1,end};
    optimalEVsBMU{1,iRunTMCMC} = calcEigenvalues(optimalDVsBMU{1,iRunTMCMC});

end
save(fullfile(strResultFolder, 'optimalDVsBMU'), 'optimalDVsBMU')
save(fullfile(strResultFolder, 'optimalEVsBMU'), 'optimalEVsBMU')


%% Repeated deterministic model updating (RDMU)

% Sampling method 
% -> choices are Halton sequence (HS), Sobol sequence (SS), Monte Carlo (MC)
strSamplerRDMU = 'HS'; 
% Optimization algorithm 
% -> choices are global pattern search (GPS), genetic algorithm (GA), evolution strategy (ES)
strOptimizerRDMU = 'GPS';
fprintf('\nRDMU calculation using %s%s\n', strSamplerRDMU, strOptimizerRDMU);

fprintf('\n  - Generate samples\n');

% Set number of samples for sample generation (SG)
nSamplesRDMU = 500;

% Generate eigenvalue samples (dependent on choice of sampling method)
EVsamplesRDMU = generateEVsamples(strSamplerRDMU, nDVs, nSamplesRDMU, ...
    strInputDist, meanEigenvalues, sigmaEigenvalues, scaleWB, shapeWB);
save(fullfile(strResultFolder,'EVsamplesRDMU'), 'EVsamplesRDMU')


fprintf('\n  - Updating\n');
% Settings for optimization runs (dependent on choice of optimization algorithm)
[optimizerRDMU, optParamsRDMU] = optSettings(strOptimizerRDMU);
optOptionsRDMU = struct('maxEvals', 1000, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);

% Model updating for each sample
optimalDVsRDMU = [];
parfor iSampleRDMU = 1:nSamplesRDMU 
    currentEVsampleRDMU = EVsamplesRDMU(iSampleRDMU, :);

    % Set model updating problem
    problemRDMU = class2DOF(currentEVsampleRDMU, [], [], 'RDMU', []);
    objFunRDMU = @(x, index)problemRDMU.calcObjValue(x, index);

    % Run optimization
    [optimalDVsRDMU(iSampleRDMU, :), ~, ~, ~, ~] = optimizerRDMU.optimize(objFunRDMU, [], lowerBound, upperBound, ...
        optOptionsRDMU, optParamsRDMU);
    
end
optimalEVsRDMU = calcEigenvalues(optimalDVsRDMU);
save(fullfile(strResultFolder, 'optimalDVsRDMU'), 'optimalDVsRDMU')
save(fullfile(strResultFolder, 'optimalEVsRDMU'), 'optimalEVsRDMU')


%% Repeated deterministic meta-model updating (metaRDMU)
% Sampling method 
% -> choices are Halton sequence (HS), Sobol sequence (SS), Monte Carlo (MC)
strSamplerMetaRDMU = 'HS'; 
% Optimization algorithm 
% -> choices are global pattern search (GPS), genetic algorithm (GA), evolution strategy (ES)
strOptimizerMetaRDMU = 'GPS';
fprintf('\nMetaRDMU calculation using meta%s%s\n', strSamplerMetaRDMU, strOptimizerMetaRDMU);

fprintf('\n  - Set up meta-model\n');
% Settings for optimization run building the meta-model 
% -> here using the global pattern search optimization algorithm
optimizerMetaModel = GlobalPattern();
optParamsMetaModel  = struct('nTrack' , 50);
optOptionsMetaModel = struct('maxEvals', 250, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);

% Create data queue
queueMetaModel  = parallel.pool.PollableDataQueue;

% Set model updating problem based on mean eigenvalues
problemMetaModel = class2DOF(meanEigenvalues, [], queueMetaModel, 'RDMU', []);
objFunMetaModel = @(x, index)problemMetaModel.calcObjValue(x, index);

% Run optimization based on mean eigenvalues
[~, ~, ~, ~, ~] = optimizerMetaModel.optimize(objFunMetaModel, [], lowerBound, upperBound, ...
    optOptionsMetaModel, optParamsMetaModel);

% Get results from queue
queueResultsMetaModel  = [];
while true
    [dataMetaModel, bDataMetaModel] = poll(queueMetaModel);
    queueResultsMetaModel = [queueResultsMetaModel; dataMetaModel];
    if ~bDataMetaModel
        break
    end
end
DVsamplesMetaModel = queueResultsMetaModel(:,1:2);
EVsamplesMetaModel = queueResultsMetaModel(:,3:end);

fprintf('\n  - Build meta-model\n');
outputMetaModel = scatterModel(DVsamplesMetaModel, EVsamplesMetaModel, lowerBound, upperBound);

fprintf('\n  - Generate samples\n');
% Set number of samples for sample generation (SG) 
% -> using the Halton sequence sampling method
nSamplesMetaRDMU = 500;

% Generate eigenvalue samples 
EVsamplesMetaRDMU = generateEVsamples(strSamplerMetaRDMU, nDVs, nSamplesMetaRDMU, ...
    strInputDist, meanEigenvalues, sigmaEigenvalues, scaleWB, shapeWB);
save(fullfile(strResultFolder,'EVsamplesMetaRDMU'), 'EVsamplesMetaRDMU')


% 4) Updating on meta-model
fprintf('\n  - Updating on meta-model\n');

% Settings for optimization run on meta-model
[optimizerMetaRDMU, optParamsMetaRDMU] = optSettings(strOptimizerMetaRDMU);
optOptionsMetaRDMU = struct('maxEvals', 1000, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);

% Model updating for each sample on meta-model
optimalDVsMetaRDMU = [];
for iSampleMetaRDMU = 1:nSamplesMetaRDMU
    currentEVsampleMetaRDMU = EVsamplesMetaRDMU(iSampleMetaRDMU, :);

    % Set model updating problem
    problemMetaRDMU = class2DOF(currentEVsampleMetaRDMU, [], [], 'MetaRDMU', @outputMetaModel.predict);
    objFunMetaRDMU = @(x, index)problemMetaRDMU.calcObjValue(x, index);

    % Run optimization
    [optimalDVsMetaRDMU(iSampleMetaRDMU, :), ~, ~, ~, ~] = ...
        optimizerMetaRDMU.optimize(objFunMetaRDMU, [], lowerBound, upperBound, ...
        optOptionsMetaRDMU, optParamsMetaRDMU);
    
end
optimalEVsMetaRDMU = calcEigenvalues(optimalDVsMetaRDMU);
save(fullfile(strResultFolder, 'optimalDVsMetaRDMU'), 'optimalDVsMetaRDMU')
save(fullfile(strResultFolder, 'optimalEVsMetaRDMU'), 'optimalEVsMetaRDMU')



%% Plot results
close all
strMethods = {'BMU (TMCMC)', ...
    ['RDMU (', strSamplerRDMU, strOptimizerRDMU, ')'], ...
    ['MetaRDMU (', strSamplerMetaRDMU, strOptimizerMetaRDMU, ')']};

% Plot colors
blueColor = [0 0.4470 0.7410];
redColor = [0.8500 0.3250 0.0980];
yellowColor = [0.9290 0.6940 0.1250];


% Plot comparison samples
for iDV = 1:nDVs
    figure
    hold on
    
    plot1 = cdfplot(optimalDVsBMU{1,1}(:, iDV));    % plot first TMCMC run
    set(plot1, 'LineWidth', 2, 'Color', blueColor)

    plot2 = cdfplot(optimalDVsRDMU(:, iDV));
    set(plot2, 'LineWidth', 2.5, 'Color', redColor)
    
    plot3 = cdfplot(optimalDVsMetaRDMU(:, iDV));
    set(plot3, 'LineWidth', 2.5, 'LineStyle', '--', 'Color', yellowColor)
    
    if isequal(strInputDist,'Gaussian')
        xline(correct_k(iDV), '--k', 'LineWidth', 2)
    elseif isequal(strInputDist,'Weibull')
        plot4 = cdfplot(example_k(:, iDV));
        set(plot4, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
    end

    for iRunTMCMC = 2:nRunsTMCMC    % plot other TMCMC run
        plot1 = cdfplot(optimalDVsBMU{1,iRunTMCMC}(:, iDV));
        set(plot1, 'LineWidth', 2, 'Color', blueColor)
    end


    xlabel(['{\it k}_', num2str(iDV)])
    ylabel(['F({\itk}_', num2str(iDV), ')'])
    xlim([lowerBound(iDV), upperBound(iDV)])
    title('')
    if iDV == 1
        xlim([0.1,0.9])
    elseif iDV == 2
        legend([strMethods, 'Correct solution'])
        xlim([1.3,1.7])
    end
    ax = gca;
    ax.FontSize = 12;
    ax.FontName = 'Times';
    grid on
    box on

    set(gcf, 'Name', 'Comparison cdfs')
    savefig(fullfile(strResultFolder,['fig2DOF_ComparisonCDFs_k', num2str(iDV)]))
end


% Plot samples
nPlotSamples = 200;
figure
hold on

plot(optimalDVsBMU{1,1}(1:nPlotSamples,1), optimalDVsBMU{1,1}(1:nPlotSamples,2), 'd', ...
    'LineWidth', 1.3, 'Color', blueColor)       % plot first TMCMC run
plot(optimalDVsRDMU(1:nPlotSamples,1), optimalDVsRDMU(1:nPlotSamples,2), 's', ...
    'LineWidth', 1.3, 'Color', redColor, 'MarkerSize', 7)
plot(optimalDVsMetaRDMU(1:nPlotSamples,1), optimalDVsMetaRDMU(1:nPlotSamples,2), 'o', ...
    'LineWidth', 1.3, 'Color', yellowColor, 'MarkerSize', 9)
if isequal(strInputDist,'Gaussian')
    plot(correct_k(1), correct_k(2), 'k +', 'MarkerSize', 10, 'LineWidth', 2)
    legend([strMethods, 'Correct solution'])
elseif isequal(strInputDist,'Weibull')
    legend(strMethods)
end

xlabel('{\it k}_1')
ylabel('{\it k}_2')
xlim([lowerBound(1), upperBound(1)])
ylim([lowerBound(2), upperBound(2)])
title('')
grid on
box on
xlim([0.1,0.9])
ylim([1.3,1.7])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Name', 'Comparison samples')
savefig(fullfile(strResultFolder,'fig2DOF_DVspace'))


% Plot eigenvalues
nPlotSamples = 200;
figure
hold on

plot(optimalEVsBMU{1,1}(1:nPlotSamples,1), optimalEVsBMU{1,1}(1:nPlotSamples,2), 'd', ...
    'LineWidth', 1.3, 'Color', blueColor)       % plot first TMCMC run
plot(optimalEVsRDMU(1:nPlotSamples,1), optimalEVsRDMU(1:nPlotSamples,2), 's', ...
    'LineWidth', 1.3, 'Color', redColor, 'MarkerSize', 7)
plot(optimalEVsMetaRDMU(1:nPlotSamples,1), optimalEVsMetaRDMU(1:nPlotSamples,2), 'o', ...
    'LineWidth', 1.3, 'Color', yellowColor, 'MarkerSize', 9)
plot(meanEigenvalues(1), meanEigenvalues(2), 'k +', 'MarkerSize', 10, 'LineWidth', 2)

t = -pi:0.01:pi;
x1 = meanEigenvalues(1) + sigmaEigenvalues(1) * cos(t);
y1 = meanEigenvalues(2) + sigmaEigenvalues(2) * sin(t);
x2 = meanEigenvalues(1) + 2*sigmaEigenvalues(1) * cos(t);
y2 = meanEigenvalues(2) + 2*sigmaEigenvalues(2) * sin(t);
plot(x1, y1, 'k', 'LineWidth', 1.3)
plot(x2, y2, '--k', 'LineWidth', 1.3)

xlabel('\lambda_1')
ylabel('\lambda_2')
xlim([2.9,3.7])
ylim([0.05,0.4])
title('')
grid on
box on
legend([strMethods, 'Mean eigenvalue', '1\bf{\sigma}_\lambda', '2\bf{\sigma}_\lambda'])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';

set(gcf, 'Name', 'Comparison eigenvalues')
savefig(fullfile(strResultFolder,'fig2DOF_EVspace'))

