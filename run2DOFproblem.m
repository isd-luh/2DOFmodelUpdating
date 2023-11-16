clear
close all
clc

% Add part of EngiO framework to path
% Reference:
% Berger et al. (2021), doi: 10.1016/j.advengsoft.2020.102959
addpath(genpath(fullfile('EngiO')), '-frozen');

% Add problem folder to path
addpath(genpath(fullfile('2DOF')), '-frozen');

% Create folder for saving results
strResultFolder = 'ResultFolder2DOF';
mkdir(strResultFolder)

%% Problem definition 2DOF system
% Reference:
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% String design variables
cstrDVs = {'k_1', 'k_2'};

% Correct design variables (= spring stiffnesses)
afDVsCorrect = [0.5, 1.5;...
    3,   0.25];

% Standard distribution noise
afSigmaNoise = [0.1, 0.05];

% Calculate correct eigenvalues
afEigenvaluesCorrect = calcEigenvalues(afDVsCorrect(1,:));
nDVs = size(afDVsCorrect,2);

% Lower and upper bounds
afLB = [0.01, 0.01];
afUB = [4,    4   ];

% Run model updating with implemented methods 'TMCMC', 'MCGO' and 'meta-MCGO'
%% Transitional Markov chain Monte Carlo (TMCMC)
fprintf('\nCalculation using TMCMC\n')

% Create data queue
mfQueueTMCMC = parallel.pool.PollableDataQueue;

% Set model updating problem
problemTMCMC = class2DOF(afEigenvaluesCorrect, afSigmaNoise, mfQueueTMCMC, 'TMCMC', []);
probTMCMC = @(x, index)problemTMCMC.calcObjValue(x, index);

% Set number of Monte Carlo samples
nMCSamplesTMCMC = 1000;

% Run TMCMC
[cSamplesTMCMC, ~] = calcTMCMC(probTMCMC, afLB, afUB, nMCSamplesTMCMC, 0);
mOptimalDVsTMCMC = cSamplesTMCMC{1,end};
mfEigenvaluesOptimalTMCMC = calcEigenvalues([mOptimalDVsTMCMC(:,1), mOptimalDVsTMCMC(:,2)]);
save(fullfile(strResultFolder, 'mOptimalDVsTMCMC'), 'mOptimalDVsTMCMC')
save(fullfile(strResultFolder, 'mfEigenvaluesOptimalTMCMC'), 'mfEigenvaluesOptimalTMCMC')

% Get results from queue
mfQueueResultsTMCMC  = [];
while true
    [dataTMCMC, bDataTMCMC] = poll(mfQueueTMCMC);
    mfQueueResultsTMCMC = [mfQueueResultsTMCMC; dataTMCMC];
    if ~bDataTMCMC
        break
    end
end
mfSamplesTMCMC = mfQueueResultsTMCMC(:,1:2);
mfEigenvaluesTMCMC = mfQueueResultsTMCMC(:,3:end);


%% Monte Carlo global optimization (MCGO)
fprintf('\nCalculation using MCGO\n');

% Settings for optimization runs
optimizerMCGO = GlobalPattern();
optOptionsMCGO = struct('maxEvals', 500, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);
optParamsMCGO  = struct('nTrack', 20);

% Set number of Monte Carlo samples
nMCSamplesMCGO = 1000;

% Calculate noisy eigenvalues (= input problem)
mfEigenvaluesNoisyMCGO = zeros(nMCSamplesMCGO, nDVs);
for iDV = 1:nDVs
    mfEigenvaluesNoisyMCGO(:,iDV) = normrnd(afEigenvaluesCorrect(iDV), ...
        afSigmaNoise(iDV), nMCSamplesMCGO, 1);
end

% Model updating for each sample
mOptimalDVsMCGO = [];
cSamplesMCGO = {};
cEigenvaluesMCGO = {};
parfor iMCSampleMCGO = 1:nMCSamplesMCGO
    afEigenvaluesNoisyMCGO = mfEigenvaluesNoisyMCGO(iMCSampleMCGO, :);

    % Create data queue
    mfQueueMCGO = parallel.pool.PollableDataQueue;

    % Set model updating problem
    problemMCGO = class2DOF(afEigenvaluesNoisyMCGO, [], mfQueueMCGO, 'MCGO', []);
    probMCGO = @(x, index)problemMCGO.calcObjValue(x, index);

    % Run optimization
    [x_opt_MCGO, ~, ~, ~, ~] = optimizerMCGO.optimize(probMCGO, [], afLB, afUB, ...
        optOptionsMCGO, optParamsMCGO);
    mOptimalDVsMCGO(iMCSampleMCGO, :) = x_opt_MCGO;
    
    % Get results from queue
    mfQueueResultsMCGO  = [];
    while true
        [dataMCGO, bDataMCGO] = poll(mfQueueMCGO);
        mfQueueResultsMCGO = [mfQueueResultsMCGO; dataMCGO];
        if ~bDataMCGO
            break
        end
    end
    cSamplesMCGO{iMCSampleMCGO, 1} = mfQueueResultsMCGO(:, 1:2);
    cEigenvaluesMCGO{iMCSampleMCGO, 1} = mfQueueResultsMCGO(:, 3:end);

end
mfEigenvaluesOptimalMCGO = calcEigenvalues([mOptimalDVsMCGO(:,1), mOptimalDVsMCGO(:,2)]);
save(fullfile(strResultFolder, 'mfEigenvaluesOptimalMCGO'), 'mfEigenvaluesOptimalMCGO')
save(fullfile(strResultFolder, 'mOptimalDVsMCGO'), 'mOptimalDVsMCGO')


%% Meta-model Monte Carlo global optimization (meta-MCGO)
fprintf('\nCalculation using meta-MCGO\n');

% Settings for optimization run building the meta-model
optimizerMeta = GlobalPattern();
optOptionsMeta = struct('maxEvals', 500, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);
optParamsMeta  = struct('nTrack' , 100);

% Create data queue
mfQueueMeta  = parallel.pool.PollableDataQueue;

% Set model updating problem
problemMeta = class2DOF(afEigenvaluesCorrect, [], mfQueueMeta, 'MCGO', []);
probMeta = @(x, index)problemMeta.calcObjValue(x, index);

% Run optimization
[~, ~, ~, ~, ~] = optimizerMeta.optimize(probMeta, [], afLB, afUB, ...
    optOptionsMeta, optParamsMeta);

% Get results from queue
mfQueueResultsMeta  = [];
while true
    [dataMeta, bDataMeta] = poll(mfQueueMeta);
    mfQueueResultsMeta = [mfQueueResultsMeta; dataMeta];
    if ~bDataMeta
        break
    end
end
mfSamplesMeta = mfQueueResultsMeta(:,1:2);
mfEigenvaluesMeta = mfQueueResultsMeta(:,3:end);

% Build meta-model (using interpolation)
outputMeta = scatterModel(mfSamplesMeta, mfEigenvaluesMeta, afLB, afUB);

% Set number of Monte Carlo samples
nMCSamplesMetaMCGO = 1000;

% Generate noisy eigenvalues (= input problem)
mfEigenvaluesNoisyMetaMCGO = zeros(nMCSamplesMetaMCGO, nDVs);
for iDV = 1:nDVs
    mfEigenvaluesNoisyMetaMCGO(:,iDV) = normrnd(afEigenvaluesCorrect(iDV), ...
        afSigmaNoise(iDV), nMCSamplesMetaMCGO, 1);
end

% Settings for optimization run on meta-model
optimizerMetaMCGO = GlobalPattern();
optOptionsMetaMCGO = struct('maxEvals', 500, 'numWorkers', 1, ...
    'saveStates', false, 'outputStatus', false);
optParamsMetaMCGO  = struct('nTrack' , 20);

% Model updating for each sample on meta-model
mOptimalDVsMetaMCGO = [];
cSamplesMetaMCGO = {};
cEigenvaluesMetaMCGO = {};
parfor iMCSampleMetaMCGO = 1:size(mfEigenvaluesNoisyMCGO,1)
    afEigenvaluesNoisyMetaMCGO = mfEigenvaluesNoisyMetaMCGO(iMCSampleMetaMCGO, :);

    % Create data queue
    mfQueueMetaMCGO = parallel.pool.PollableDataQueue;

    % Set model updating problem
    problemMetaMCGO = class2DOF(afEigenvaluesNoisyMetaMCGO, [], ...
        mfQueueMetaMCGO, 'meta-MCGO', @outputMeta.predict);
    probMetaMCGO = @(x)problemMetaMCGO.calcObjValue(x);

    % Run optimization
    [x_opt_MetaMCGO, ~, ~, ~, ~] = optimizerMetaMCGO.optimize(probMetaMCGO, [], ...
        afLB, afUB, optOptionsMetaMCGO, optParamsMetaMCGO);
    mOptimalDVsMetaMCGO(iMCSampleMetaMCGO, :) = x_opt_MetaMCGO;
    
    % Get results from queue
    mfQueueResultsMetaMCGO  = [];
    while true
        [dataMetaMCGO, bDataMetaMCGO] = poll(mfQueueMetaMCGO);
        mfQueueResultsMetaMCGO = [mfQueueResultsMetaMCGO; dataMetaMCGO];
        if ~bDataMetaMCGO
            break
        end
    end
    cSamplesMetaMCGO{iMCSampleMetaMCGO} = mfQueueResultsMetaMCGO(:,1:2);
    cEigenvaluesMetaMCGO{iMCSampleMetaMCGO} = mfQueueResultsMetaMCGO(:,3:end);

end
mfEigenvaluesOptimalMetaMCGO = calcEigenvalues([mOptimalDVsMetaMCGO(:,1), mOptimalDVsMetaMCGO(:,2)]);
save(fullfile(strResultFolder, 'mfEigenvaluesOptimalMetaMCGO'), 'mfEigenvaluesOptimalMetaMCGO')
save(fullfile(strResultFolder, 'mOptimalDVsMetaMCGO'), 'mOptimalDVsMetaMCGO')


%% Plot results
cstrMethods = {'TMCMC', 'MCGO', 'Meta-MCGO'};

plotColors = {[0 0.4470 0.7410], ...
    [0.8500 0.3250 0.0980], ...
    [0.9290 0.6940 0.1250]};


%% Plot comparison samples
for iDV = 1:nDVs
    figure
    hold on
    
    plot1 = cdfplot(mOptimalDVsTMCMC(:, iDV));
    set(plot1, 'LineWidth', 2, 'Color', plotColors{1})
    
    plot2 = cdfplot(mOptimalDVsMCGO(:, iDV));
    set(plot2, 'LineWidth', 2, 'Color', plotColors{2})
    
    plot3 = cdfplot(mOptimalDVsMetaMCGO(:, iDV));
    set(plot3, 'LineWidth', 2, 'Color', plotColors{3})
    
    xline(afDVsCorrect(1, iDV), '--k', 'LineWidth', 2)
    xline(afDVsCorrect(2, iDV), '--k', 'LineWidth', 2)

    xlabel(['{\it k}_', num2str(iDV)])
    ylabel(['F({\itk}_', num2str(iDV), ')'])
    xlim([afLB(iDV), afUB(iDV)])
    title('')
    if iDV == 2
        legend([cstrMethods, 'Correct values'])
    end
    ax = gca;
    ax.FontSize = 12;
    ax.FontName = 'Times';
    grid on
    box on

    set(gcf, 'Name', 'Comparison cdfs')
    savefig(fullfile(strResultFolder,['fig2DOF_ComparisonCDFs_k', num2str(iDV)]))
end


%% Plot samples
nPlotSamples = 200;
figure
hold on

scatter(mOptimalDVsTMCMC(1:nPlotSamples,1), mOptimalDVsTMCMC(1:nPlotSamples,2), 'd', ...
    'LineWidth', 1.3, 'Color', plotColors{1})
scatter(mOptimalDVsMCGO(1:nPlotSamples,1), mOptimalDVsMCGO(1:nPlotSamples,2), 's', ...
    'LineWidth', 1.3, 'Color', plotColors{2})
scatter(mOptimalDVsMetaMCGO(1:nPlotSamples,1), mOptimalDVsMetaMCGO(1:nPlotSamples,2), 'o', ...
    'LineWidth', 1.3, 'Color', plotColors{3})
plot(afDVsCorrect(1,1), afDVsCorrect(1,2), 'k +', 'MarkerSize', 10, 'LineWidth', 2)
plot(afDVsCorrect(2,1), afDVsCorrect(2,2), 'k +', 'MarkerSize', 10, 'LineWidth', 2)

xlabel('{\it k}_1')
ylabel('{\it k}_2')
xlim([afLB(1), afUB(1)])
ylim([afLB(2), afUB(2)])
title('')
grid on
box on
legend([cstrMethods, 'Correct values'])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Name', 'Comparison samples')
savefig(fullfile(strResultFolder,'fig2DOF_DVspace'))


%% Plot eigenvalues
nPlotSamples = 200;
figure
hold on

scatter(mfEigenvaluesOptimalTMCMC(1:nPlotSamples,1), mfEigenvaluesOptimalTMCMC(1:nPlotSamples,2), 'd', ...
    'LineWidth', 1.3, 'Color', plotColors{1})
scatter(mfEigenvaluesOptimalMCGO(1:nPlotSamples,1), mOptimalDVsMCGO(1:nPlotSamples,2), 's', ...
    'LineWidth', 1.3, 'Color', plotColors{2})
scatter(mfEigenvaluesOptimalMetaMCGO(1:nPlotSamples,1), mOptimalDVsMetaMCGO(1:nPlotSamples,2), 'o', ...
    'LineWidth', 1.3, 'Color', plotColors{3})
plot(afEigenvaluesCorrect(1), afEigenvaluesCorrect(2), 'k +', 'MarkerSize', 10, 'LineWidth', 2)

t = -pi:0.01:pi;
x1 = afEigenvaluesCorrect(1) + afSigmaNoise(1) * cos(t);
y1 = afEigenvaluesCorrect(2) + afSigmaNoise(2) * sin(t);
x2 = afEigenvaluesCorrect(1) + 2*afSigmaNoise(1) * cos(t);
y2 = afEigenvaluesCorrect(2) + 2*afSigmaNoise(2) * sin(t);
plot(x1, y1, 'k', 'LineWidth', 1.3)
plot(x2, y2, '--k', 'LineWidth', 1.3)

xlabel('\lambda_1')
ylabel('\lambda_2')
xlim([2.9,3.7])
ylim([0.05,0.4])
title('')
grid on
box on
legend([cstrMethods, 'Correct value', '1\bf{\sigma}_\lambda', '2\bf{\sigma}_\lambda'])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';

set(gcf, 'Name', 'Comparison eigenvalues')
savefig(fullfile(strResultFolder,'fig2DOF_OFVspace'))

