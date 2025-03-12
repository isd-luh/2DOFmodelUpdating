function EVsamplesRDMU = generateEVsamples(strSampler, nDVs, nSamplesRDMU, ...
    strInputDist, meanEigenvalues, sigmaEigenvalues, scaleWB, shapeWB)

% Choose sampling for deterministic sampling methods
if isequal(strSampler,'HS')
    sampling = net(haltonset(nDVs * 2, 'Skip',1), nSamplesRDMU);
elseif isequal(strSampler,'SS')
    sampling = net(sobolset(nDVs * 2, 'Skip',1), nSamplesRDMU);
end

% Generate eigenvalue samples (dependent on choice of sampling method)
if isequal(strSampler,'HS') || isequal(strSampler,'SS')
    if isequal(strInputDist,'Gaussian')
        EVsamplesRDMU = norminv(sampling(:, 1:nDVs), meanEigenvalues, sigmaEigenvalues);
    elseif isequal(strInputDist,'Weibull')
        EVsamplesRDMU = wblinv(sampling(:, 1:nDVs), scaleWB, shapeWB);
    end
elseif isequal(strSampler,'MC')
    if isequal(strInputDist,'Gaussian')
        EVsamplesRDMU = normrnd(repmat(meanEigenvalues, nSamplesRDMU, 1), ...
            repmat(sigmaEigenvalues, nSamplesRDMU, 1));
    elseif isequal(strInputDist,'Weibull')
        EVsamplesRDMU = wblrnd(repmat(scaleWB, nSamplesRDMU, 1), ...
            repmat(shapeWB, nSamplesRDMU, 1));
    end
end
end