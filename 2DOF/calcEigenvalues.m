function [mfEigenvalues] = calcEigenvalues(afDVs)
% This function calculates the eigenvalues of the 2DOF system 
% corresponding to a given set of design variables
%
% Reference: 
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% INPUT: 
% - Design variables: afDVs

% OUTPUT: 
% - Eigenvalues: mfEigenvalues

%%
% Calculate eigenvalues
mfEigenvalues(:,1) = 0.5 * ((afDVs(:,1) + 2*afDVs(:,2)) + (afDVs(:,1).^2 + 4*(afDVs(:,2).^2)).^0.5);
mfEigenvalues(:,2) = 0.5 * ((afDVs(:,1) + 2*afDVs(:,2)) - (afDVs(:,1).^2 + 4*(afDVs(:,2).^2)).^0.5);

end