%%
% get GLM Canonical calcium response function
% returns time range (t), CRF time-series (crf)
% input:
%  dt                 time resolution (sec) (default:0.01)
%  gammaA             a value of gamma function (default:1.5)
%  gammaB             HRF scale (gamma b) (default: 0.9)
%  kernelSec          kernel time length (sec) (default: 5)

function [t, hrf] = getGlmCRF(dt, gammaA, gammaB, kernelSec)
    if nargin < 4, kernelSec = 5; end
    if nargin < 3, gammaB = 0.9; end
    if nargin < 2, gammaA = 1.5; end
    if nargin < 1, dt = 0.01; end
    
    t = 0:dt:kernelSec;
    hrf = gampdf(t,gammaA,gammaB);
    hrf = hrf'/sum(hrf);
end
