% SCRIPT THAT REPRODUCES FIGURE 1 AND TABLE 2
clear
close all
clc

addpath(genpath('../src'))

lambdaPoisson = 1;
lambdaLaplace = 1;
kappaOU = 0.5;
x_0 = 10;

bgouSymmetric = GOUBGOU(lambdaPoisson, lambdaLaplace, kappaOU, x_0);

timeStep = 1/4;
timeGrid = 0:timeStep:1;
shape = lambdaPoisson / (2 * kappaOU);
a = exp(-kappaOU * timeGrid(2:end));
numberOfCoefficients = 40;
% numberOfSimulations =  [1000, 10000, 50000, 100000, 500000, 1000000]';
numberOfSimulations = 10; [10000, 40000, 160000, 640000, 2560000]';
mean_=  x_0 * a;
variance_ = 2 * (1-a.^2) * shape / lambdaLaplace^2;
skewness_ = 0;
kurtosis_ = (1+a.^2) ./ (1-a.^2) * 3 / shape + 3;
exponentialMean = ((lambdaLaplace.^2 + a.^2) ./ (lambdaLaplace.^2 + 1)) .^(lambdaPoisson / (2 * kappaOU));

timeAlgorithmContTankov = zeros(size(numberOfSimulations));
timeAlgorithmCufaroSabino = zeros(size(numberOfSimulations));
timeAlgorithmCufaroSabinoRejection = zeros(size(numberOfSimulations));
timeAlgorithmQDZ = zeros(size(numberOfSimulations));

ratioTimes_CT_CS = zeros(size(numberOfSimulations));
ratioTimes_CT_QDZ = zeros(size(numberOfSimulations));
ratioTimes_QDZ_CS = zeros(size(numberOfSimulations));
ratioTimes_CSRejection_QDZ = zeros(size(numberOfSimulations));

meanAlgorithm_CT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
varianceAlgorithm_CT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
skewnessAlgorithm_CT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
kurtosisAlgorithm_CT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
exponentialMeanAlgorithm_CT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
rmseExpMeanCT = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

meanAlgorithm_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
varianceAlgorithm_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
skewnessAlgorithm_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
kurtosisAlgorithm_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
exponentialMeanAlgorithm_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
rmseExpMeanCS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

meanAlgorithm_CSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
varianceAlgorithm_CSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
skewnessAlgorithm_CSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
kurtosisAlgorithm_CSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
exponentialMeanAlgorithm_CSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
rmseExpMeanCSRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

meanAlgorithm_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
varianceAlgorithm_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
skewnessAlgorithm_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
kurtosisAlgorithm_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
exponentialMeanAlgorithm_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
rmseExpMeanQDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

for n = 1:size(numberOfSimulations)
    % Algorithm_CT
    tic;
    algorithm_CT = bgouSymmetric.SimulateTrajectory(timeGrid, numberOfSimulations(n), "bgousymmetricincrementconttankov", numberOfCoefficients);
    algorithm_CT(1, :) = [];
    timeAlgorithmContTankov(n) = toc;
    
    meanAlgorithm_CT(n, :) = mean(algorithm_CT, 2);
    varianceAlgorithm_CT(n, :) = var(algorithm_CT, 0, 2);
    skewnessAlgorithm_CT(n, :) = skewness(algorithm_CT, 0, 2);
    kurtosisAlgorithm_CT(n, :) = kurtosis(algorithm_CT, 0, 2);
%     exponentialMeanAlgorithm_CT(n) = mean(exp(algorithm_CT(end, :)));

    clear algorithm_CT
    
    % Algorithm_CS
    tic
    algorithm_CS = bgouSymmetric.SimulateTrajectory(timeGrid, numberOfSimulations(n), "bgousymmetricincrementcufarosabino", numberOfCoefficients);
        algorithm_CS(1, :) = [];
    timeAlgorithmCufaroSabino(n) = toc;
    
    meanAlgorithm_CS(n, :) = mean(algorithm_CS, 2);
    varianceAlgorithm_CS(n, :) = var(algorithm_CS, 0, 2);
    skewnessAlgorithm_CS(n, :) = skewness (algorithm_CS, 0, 2);
    kurtosisAlgorithm_CS(n, :) = kurtosis(algorithm_CS, 0, 2);
%     exponentialMeanAlgorithm_CS(n, :) = mean(exp(algorithm_CS), 2);

    clear algorithm_CS
    
    % Algorithm Acceptance Rejection
    tic
    algorithm_CS_Rejection = bgouSymmetric.SimulateTrajectory(timeGrid, numberOfSimulations(n), "bgousymmetricincrementcufarosabinorejection", numberOfCoefficients);
    algorithm_CS_Rejection(1, :) = [];
    timeAlgorithmCufaroSabinoRejection(n) = toc;
    
    meanAlgorithm_CSRejection(n, :) = mean(algorithm_CS_Rejection, 2);
    varianceAlgorithm_CSRejection(n, :) = var(algorithm_CS_Rejection, 0, 2);
    skewnessAlgorithm_CSRejection(n, :) = skewness(algorithm_CS_Rejection, 0, 2);
    kurtosisAlgorithm_CSRejection(n, :) = kurtosis(algorithm_CS_Rejection, 0, 2);
%     exponentialMeanAlgorithm_CSRejection(n, :) = mean(exp(algorithm_CS_Rejection, 2)););

    clear algorithm_CS_Rejection
    
    % Algorithm_QDZ
    tic
    algorithm_QDZ = bgouSymmetric.SimulateTrajectory(timeGrid, numberOfSimulations(n), "bgousymmetricincrementqdz", numberOfCoefficients);
   algorithm_QDZ(1, :) = [];
    timeAlgorithmQDZ(n) = toc;
    
    meanAlgorithm_QDZ(n, :) = mean(algorithm_QDZ, 2);
    varianceAlgorithm_QDZ(n, :) = var(algorithm_QDZ, 0, 2);
    skewnessAlgorithm_QDZ(n, :) = skewness(algorithm_QDZ, 0, 2);
    kurtosisAlgorithm_QDZ(n, :) = kurtosis(algorithm_QDZ, 0, 2);
%     exponentialMeanAlgorithm_QDZ(n, :) = mean(exp(algorithm_QDZ, 2));    
    
    clear algorithm_QDZ
    
    ratioTimes_CT_CS(n) = timeAlgorithmContTankov(n) / timeAlgorithmCufaroSabino(n);
    ratioTimes_CT_QDZ(n) = timeAlgorithmContTankov(n) / timeAlgorithmQDZ(n);
    ratioTimes_QDZ_CS(n) = timeAlgorithmCufaroSabino(n) / timeAlgorithmQDZ(n);
    ratioTimes_CSRejection_QDZ(n) = timeAlgorithmCufaroSabinoRejection(n) / timeAlgorithmQDZ(n);
end

rmpath(genpath('../src'))



