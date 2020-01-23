% SCRIPT THAT REPRODUCES FIGURE 1 AND TABLE 2
clear
close all
clc

addpath(genpath('../src'))

lambdaPoisson = 1;
lambdaExponential = 1;
kappaOU = 0.5;
x_0 = 10;
timeStep = 1/4;
timeGrid = 0:timeStep:1;
shape = lambdaPoisson / kappaOU;
a = exp(-kappaOU * timeGrid(2:end));

gou = GOUBGOU(lambdaPoisson, lambdaExponential, kappaOU, x_0);
     
numberOfCoefficients = 40;
% numberOfSimulations = 100000; [1000, 10000, 50000, 100000, 500000, 1000000]';
numberOfSimulations =10;[10000, 40000, 160000, 640000, 2560000]';
mean_= (1-a) * shape / lambdaExponential + x_0 * a;
variance_ = (1-a.^2) * shape / lambdaExponential^2;
skewness_ = (1-a.^3) ./ (1-a.^2).^(1.5) * 2 / sqrt(shape);
kurtosis_ = (1+a.^2) ./ (1-a.^2) * 6 / shape + 3;
exponentialMean = exp(-x_0 * a) .* ((lambdaExponential - a) / (lambdaExponential - 1)) .^ (lambdaPoisson / kappaOU);

timeAlgorithmContTankov = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
timeAlgorithmCufaroSabino = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
timeAlgorithmCufaroSabinoRejection = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
timeAlgorithmQDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

ratioTimes_CT_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
ratioTimes_CT_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
ratioTimes_QDZ_CS = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);
ratioTimes_CSRejection_QDZ = zeros(numel(numberOfSimulations), numel(timeGrid) - 1);

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

for n = 1:numel(numberOfSimulations)
    % Algorithm_CT
    tic;
    algorithm_CT = gou.SimulateTrajectory(timeGrid, numberOfSimulations(n), "gouincrementconttankov");
    algorithm_CT(1, :) = [];
    timeAlgorithmContTankov(n) = toc;
    
    meanAlgorithm_CT(n, :) = mean(algorithm_CT, 2);
    varianceAlgorithm_CT(n, :) = var(algorithm_CT, 0, 2);
    skewnessAlgorithm_CT(n, :) = skewness(algorithm_CT, 0, 2);
    kurtosisAlgorithm_CT(n, :) = kurtosis(algorithm_CT, 0, 2);
%     exponentialMeanAlgorithm_CT(n, :) = mean(exp(algorithm_CT, 2));
  
    clear algorithm_CT
    
    % Algorithm_CS
    tic
    algorithm_CS = gou.SimulateTrajectory(timeGrid, numberOfSimulations(n), "gouincrementcufarosabino");
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
    algorithm_CS_Rejection = gou.SimulateTrajectory(timeGrid, numberOfSimulations(n), "gouincrementcufarosabinorejection", numberOfCoefficients);
    algorithm_CS_Rejection(1, :) = [];
    timeAlgorithmCufaroSabinoRejection(n) = toc;
    
    meanAlgorithm_CSRejection(n, :) = mean(algorithm_CS_Rejection, 2);
    varianceAlgorithm_CSRejection(n, :) = var(algorithm_CS_Rejection, 0, 2);
    skewnessAlgorithm_CSRejection(n, :) = skewness(algorithm_CS_Rejection, 0, 2);
    kurtosisAlgorithm_CSRejection(n, :) = kurtosis(algorithm_CS_Rejection, 0, 2);
%     exponentialMeanAlgorithm_CSRejection(n, :) = mean(exp(algorithm_CS_Rejection, 2));
% 
    clear algorithm_CS_Rejection
    
    % Algorithm_QDZ
    tic
    algorithm_QDZ = gou.SimulateTrajectory(timeGrid, numberOfSimulations(n), "gouincrementqdz");
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

