% SCRIPT THAT REPRODUCES FIGURES AND TABLES
clear
close all
clc

addpath(genpath('../src'))

lambdaPoisson = 10;
lambdaLaplace = 3;
kappaOU = 36;
timeStep = 1/365;
shape = lambdaPoisson / (2 * kappaOU);
a = exp(-kappaOU * timeStep);
numberOfCoefficients = 40;
% numberOfSimulations =  [1000, 10000, 50000, 100000, 500000, 1000000]';
numberOfSimulations =  10;[10000, 40000, 160000, 640000, 2560000]';
mean_= 0;
variance_ = 2 * (1-a^2) * shape / lambdaLaplace^2;
skewness_ = 0;
kurtosis_ = (1+a^2) / (1-a^2) * 3 / shape + 3;
exponentialMean = ((lambdaLaplace^2 + a^2) / (lambdaLaplace^2 + 1)) ^(lambdaPoisson / (2 * kappaOU));

timeAlgorithmContTankov = zeros(size(numberOfSimulations));
timeAlgorithmCufaroSabino = zeros(size(numberOfSimulations));
timeAlgorithmCufaroSabinoRejection = zeros(size(numberOfSimulations));
timeAlgorithmQDZ = zeros(size(numberOfSimulations));

ratioTimes_CT_CS = zeros(size(numberOfSimulations));
ratioTimes_CT_QDZ = zeros(size(numberOfSimulations));
ratioTimes_QDZ_CS = zeros(size(numberOfSimulations));
ratioTimes_CSRejection_QDZ = zeros(size(numberOfSimulations));

meanAlgorithm_CT = zeros(size(numberOfSimulations));
varianceAlgorithm_CT = zeros(size(numberOfSimulations));
skewnessAlgorithm_CT = zeros(size(numberOfSimulations));
kurtosisAlgorithm_CT = zeros(size(numberOfSimulations));
moment5Algorithm_CT = zeros(size(numberOfSimulations));
exponentialMeanAlgorithm_CT = zeros(size(numberOfSimulations));
rmseExpMeanCT = zeros(size(numberOfSimulations));

meanAlgorithm_CS = zeros(size(numberOfSimulations));
varianceAlgorithm_CS = zeros(size(numberOfSimulations));
skewnessAlgorithm_CS = zeros(size(numberOfSimulations));
kurtosisAlgorithm_CS = zeros(size(numberOfSimulations));
moment5Algorithm_CS = zeros(size(numberOfSimulations));
exponentialMeanAlgorithm_CS = zeros(size(numberOfSimulations));
rmseExpMeanCS = zeros(size(numberOfSimulations));

meanAlgorithm_CSRejection = zeros(size(numberOfSimulations));
varianceAlgorithm_CSRejection = zeros(size(numberOfSimulations));
skewnessAlgorithm_CSRejection = zeros(size(numberOfSimulations));
kurtosisAlgorithm_CSRejection = zeros(size(numberOfSimulations));
moment5Algorithm_CSRejection = zeros(size(numberOfSimulations));
exponentialMeanAlgorithm_CSRejection = zeros(size(numberOfSimulations));
rmseExpMeanCSRejection = zeros(size(numberOfSimulations));

meanAlgorithm_QDZ = zeros(size(numberOfSimulations));
varianceAlgorithm_QDZ = zeros(size(numberOfSimulations));
skewnessAlgorithm_QDZ = zeros(size(numberOfSimulations));
kurtosisAlgorithm_QDZ = zeros(size(numberOfSimulations));
moment5Algorithm_QDZ = zeros(size(numberOfSimulations));
exponentialMeanAlgorithm_QDZ = zeros(size(numberOfSimulations));
rmseExpMeanQDZ = zeros(size(numberOfSimulations));

for n = 1:size(numberOfSimulations)
    % Algorithm_CT
    tic;
    algorithm_CT = BGOUSymmetricIncrementContTankov.LaplaceCompoundPoissonOU(lambdaPoisson, lambdaLaplace, kappaOU, timeStep, numberOfSimulations(n));
    
    timeAlgorithmContTankov(n) = toc;
    
%     meanAlgorithm_CT(n) = mean(algorithm_CT);
    varianceAlgorithm_CT(n) = var(algorithm_CT);
%     skewnessAlgorithm_CT(n) = skewness(algorithm_CT);
    kurtosisAlgorithm_CT(n) = kurtosis(algorithm_CT);
    exponentialMeanAlgorithm_CT(n) = mean(exp(algorithm_CT));

    clear algorithm_CT
    
    % Algorithm_CS
    tic
    algorithm_CS = BGOUSymmetricIncrementCufaroSabino.PolyaMixtureErlangDifference(shape, lambdaLaplace, a, numberOfSimulations(n));
    
    timeAlgorithmCufaroSabino(n) = toc;
    
%     meanAlgorithm_CS(n) = mean(algorithm_CS);
    varianceAlgorithm_CS(n) = var(algorithm_CS);
%     skewnessAlgorithm_CS(n) = skewness(algorithm_CS);
    kurtosisAlgorithm_CS(n) = kurtosis(algorithm_CS);
    exponentialMeanAlgorithm_CS(n) = mean(exp(algorithm_CS));

    clear algorithm_CS
    
    % Algorithm Acceptance Rejection
    tic
    algorithm_CS_Rejection = BGOUSymmetricIncrementCufaroSabinoRejection.randAdditionalSDLaplace(a, shape, lambdaLaplace, numberOfSimulations(n), numberOfCoefficients);
    
    timeAlgorithmCufaroSabinoRejection(n) = toc;
    
%     meanAlgorithm_CSRejection(n) = mean(algorithm_CS_Rejection);
    varianceAlgorithm_CSRejection(n) = var(algorithm_CS_Rejection);
%     skewnessAlgorithm_CSRejection(n) = skewness(algorithm_CS_Rejection);
    kurtosisAlgorithm_CSRejection(n) = kurtosis(algorithm_CS_Rejection);
    exponentialMeanAlgorithm_CSRejection(n) = mean(exp(algorithm_CS_Rejection));

    clear algorithm_CS_Rejection
    
    % Algorithm_QDZ
    tic
    algorithm_QDZ = BGOUSymmetricIncrementQDZ.LaplaceCompoundPoissonQDZ(lambdaPoisson, lambdaLaplace, kappaOU,  timeStep, numberOfSimulations(n));
    
    timeAlgorithmQDZ(n) = toc;
    
%     meanAlgorithm_QDZ(n) = mean(algorithm_QDZ);
    varianceAlgorithm_QDZ(n) = var(algorithm_QDZ);
%     skewnessAlgorithm_QDZ(n) = skewness(algorithm_QDZ);
    kurtosisAlgorithm_QDZ(n) = kurtosis(algorithm_QDZ);
    exponentialMeanAlgorithm_QDZ(n) = mean(exp(algorithm_QDZ));
    clear algorithm_QDZ
    
    ratioTimes_CT_CS(n) = timeAlgorithmContTankov(n) / timeAlgorithmCufaroSabino(n);
    ratioTimes_CT_QDZ(n) = timeAlgorithmContTankov(n) / timeAlgorithmQDZ(n);
    ratioTimes_QDZ_CS(n) = timeAlgorithmQDZ(n) / timeAlgorithmCufaroSabino(n);
    ratioTimes_CSRejection_QDZ(n) = timeAlgorithmCufaroSabinoRejection(n) / timeAlgorithmQDZ(n);
end

rmpath(genpath('../src'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1

close all
figure(1)

% Right

plot(numberOfSimulations, ...
    timeAlgorithmCufaroSabino, ...
    '-.*k', ...
    numberOfSimulations, ...
    timeAlgorithmCufaroSabinoRejection, ...
    '-*m', ...
    numberOfSimulations, ...
    timeAlgorithmContTankov, ...
    '--or', ...
    numberOfSimulations, ...
    timeAlgorithmQDZ, ...
    ':*b');

title('CPU Times', 'Color','black')
xlabel('Number of Simulations', 'Color','black') % x-axis label
ylabel('CPU time (sec)', 'Color','black') % y-axis labell
legend('Alg 1', 'Alg 3', 'Alg 4', 'Alg 4', 'Location','northwest')

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XColor = 'black';
ax.YColor = 'black';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2

figure(2)

% Right

plot(numberOfSimulations, ...
    ratioTimes_CT_CS, ...
    '--or', ...
    numberOfSimulations, ...
    1./ratioTimes_CT_QDZ, ...
    ':*b', ...
    numberOfSimulations, ...
    ratioTimes_QDZ_CS, ...
    '-.*k');

title('Ratios of CPU Times', 'Color','black')
xlabel('Number of Simulations', 'Color','black') % x-axis label
ylabel('CPU time ratio', 'Color','black') % y-axis labell
legend('Alg 4 / Alg 1', 'Alg 4 / Alg 5', 'Alg 5 / Alg 1', 'Location','east')

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XColor = 'black';
ax.YColor = 'black';


