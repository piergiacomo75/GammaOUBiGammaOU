classdef BGOUSymmetricIncrementContTankov < GOUIncrementMethod
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaJumps, kappaOU, timeStep, numberOfSimulations)
            
            increment = BGOUSymmetricIncrementContTankov.LaplaceCompoundPoissonOU(lambdaPoisson, lambdaJumps, kappaOU, timeStep, numberOfSimulations);
        end
    end
    
    methods (Static)
        function laplaceCompoundPoissonOU = LaplaceCompoundPoissonOU(lambdaPoisson, ...
                lambdaLaplace, ...
                kappaOU, ...
                timeStep, ...
                numberOfSimulations)
            
            laplaceCompoundPoissonOU = zeros(numberOfSimulations, 1);
            % Generate a N(t)
            poissonRV = poissrnd(lambdaPoisson * timeStep,  [numberOfSimulations, 1]);
            for n = 1:numberOfSimulations
                % Given N(t) generate N(t) uniforms in [0, t]
                uniformRVs = rand(poissonRV(n, 1), 1) * timeStep;
                % Sort the uniforms
                jumpTimes = sort(uniformRVs);
                % Generate the exponential jump sizes
                laplaceRVs = exprnd(1/lambdaLaplace, [poissonRV(n, 1), 1]) - exprnd(1/lambdaLaplace, [poissonRV(n, 1), 1]);
                laplaceCompoundPoissonOU(n, 1) =  exp(-kappaOU * timeStep) * dot(exp(kappaOU * jumpTimes), laplaceRVs);
            end
        end
    end
end
    
