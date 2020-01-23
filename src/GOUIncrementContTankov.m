classdef GOUIncrementContTankov < GOUIncrementMethod    
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations)
            
            increment = GOUIncrementContTankov.ExponentialCompoundPoissonOU(lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations);
        end
    end
    
    methods (Static)
        function exponentialCompoundPoissonOU = ExponentialCompoundPoissonOU(lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations)
            
            exponentialCompoundPoissonOU = zeros(numberOfSimulations, 1);
            % Generate a N(t)
            poissonRV = poissrnd(lambdaPoisson * timeStep,  [numberOfSimulations, 1]);
            for n = 1:numberOfSimulations
                % Given N(t) generate N(t) uniforms in [0, t]
                uniformRVs = rand(poissonRV(n, 1), 1) * timeStep;
                % Sort the uniforms
                jumpTimes = sort(uniformRVs);
                % Generate the exponential jump sizes
                exponentialRVs = exprnd(1/lambdaExponential, [poissonRV(n, 1), 1]);
                exponentialCompoundPoissonOU(n, 1) =  exp(-kappaOU * timeStep) * dot(exp(kappaOU * jumpTimes), exponentialRVs);
            end
        end
    end
end

