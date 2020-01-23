classdef GOUIncrementQDZ < GOUIncrementMethod
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations)
            
            increment = GOUIncrementQDZ.ExponentialCompoundPoissonQDZ(lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations);
        end
    end
    
    methods (Static)
        function exponentialCompoundPoissonQDZ = ExponentialCompoundPoissonQDZ(lambdaPoisson, ...
                lambdaExponential, ...
                kappaOU, ...
                timeStep, ...
                numberOfSimulations)
                        
            exponentialCompoundPoissonQDZ = zeros(numberOfSimulations, 1);
            % Generate a N(t)
            poissonRV = poissrnd(lambdaPoisson * timeStep,  [numberOfSimulations, 1]);
            for n = 1:numberOfSimulations
                % Given N(t)=n generate n uniforms for the random mean
                randomMean =  exp(-kappaOU * timeStep * rand(poissonRV(n, 1), 1)) / lambdaExponential;
                % Generate the exponential jump sizes and sum
                exponentialCompoundPoissonQDZ(n, 1) = sum(exprnd(randomMean, [poissonRV(n, 1), 1]));
            end
        end
    end
end

