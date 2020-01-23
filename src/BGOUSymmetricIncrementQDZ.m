classdef BGOUSymmetricIncrementQDZ < GOUIncrementMethod
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaLaplace, kappaOU, timeStep, numberOfSimulations)
            
            increment = BGOUSymmetricIncrementQDZ.LaplaceCompoundPoissonQDZ(lambdaPoisson, lambdaLaplace, kappaOU, timeStep, numberOfSimulations);
        end
    end
    
    methods (Static)
        function laplaceCompoundPoissonQDZ = LaplaceCompoundPoissonQDZ(lambdaPoisson, ...
                lambdaLaplace, ...
                kappaOU, ...
                timeStep, ...
                numberOfSimulations)
            
            laplaceCompoundPoissonQDZ = zeros(numberOfSimulations, 1);
            % Generate a N(t)
            poissonRV = poissrnd(lambdaPoisson * timeStep,  [numberOfSimulations, 1]);
            for n = 1:numberOfSimulations
                % Given N(t)=n generate n uniforms for the random mean
                randomMean_up =  exp(-kappaOU * timeStep * rand(poissonRV(n, 1), 1)) / lambdaLaplace;
                randomMean_down =  exp(-kappaOU * timeStep * rand(poissonRV(n, 1), 1)) / lambdaLaplace;
                jumps_up = exprnd(randomMean_up, [poissonRV(n, 1), 1]);
                jumps_down = exprnd(randomMean_down, [poissonRV(n, 1), 1]);
                jumps = jumps_up - jumps_down;
                % Generate the exponential jump sizes and sum
                laplaceCompoundPoissonQDZ(n, 1) = sum(jumps);
            end
        end
    end
end


