classdef GOUIncrementCufaroSabino < GOUIncrementMethod
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations)
            rate = lambdaExponential;
            shape = lambdaPoisson / kappaOU;
            a = exp(-kappaOU * timeStep);
            increment = GOUIncrementCufaroSabino.PolyaMixture(shape, rate, a, numberOfSimulations);
        end
    end
    
    methods (Static)
        
        function gammaRemainder = PolyaMixture(shape, rate, a, numberOfSimulations)
            
            % Note that in the signature of rnbinom the probability is defines as 1 - p wrt the paper
            negative_binomial = nbinrnd(shape, a,  [numberOfSimulations, 1]);
            
            gammaRemainder = gamrnd(negative_binomial, a/rate, [numberOfSimulations, 1]);
        end
    end
end
    
