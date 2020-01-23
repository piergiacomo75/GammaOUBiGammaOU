classdef BGOUSymmetricIncrementCufaroSabino < GOUIncrementMethod
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaLaplace, kappaOU, timeStep, numberOfSimulations)
            rate = lambdaLaplace;
            shape = lambdaPoisson / (2 * kappaOU);
            a = exp(-kappaOU * timeStep);
            increment = BGOUSymmetricIncrementCufaroSabino.PolyaMixtureErlangDifference(shape, rate, a, numberOfSimulations);
        end
    end
    
    methods (Static)
        
        function laplaceRemainder = PolyaMixtureErlangDifference(shape, lambdaLaplace, a, numberOfSimulations)
            
            % Note that in the signature of rnbinom the probability is defines as 1 - p wrt the paper
            negative_binomial = nbinrnd(shape, a^2,  [numberOfSimulations, 1]);
            
            laplaceRemainder = gamrnd(negative_binomial, a/lambdaLaplace, [numberOfSimulations, 1]) - gamrnd(negative_binomial, a/lambdaLaplace, [numberOfSimulations, 1]);
            
        end
    end
end
    
