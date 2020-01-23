classdef BGOUSymmetricIncrementCufaroSabinoRejection < GOUIncrementCufaroSabinoRejection
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaLaplace, kappaOU, timeStep, numberOfSimulations)
            rate = lambdaLaplace;
            shape = lambdaPoisson / (2 * kappaOU);
            a = exp(-kappaOU * timeStep);
            increment = BGOUSymmetricIncrementCufaroSabinoRejection.randAdditionalSDLaplace(a, shape, rate, numberOfSimulations, obj.numberOfCoefficients);
        end
    end
    
    methods(Static, Access = protected)
        
        negativeMixture = SDNegativeMixtureBiGamma(x, a, shape, numberOfCoefficients)
        positiveMixture = SDPositiveMixtureBiGamma(x, a, shape, numberOfCoefficients)
        conditionalBigamma = randSDPositiveMixtureBiGamma(a, shape, numCoefficients, numberOfSimulations)
    end
    
    methods (Static)
        
        %' randAdditionalSD
        %'
        %' Simulates the "additional" random variable entering the self-decomposiability of the gamma law
        %'
        %'  a 0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
        %'  shape  positive number equivalent to the shape parameter of the gamma distribution
        %'  numCoefficients maximum number of terms in the formal mixture
        %'  numSimulations number of Monte Carlo simulations
        %'
        %' return vector of simulated random variables
        
        function returnVector = randAdditionalSDLaplace(a, shape, lambdaLaplace, numberOfSimulations, numberOfCoefficients)
            returnVector = zeros(numberOfSimulations, 1);
            
            for i = 1:numberOfSimulations
                accept = false;
                while ~accept
                    conditionalBiErlang = BGOUSymmetricIncrementCufaroSabinoRejection.randSDPositiveMixtureBiGamma(a, shape, 1, numberOfCoefficients);
                    if(conditionalBiErlang == 0)
                        returnVector(i) = 0;
                        break
                    end
                    
                    uniformRV2 = rand(1);
                    accept = (uniformRV2 <= 1 - ...
                        BGOUSymmetricIncrementCufaroSabinoRejection.SDNegativeMixtureBiGamma(conditionalBiErlang, a, shape, numberOfCoefficients) / ...
                        BGOUSymmetricIncrementCufaroSabinoRejection.SDPositiveMixtureBiGamma(conditionalBiErlang, a, shape, numberOfCoefficients));
                    if(accept)
                        returnVector(i) = conditionalBiErlang / lambdaLaplace;
                        break
                    end
                end
            end
        end
    end
end