classdef GOUIncrementCufaroSabinoRejection < GOUIncrementMethod
    
    properties
        numberOfCoefficients
    end
    
    methods
        function obj = GOUIncrementCufaroSabinoRejection(numberOfCoefficients)
            obj.numberOfCoefficients = numberOfCoefficients;
        end
    end
    
    methods
        function increment = GenerateIncrement(obj, lambdaPoisson, lambdaExponential, kappaOU, timeStep, numberOfSimulations)
            
            shape = lambdaPoisson / kappaOU;
            a = exp(-kappaOU * timeStep);
            increment = GOUIncrementCufaroSabinoRejection.randAdditionalSD(a, shape, lambdaExponential, numberOfSimulations, obj.numberOfCoefficients);
        end
    end
    
    methods(Static, Access = protected)
        
        negativeMixture = SDNegativeMixture(x, a, shape, numCoefficients)
        positiveMixture = SDPositiveMixture(x, a, shape, numCoefficients)
        conditionalErlang = randSDPositiveMixture(a, shape, numCoefficients, numSimulations)
        mixtureCoefficients = CalculateMixtureCoefficients(a, shape, numCoefficients)
        negativeCoefficients = CalculateNegativeMixtureCoefficients(a, shape, numCoefficients)
        positiveCoefficients = CalculatePositiveMixtureCoefficients(a, shape, numCoefficients)
        
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
        
        function returnVector = randAdditionalSD(a, shape, lambdaExponential, numberOfSimulations, numberOfCoefficients)
            
            returnVector = zeros(numberOfSimulations, 1);
            conditionalErlang = GOUIncrementCufaroSabinoRejection.randSDPositiveMixture(a, shape, numberOfCoefficients, numberOfSimulations);
            acceptIndicesZeros = (conditionalErlang == 0);
            returnVector(acceptIndicesZeros) = 0;
            
            if all(acceptIndicesZeros)
                return
            end
            numberAcceptanceRejectionSimulations    = sum(~acceptIndicesZeros);
            numberOfSimulationPerLoop               = numberAcceptanceRejectionSimulations;
            
            accept                      = false;
            acceptIndices               = acceptIndicesZeros;
            newAcceptRejectionIndices   = acceptIndicesZeros;
            
            conditionalErlangNonZero = conditionalErlang(~acceptIndices);
            
            while accept == false
                
                uniformRV                   = rand(numberOfSimulationPerLoop,1);
                ratio                       = 1 - ...
                    GOUIncrementCufaroSabinoRejection.SDNegativeMixture(conditionalErlangNonZero, a, shape, numberOfCoefficients) ./ ...
                    GOUIncrementCufaroSabinoRejection.SDPositiveMixture(conditionalErlangNonZero, a, shape, numberOfCoefficients);

                % If conditionalErlangNonZero is zero the ratio is nan but
                % then is acceptable. This condition may be
                % verified within the while if accept is not verified at
                % the first step
                ratio(isnan(ratio)) = 1; 
                acceptRejectionIndices = (uniformRV <= ratio);
                isFalseAcceptIndices = ~acceptIndices; % equivalent to acceptIndices == false;
                newAcceptRejectionIndices(isFalseAcceptIndices) = acceptRejectionIndices;
                                
                idxReturnVector = newAcceptRejectionIndices & (~acceptIndices); % equivalent to newAcceptRejectionIndices & (acceptIndices == false)
                returnVector(idxReturnVector) = conditionalErlangNonZero(acceptRejectionIndices) / lambdaExponential;
                
                % Update new indices for accept
                acceptIndices               = acceptIndices | newAcceptRejectionIndices;                                
                
                totalNumberOfAccepted       = sum(acceptIndices);
                numberOfSimulationPerLoop   = numberOfSimulations - totalNumberOfAccepted;
                
                accept                      = (totalNumberOfAccepted >= numberOfSimulations);                
                newAcceptRejectionIndices   = false(numberOfSimulations, 1);
                conditionalErlangNonZero    = GOUIncrementCufaroSabinoRejection.randSDPositiveMixture(a, shape, numberOfCoefficients, numberOfSimulationPerLoop);
                %     conditionalErlangNonZero    = conditionalErlang(conditionalErlang > 0);
            end
        end
        
        %' randAdditionalSDSlow
        %'
        %' Simulates the "additional" random variable entering the self-decomposiability of the gamma law
        %'
        %'  a 0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
        %'  shape  positive number equivalent to the shape parameter of the gamma distribution
        %'  numCoefficients maximum number of terms in the formal mixture
        %'  numSimulations number of Monte Carlo simulations
        %'
        %' return vector of simulated random variables
        
        function returnVector = randAdditionalSDSlow(a, shape, lambdaExponential, numSimulations, numCoefficients)
            returnVector = zeros(numSimulations, 1);
            
            for i = 1:numSimulations
                accept = false;
                while ~accept
                    conditionalErlang = GOUIncrementCufaroSabinoRejection.randSDPositiveMixture(a, shape, numCoefficients, 1);
                    if(conditionalErlang == 0)
                        returnVector(i) = 0;
                        break
                    end
                    
                    uniformRV2 = rand(1);
                    accept = (uniformRV2 <= 1 - ...
                        GOUIncrementCufaroSabinoRejection.SDNegativeMixture(conditionalErlang, a, shape, numCoefficients) / ...
                        GOUIncrementCufaroSabinoRejection.SDPositiveMixture(conditionalErlang, a, shape, numCoefficients));
                    if(accept)
                        returnVector(i) = conditionalErlang / lambdaExponential;
                        break
                    end
                end
            end
        end
    end
end
