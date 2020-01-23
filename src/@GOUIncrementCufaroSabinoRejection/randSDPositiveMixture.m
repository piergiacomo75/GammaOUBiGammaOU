%' rSDPositiveMixture
%'
%' Simulates a random variable having the (normalized) positive mixture as density
%'
%'  a 0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
%'  shape  positive number equivalent to the shape parameter of the gamma distribution
%'  numCoefficients maximum number of terms in the formal mixture
%'  numSimulations number of Monte Carlo simulations
%'
%' return vector of simulated random variables

function conditionalErlang = randSDPositiveMixture(a, shape, numCoefficients, numSimulations)

positiveMixtureCoefficients = GOUIncrementCufaroSabinoRejection.CalculatePositiveMixtureCoefficients(a, shape, numCoefficients);
cdfDiscrete = cumsum(positiveMixtureCoefficients) / sum(positiveMixtureCoefficients);
uniformRV1 = rand(numSimulations, 1);
 
discreteRV = arrayfun( @(x) find(x <= cdfDiscrete, 1, 'first') - 1, uniformRV1);
conditionalErlang = gamrnd(discreteRV, 1, [numSimulations, 1]);

end
