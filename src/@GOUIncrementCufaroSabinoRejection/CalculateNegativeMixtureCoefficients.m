%' CalculateNegativeMixtureCoefficients
%'
%' %' Calculate the negative coefficients of the formal mixture
%'
%'  a 0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
%'  shape  positive number equivalent to the shape parameter of the gamma distribution
%'  numCoefficients maximum number of terms in the mixture
%'
%' @return vector of the positive coefficients (terms) of the overall mixture

function negativeCoefficients = CalculateNegativeMixtureCoefficients(a, shape, numCoefficients)

  mixtureCoefficients = GOUIncrementCufaroSabinoRejection.CalculateMixtureCoefficients(a, shape, numCoefficients);
  negativeCoefficients = -min(mixtureCoefficients,  0);
  
end
