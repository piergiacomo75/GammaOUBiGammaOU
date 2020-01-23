%' SDPositiveMixture
%'
%' Positive part of the mixture
%'
%'  x argument of the positive part of the formal mixture
%'  a 0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
%'  shape  positive number equivalent to the shape parameter of the gamma distribution
%'  numCoefficients maximum number of terms in the mixture
%'
%' @return value of the positive part of the formal mixture

function positiveMixture = SDPositiveMixtureBiGamma (x, a, shape, numCoefficients)
  mixtureCoefficients = GOUIncrementCufaroSabinoRejection.CalculatePositiveMixtureCoefficients(a, shape, numCoefficients);
  indices = 0:(numCoefficients - 1);
  fun = @(y, z) abs(y).^(z-1/2) .* exp(-y) .* besselk(z-1/2, abs(y)) .* sqrt(2/pi) ./ gamma(z);
  
  positiveErlang = bsxfun(fun, x, indices);
  positiveMixture =  positiveErlang * mixtureCoefficients';

end
