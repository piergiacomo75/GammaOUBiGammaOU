%' CalculateMixtureCoefficients
%'
%' Calculate the coefficients of the formal mixture
%'
%'
%'  a $0\le a \le 1$ it comes from the self-decomposability of the gamma random variables
%'  shape  positive number equivalent to the shape parameter of the gamma distribution
%'  numCoefficients maximum number of terms in the mixture
%'
%' return vector of the coefficients (terms) of the mixture
%'

function mixtureCoefficients = CalculateMixtureCoefficients(a, shape, numCoefficients)
  k = 0:numCoefficients - 1;
  term1 = gamma(shape + 1);
  term2 = gamma(shape - k + 1);
  term3 = gamma(k + 1);
  term4 = a.^(shape - k);
  term5 = (1-a).^k;
  mixtureCoefficients = term1 .* term4 .* term5 ./ (term2 .* term3);
end
