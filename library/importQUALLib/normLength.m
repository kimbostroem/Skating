function NormVec = normLength(Vector,Dim)
% This function will return a vector of unit length.
% Vector can be multidimensional. Per default, the normalization is made
% over the first dimension, so if Vector is of the form [3,ns].
% Marc de Lussanet, WWU Muenster, 20200530

if nargin > 1
    NormVec = Vector./sqrt(sum(Vector.^2, Dim));
else
    NormVec = Vector./sqrt(sum(Vector.^2));
end

end
