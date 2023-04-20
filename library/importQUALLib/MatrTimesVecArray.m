function OutArray = MatrTimesVecArray(Marray,Varray)
% multiply arrays of Matrix(n,m,nsamp) * Vector(n, nsamp)
% 3.6.2020 MdL, WWU Muenster

SzM = size(Marray);
SzV = size(Varray);
nSamp = SzV(end);
Dim   = SzM(1);

% Check Dimensions
if length(SzM)~=3 || length(SzV)~=2 || nSamp ~= SzM(3) || Dim~=SzM(1) || Dim~=SzM(1) || Dim~=SzV(1)
    error('Inconsistent matrix dimensions Marray %d %d %d or Varray %d %d',SzM(1),SzM(2),SzM(3),SzV(1),SzV(2));
end

% preallocate the output
OutArray = Varray;

% multiply
for j = 1:Dim
    OutArray(j,:) = sum(squeeze(Marray(j,:,:)).*Varray);
end

end
