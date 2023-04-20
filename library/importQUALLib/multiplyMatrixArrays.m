function ABvector = multiplyMatrixArrays(Aarray,Barray)
% vector-multiply two arrays of matrices
% Oct 22 2013 at 15:06 by Ruben

[mA,nA,oA] = size(Aarray);
[mB,nB,oB] = size(Barray);

% preallocate the output
ABvector = zeros(mA,nB,oA);

%error message if the dimensions are not appropriate
if nA ~= mB || oA ~= oB
    error('Inconsistent matrix dimensions na%d~=mb%d or oa%d~=ob%d',nA,mB,oA,oB);
end

% if statement minimizes for loops by looping the smallest matrix dimension
if mA > nB
    Bp = ABvector;
    for j = 1:nB
        Bp(j,:,:) = Barray(:,j,:);
        ABvector(:,j,:) = sum(Aarray.*repmat(Bp(j,:,:),[mA,1]),2);
    end
else
    Ap = ABvector;
    for i = 1:mA
        Ap(:,i,:) = Aarray(i,:,:);
        ABvector(i,:,:) = sum(repmat(Ap(:,i,:),[1,nB]).*Barray,1);
    end
end

end
