function Distance = point2LineDistance(Point,Base,Direction)
%% Return the distance of Point to the line defined by Base,Direction
% Vectors are defined as [Coord, NSamp] (Dim=1) or [NSamp, Coord]

% Check dimensions
S1=size(Point);
S2=size(Base);
S3=size(Direction);
if any(S1-S2) || any(S1-S3), error('vectors should have same size'); end
if S1(1)==1 || S1(1)>S1(2)
    Dim = 2;
else
    Dim = 1;
end

% Line defined by Base and Point2
Point2 = Base+Direction;
% The equation for getting the distance
Numerator = cross((Point-Base), (Point-Point2));
NumeratorLen = sqrt(sum(Numerator.^2,Dim));
DenominatorLen = sqrt(sum((Point2-Base).^2,Dim));
Distance = NumeratorLen ./ DenominatorLen;
end
