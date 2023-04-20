function CS = createCS(Origin,Vector1,Vector2,Order)
%% Create a CS (oordinate system) on the basis of two vector time-series and the vector origin.
% vector 1 aligns with the first axis listed in Order. 
% The second coordinate of Order is the norm vector to the two vectors, and
% the third one follows from the first two.
%
% input:
%    origin  [3,ns]  origin/base of the segment
%    vector1 [3,ns]  vector time series aligning with the first axis listed in Origin
%    vector2 [3,ns]  vector time series
%    order   [string] 'xyz', 'xzy', etc order of the vectors
%
% output:
%    CS(1:3,3,1:NS) vectors for axes 1 2 3
%    CS(4  ,3,1:NS) vector to origin of segment
%
% Example:
%     ns     = 100;
%     origin = randi(25,3,ns);
%     vector1= randi(25,3,ns);
%     vector2= randi(25,3,ns);
%     outCS  = cs_create(origin, v1, v2,'yzx');
%
% (c) 2020 Predimo GmbH
% 200529 (MdL) Removed useless parts; implemented for time series of vectors
% 211206 (MdL) minor cleaning, documentation
% 220314 (MdL) deal with data of just one sample

% checks for the dimensionality of the input:
SzO =size(Origin);
SzP1=size(Vector1);
SzP2=size(Vector2);
% if the data is just one sample, it might be transposed
IsOneSample = all(SzO==[1,3]) && all(SzO==SzP1) && all(SzP1==SzP2);
if IsOneSample
    Origin  = Origin(:);
    Vector1 = Vector1(:);
    Vector2 = Vector2(:);
    SzO  = size(Origin);
    SzP1 = size(Vector1);
    SzP2 = size(Vector2);
end
CS = zeros(4,3,SzO(2));
if isempty(Origin) || isempty(Vector1) || isempty(Vector2)
    return; % if markers are not present, the inputs are empty
end
if sum(SzO-SzP1) || sum(SzO-SzP2) || SzO(1)~=3
    error('unexpected dimensionality');
end

% compute three orthogonal unit vectors
Rot_1 = normLength(Vector1);
Rot_2 = normLength(cross(Vector1,Vector2,1));
Rot_3 = normLength(cross(Rot_1,Rot_2,1));

% vector to origin of segment
CS(4,:,:) = Origin;

% rotation matrix
switch Order
    case 'xyz'
        CS(1,:,:) = Rot_1;
        CS(2,:,:) = Rot_2;
        CS(3,:,:) = Rot_3;
    case 'xzy'
        CS(1,:,:) = Rot_1;
        CS(3,:,:) = Rot_2;
        CS(2,:,:) =-Rot_3;
    case 'yxz'
        CS(2,:,:) = Rot_1;
        CS(1,:,:) = Rot_2;
        CS(3,:,:) =-Rot_3;
    case 'yzx'
        CS(2,:,:) = Rot_1;
        CS(3,:,:) = Rot_2;
        CS(1,:,:) = Rot_3;
    case 'zxy'
        CS(3,:,:) = Rot_1;
        CS(1,:,:) = Rot_2;
        CS(2,:,:) = Rot_3;
    case 'zyx'
        CS(3,:,:) = Rot_1;
        CS(2,:,:) = Rot_2;
        CS(1,:,:) =-Rot_3;
    otherwise
        error('something is wrong with the order "%s"',Order);
end
end
