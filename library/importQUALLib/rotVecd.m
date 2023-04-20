function out = rotVecd(v,axis,angdeg)
% rotate vector v about the vector "axis" by an angle angdeg (degrees)
% generalized for time series of vectors
% v and axis may be of size [3,ns], [1,3] or [3,1]
% 
% (c) 2020 by Predimo
% 200529 (MdL)

% check dimensionality
szv=size(v); sza=size(axis);
if sum(szv-sza) || (sza(1)~=3 && sza(2)~=3)
    error('unexpected dimensionality');
end

% get the axes and rotations
axis = normLength(axis);                            % MdL : normLength was normalize
% dot product gives the cos of the angle between the vectors
vn = dot(axis,v).*axis;
vx = v - vn;
vy = cross(axis,vx);
out = vn + vx*cosd(angdeg) + vy*sind(angdeg);
end
