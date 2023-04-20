function JointStream = makeMNJoint(DefBody,Name,Rot,Trans)
%% Creates a substructure from the rotation and translation matrices
% converts euler angles to rotation vector components
%
% SYNTAX
% JointStream = makeMNJoint(DefBody,Name,Rot,Trans)
%
% INPUT
%     DefBody (table structure) contains all sheets from def_body
%     Name    (char) name of the joint 
%     Rot     (double) euler representation of the joint angle
%     Trans   (double) optional vector of the translation 
%
% OUTPUT
%     JointStream  (struct) MN data structure
%
% Local functions: 
%
% See also: calcJointStream
% 
% (c) 2020 by Predimo GmbH
% Author: Marc de Lussanet
% version 230204 (Mdl & LK) Added and improved header

idx = DefBody.joints.id(matches(DefBody.joints.name,Name));
JointStream = struct('id',idx,'name',Name);

% transform the euler angles to rotation vector representation 
rotvec = kbeul2rotvec(Rot*pi/180,'xyz');
JointStream.dynamic.rot.units = 'rad';
JointStream.dynamic.rot.data = rotvec';

% fill the translation components if present
if nargin > 3
    JointStream.dynamic.trans.units = 'm';
    JointStream.dynamic.trans.data = Trans';
end
end
