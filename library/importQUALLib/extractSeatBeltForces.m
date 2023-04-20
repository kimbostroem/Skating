function Stream = extractSeatBeltForces(Stream, Data, DefBody, Segments, Info,DoInsert_IsContact)
%% get seat and belt forces from the QTM data structure
% The force2 structure is not QTM compatible, and so not desirable.
% This function is deprecated (use extractErgoForces instead). We still keep it for backwards
% compatibility
%
% SYNTAX
%     Stream = extractSeatBeltForce(Stream, Data, RefBody, Segments, Info,DoInsert_IsContact)
%
% INPUT
%     Stream   (MNData struct) Extracted Data stream to which ergometer forces are added
%     Data     (QTM struct) Imported Data possibly containing ergometer forces
%     DefBody  (Def_Body struct) struct of tables imported from Def_Body.xlsx
%     Segments (struct) Coordinate systems of the segments
%     Info     (struct) containing data information such as time
%     DoInsert_IsContact (logical) 
%            
% OUTPUT
%     Stream  (MNData struct) Imported Data stream to which ergometer forces are added
%
% See also: extractForcePlates, extractErgoForces, importQUAL
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Maarten van den Heuvel, Kim J. Bostroem
% version 230203 (MdL & LK) make names of loop iterators meaningfull
% version 230310 (MdL) handle units

narginchk(5,6)
if nargin<6 || isempty(DoInsert_IsContact)
    DoInsert_IsContact = false;
end

if isfield(Data,'Units')
    if strcmp(Data.Units,'mm')
        Scale = 0.001;
    elseif strcmp(Data.Units,'m')
        Scale = 1;
    else
        error('loadInputFile: non-defined units %s',Data.Units);
    end
else
    Scale = 1;
end
QTMTime      = Info.QTMTime;
QTMTimeSRint = Info.QTMTimeSRint;

if isfield(Data,'Force2')
    ErgoName = {'belt_force'};
    if isfield(Data.Force2,'Name')
        names = {Data.Force2.Name};
    elseif isfield(Data.Force2,'ForcePlateName')
        names = {Data.Force2.ForcePlateName};
    else
        names = {};
    end
    idx = strcmp(names,'seatbelt');
    if any(idx)
        if length(idx) > 1
            id = find(idx);
        else
            id = 1;
        end
        fprintf('Seat belt force data found\n');
        Sz=size(Data.Force2(1).COP);
        % there are different designs of the Force2 structure
        beltCOP = zeros(3,length(QTMTimeSRint));
        beltForce = zeros(3,length(QTMTimeSRint));
        for iCoordinate = 1:3
            if Sz(1) == 3
                beltCOP(iCoordinate,:)   = interp1(QTMTime,Data.Force2.COP(iCoordinate,:),  QTMTimeSRint) * Scale;
                beltForce(iCoordinate,:) = interp1(QTMTime,Data.Force2.Force(iCoordinate,:),QTMTimeSRint);
            else
                beltCOP(iCoordinate,:)   = interp1(QTMTime,Data.Force2(id).COP(:,iCoordinate),  QTMTimeSRint) * Scale;
                beltForce(iCoordinate,:) = interp1(QTMTime,Data.Force2(id).Force(:,iCoordinate),QTMTimeSRint);
            end
        end
        % make relative to thorax
        beltCOP = beltCOP - squeeze(Segments.thorax(4,:,:));
        % remove NaNs
        beltCOP(isnan(beltCOP)) = 0;
        beltForce(isnan(beltForce)) = 0;
        
        if isfield(Stream,'forces')
            st = length(Stream.forces);
        else
            st = 0;
        end
        idx = DefBody.forces.id(matches(DefBody.forces.name,ErgoName));
        Stream.forces(st+1).id = idx;
        Stream.forces(st+1).name = ErgoName;
        Stream.forces(st+1).dynamic.force.units = 'N';
        Stream.forces(st+1).dynamic.force.data = beltForce';
        Stream.forces(st+1).dynamic.COP.units = 'm';
        Stream.forces(st+1).dynamic.COP.data = beltCOP';
        if DoInsert_IsContact
            Stream.forces(st+1).dynamic.isContact.units = 'bool';
            Stream.forces(st+1).dynamic.isContact.data  = true(length(QTMTimeSRint),1);
        end
    end
end
end
