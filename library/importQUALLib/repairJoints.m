function SimInMN = repairJoints(SimInMN,DefBody,Info)
%% Check if any joints are missing or partially missing. 
% Add missing joints with zero values and fill those that are partially missing
%
% Version 230208 (MdL) do not check for joints field (this is done in importQUAL)

%% Init
IsDefaultNPose = true;
Freq = Info.QTMFreqSRint;
Ns = length(SimInMN.signals.time.data);
NJoints = length(SimInMN.signals.joints);

%% add required joints, if necessary
JointNames = {SimInMN.signals.joints.name};
DefNames = DefBody.joints.name;
Required = DefBody.joints.required;
NJointNames = length(DefNames);
Missing = false(NJointNames,1);
for iName = 1:NJointNames
    if Required(iName)
        Missing(iName) = ~any(matches(DefNames(iName),JointNames));
    end
    if Missing(iName)
        % add missing required joints with nan values
        NJoints = NJoints+1;
        SimInMN.signals.joints(NJoints).id   = DefBody.joints.id(iName);
        SimInMN.signals.joints(NJoints).name = DefNames{iName};
        SimInMN.signals.joints(NJoints).dynamic = struct;
        SimInMN.signals.joints(NJoints).dynamic.rot = struct;
        SimInMN.signals.joints(NJoints).dynamic.rot.units = 'rad';
        SimInMN.signals.joints(NJoints).dynamic.rot.data  = nan(Ns,3);
    end
end

%% replace empty (all nan) series with zeros and interpolate gaps
for iJoint = 1:NJoints % loop all joints
    % rotations
    IsGap = all(isnan(SimInMN.signals.joints(iJoint).dynamic.rot.data),2);
    JointIsAllNaN = all(IsGap);
    if JointIsAllNaN
        % set joint angle to zero
        fprintf('Setting missing joint "%s" to zero\n',SimInMN.signals.joints(iJoint).name)
        SimInMN.signals.joints(iJoint).dynamic.rot.data  = zeros(Ns,3);
        if IsDefaultNPose % in absence of arm markers set shoulder not to t-pose but to n-pose
            X=1; Y=2; Z=3; %#ok<NASGU> 
            ShoulderNames = {'right_shoulder_joint','left_shoulder_joint'};
            ShoulderNo = find(contains(JointNames,ShoulderNames));
            if iJoint==ShoulderNo(1) 
                SimInMN.signals.joints(iJoint).dynamic.rot.data(:,X) =-0.45*pi();
            elseif iJoint==ShoulderNo(2)
                SimInMN.signals.joints(iJoint).dynamic.rot.data(:,X) = 0.45*pi();
            end
        end
    elseif any(IsGap)
        warning('repairJoints : Gaps are still present, but any gaps should already have been filled at this point!');
        % linear interpolation
        DataWithGaps = SimInMN.signals.joints(iJoint).dynamic.rot.data;
        [DataFilled,Gaps]  = polyfillgaps(DataWithGaps, Freq, 1,1);
        WereAnyGapsFilled = ~isempty(Gaps);
        SimInMN.signals.joints(iJoint).dynamic.rot.data = DataFilled;
        % print a statement that gaps were filled and include the start and end samples for each gap
        if WereAnyGapsFilled
            fprintf('Filling gaps in joint "%s":\t',SimInMN.signals.joints(iJoint).name);
            % loop the gaps
            for iGap=1:size(Gaps,1)
                fprintf('%d-%d\t',Gaps(iGap,1),Gaps(iGap,2)); 
            end
            fprintf('\n');
        end
    end
    % translations (only if present)
    if isfield(SimInMN.signals.joints(iJoint).dynamic,'trans')
        JointIsAllNaN = all(isnan(SimInMN.signals.joints(iJoint).dynamic.trans.data),'all');
        if JointIsAllNaN % set to default value
            SimInMN.signals.joints(iJoint).dynamic = rmfield(SimInMN.signals.joints(iJoint).dynamic,'trans');
        elseif any(IsGap)
            % interpolate linearly and extrapolate with padding
            DataWithGaps = SimInMN.signals.joints(iJoint).dynamic.trans.data;
            DataFilled = polyfillgaps(DataWithGaps, Freq,1,1);
            % figure; plot(DataFilled);hold on;plot(DataWithGaps);
            SimInMN.signals.joints(iJoint).dynamic.trans.data = DataFilled;
        end
    end
end
end
